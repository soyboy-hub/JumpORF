#!/usr/bin/env python

import os
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import re
from Bio.Seq import Seq
from functools import partial
from concurrent.futures import ProcessPoolExecutor
import uuid
from pathlib import Path

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

def parse_dante_ltr(tmp, outpath, prefix):
    with open(f'{tmp}/{prefix}.dante_ltr.gff3', 'r') as old, \
    open(f'{outpath}/{prefix}.ltr.gff3', 'w') as new, \
    open(f'{tmp}/{prefix}.ltr.partial.gff3', 'w') as new1:
        intact_elements = 0
        count = 0
        dict_ = {}
        for line in old:
            if '#' in line:
                continue
            state = 'intact'
            if 'partial' in line:
                state = 'partial'
            line = line.strip().split('\t')
            if line[2] == 'transposable_element':
                count += 1
                superfamily = line[8].split('Final_Classification=')[1].split(';')[0].split('|')[2]
                lineage = line[8].split('Final_Classification=')[1].split(';')[0].split('|')[-1]
                ID = line[8].split('ID=')[1].split(';')[0]
                ID_ = f'LTR_TE_{count}'
                dict_[ID] = ID_
                new_line = f'{line[0]}\tdante_ltr\ttransposable_element\t{line[3]}\t{line[4]}\t.\t{line[6]}\t.\tID={ID_};Name={ID_};Class=I;Superfamily={superfamily};Lineage={lineage}\n'
                if state == 'intact':
                    intact_elements += 1
                    new.write(new_line)
                if state == 'partial':
                    new1.write(new_line)
            if line[2] == 'long_terminal_repeat':
                parent = line[8].split('Parent=')[1].split(';')[0]
                try:
                    parent_ = dict_[parent]
                except:
                    continue
                ltr_num = line[8].split('LTR=')[1].split(';')[0]
                new_line = f'{line[0]}\tdante_ltr\tlong_terminal_repeat\t{line[3]}\t{line[4]}\t.\t{line[6]}\t.\tID={parent_}:{ltr_num};Name={parent_}:{ltr_num};Parent={parent_}\n'
                if state == 'intact':
                    new.write(new_line)
                if state == 'partial':
                    new1.write(new_line)
            if line[2] == 'protein_domain':
                domain_type = line[8].split('Name=')[1].split(';')[0]
                parent = line[8].split('Parent=')[1].split(';')[0]
                try:
                    parent_ = dict_[parent]
                except:
                    continue
                new_line = f'{line[0]}\tdante_ltr\tprotein_domain\t{line[3]}\t{line[4]}\t.\t{line[6]}\t.\tID={parent_}:{domain_type};Name={parent_}:{domain_type};Parent={parent_}\n'
                if state == 'intact':
                    new.write(new_line)
                if state == 'partial':
                    new1.write(new_line)
                    
    return intact_elements

def cluster_bed_by_class_self(tmp, prefix, domain_name, max_distance=300):

    tmp_file = f'{tmp}/{prefix}_{domain_name}_tmp.dom.bed' + ".clustered"

    # run bedtools clustering (strand-aware)
    os.system(f"bedtools cluster -i {tmp}/{prefix}_{domain_name}_tmp.dom.bed -d {max_distance} -s > {tmp_file}")

    clusters = defaultdict(list)

    with open(tmp_file) as f:

        for line in f:
            line = line.strip()
            if not line:
                continue

            cols = line.split()

            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])

            name = cols[3]
            strand = cols[5]

            te_class = name.split("::")[1] if "::" in name else "NA"

            cluster_id = cols[-1]

            key = f"{chrom}|{strand}|{te_class}|{cluster_id}"

            clusters[key].append((start, end))

    with open(f'{tmp}/{prefix}_{domain_name}.dom.bed', "w") as out:

        te_counter = 1

        for key, intervals in clusters.items():

            chrom, strand, te_class, cid = key.split("|")

            merged_start = min(i[0] for i in intervals)
            merged_end = max(i[1] for i in intervals)

            count = len(intervals)

            new_name = f"{domain_name}_{te_counter}::{te_class}"

            out.write(
                f"{chrom}\t{merged_start}\t{merged_end}\t{new_name}\t0\t{strand}\n"
            )

            te_counter += 1

    #os.system(f"rm {tmp_file}")

def mask_annotated_LTR_RTE_free_INT(genome, tmp, outpath, prefix):
    os.system(f'bedtools maskfasta -fi {genome} -bed {outpath}/{prefix}.ltr.gff3 -fo {tmp}/{prefix}.masked.LTR_RTE_intact.fasta')
    os.system(f'bedtools maskfasta -fi {tmp}/{prefix}.masked.LTR_RTE_intact.fasta -bed {tmp}/{prefix}.ltr.partial.gff3 -fo {tmp}/{prefix}.masked.LTR_RTE_intact_partial.fasta')
    os.system(f'grep "INT" {tmp}/{prefix}.domain_dante.gff3 > {tmp}/{prefix}.domain_dante.INT.gff3')
    os.system(f'bedtools maskfasta -fi {tmp}/{prefix}.masked.LTR_RTE_intact_partial.fasta -bed {tmp}/{prefix}.domain_dante.INT.gff3 -fo {tmp}/{prefix}.masked.LTR_RTE.masked.fasta')
    
    return f'{tmp}/{prefix}.masked.LTR_RTE.masked.fasta'

def mask_annotated_DNA_TE(genome_masked_ltr, tmp, outpath, prefix):
    os.system(f'bedtools maskfasta -fi {genome_masked_ltr} -bed {outpath}/{prefix}.dna.gff3 -fo {tmp}/{prefix}.masked.LTR_RTE.DNA.fasta')
    
    return f'{tmp}/{prefix}.masked.LTR_RTE.DNA.fasta'

def mask_annotated_Helitron_TE(genome_masked_ltr_dna, tmp, outpath, prefix):
    os.system(f'bedtools maskfasta -fi {genome_masked_ltr_dna} -bed {outpath}/{prefix}.helitron.gff3 -fo {tmp}/{prefix}.masked.LTR_RTE.DNA.HEL.fasta')
    
    return f'{tmp}/{prefix}.masked.LTR_RTE.DNA.HEL.fasta'
    
def annotate_domain(domain_hmm, domain_name, genome, processors, tmp, prefix):
    if domain_hmm == 'rexdb-pnas':
        os.system(f'TEsorter -db rexdb-pnas -p {processors} -eval 0.0001 -cov 10 -genome -pre {tmp}/{prefix}_{domain_name}_tmp {genome} > /dev/null 2>&1')
    if domain_hmm != 'rexdb-pnas':
        os.system(f'TEsorter --db-hmm {domain_hmm} -p {processors} -eval 0.0001 -cov 10 -genome -pre {tmp}/{prefix}_{domain_name}_tmp {genome} > /dev/null 2>&1')
    os.system(f'rm {tmp}/{prefix}_{domain_name}_tmp.domtbl')

def cluster_domain_nucl(domain_name, genome, tmp, prefix, processors, min_cov, min_seq_id):
    with open(f'{tmp}/{prefix}_{domain_name}_tmp.dom.gff3', 'r') as old, \
    open(f'{tmp}/{prefix}_{domain_name}_tmp.dom.bed', 'w') as new:
        c = 0
        for line in old:
            c += 1
            line = line.split('\t')
            id_ = f'{domain_name}_{c}'
            classification = line[8].split('Classification=')[1].split(';')[0]
            new_line = f'{line[0]}\t{line[3]}\t{line[4]}\t{id_}::{classification}\t0\t{line[6]}\n'
            new.write(new_line)
    cluster_bed_by_class_self(tmp, prefix, domain_name)
    os.system(f'bedtools getfasta -fi {genome} -bed {tmp}/{prefix}_{domain_name}.dom.bed -fo {tmp}/{prefix}_{domain_name}.dom.fasta -s -nameOnly')
    os.system(f'mmseqs easy-cluster {tmp}/{prefix}_{domain_name}.dom.fasta {tmp}/{prefix}_{domain_name}.dom.cluster {tmp}/tmp_mmseqs_{domain_name} --max-seqs 1000 -c {min_cov} --cov-mode 0 --min-seq-id {min_seq_id} --threads {processors} > /dev/null 2>&1')

def parse_annotation_clustering(domain_name, tmp, prefix):
    with open(f'{tmp}/{prefix}_{domain_name}.dom.bed', 'r') as old, \
    open(f'{tmp}/{prefix}_{domain_name}.dom.cluster_cluster.tsv', 'r') as old1, \
    open(f'{tmp}/{prefix}_{domain_name}.dom.clustered.bed', 'w') as new:
        dict_clusters = {}
        dict_clusters_refined = {}
        for line in old1:
            line = line.strip().split('\t')
            line[0] = line[0].split('(')[0]
            line[1] = line[1].split('(')[0]
            cluster_id = line[0]
            cluster_member = line[1]
            if cluster_id not in dict_clusters:
                dict_clusters[cluster_id] = []
            dict_clusters[cluster_id].append(cluster_member)
        c = 0
        for cluster in dict_clusters:
            if len(dict_clusters[cluster]) < 2:
                continue
            c += 1
            dict_clusters_refined[f'cluster_{c}'] = dict_clusters[cluster]
        c = 0
        for line in old:
            line = line.strip().split('\t')
            id_ = line[3].split('::')[0]
            classification = line[3].split('::')[1]
            for cluster in dict_clusters_refined:
                if line[3] in dict_clusters_refined[cluster]:
                    c += 1
                    new_id = f'{domain_name}_{c}'
                    new_line = f'{line[0]}\t{line[1]}\t{line[2]}\t{new_id}::{cluster}::{classification}\t0\t{line[5]}\n'
                    new.write(new_line)

def write_fasta(records, filepath):
    with open(filepath, 'w') as f:
        for seq_id, seq in records:
            f.write(f'>{seq_id}\n{seq}\n')

def chunk_records(records, chunk_size=50):
    for i in range(0, len(records), chunk_size):
        yield records[i:i + chunk_size]

def merge_boundary_dicts(dicts):
    merged = {}
    for d in dicts:
        merged.update(d)
    return merged

def merge_intervals(data, gap=200):
    result = {}

    for key, intervals in data.items():
        intervals = sorted(intervals, key=lambda x: x[0])

        merged = []
        current_start, current_end = intervals[0]

        for s2, e2 in intervals[1:]:
            s1, e1 = current_start, current_end

            if s2 <= e1:
                current_end = max(e1, e2)
            elif (s2 - e1) < gap:
                current_end = e2
            else:
                merged.append([current_start, current_end])
                current_start, current_end = s2, e2

        merged.append([current_start, current_end])
        result[key] = merged

    return result


def filter_intervals(d1, d2):
    result = {}

    for key, intervals in d1.items():
        if key not in d2:
            continue

        s2, e2 = d2[key]

        filtered = []
        for s1, e1 in intervals:
            if s1 <= s2 and e1 >= e2:
                filtered.append(s2 - s1)
                filtered.append(e1 - e2)

        if filtered:
            result[key] = filtered

    return result

def run_blastn_merge_boundaries(fasta, reg_to_adj, processors):

    domain_position = {}
    for seq in SeqIO.parse(fasta, 'fasta'):
        seq_len = len(seq.seq)
        domain_position[seq.id] = [reg_to_adj, seq_len - reg_to_adj]

    os.system(f"makeblastdb -dbtype nucl -in {fasta} > /dev/null 2>&1")
    os.system(f"blastn -query {fasta} -db {fasta} "
              f"-task megablast "
              f"-out {fasta}.m8 -evalue 1e-20 "
              f"-outfmt '6 qseqid sseqid qstart qend' "
              f"-dust no -soft_masking false -perc_identity 90 "
              f"-num_threads {processors}")

    hsps = {}
    with open(f"{fasta}.m8") as f:
        for line in f:
            q, s, start, end = line.strip().split("\t")
            if q == s:
                continue
            hsps.setdefault(q, []).append([int(start), int(end)])

    merged = merge_intervals(hsps)
    return filter_intervals(merged, domain_position)

def adjast_align_clusters(genome, tmp, prefix, domain_name, ref_to_adj, processors):
    dict_genome_lens = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        dict_genome_lens[seq.id] = len(seq.seq)
        
    clusters = {}
    
    with open(f'{tmp}/{prefix}_{domain_name}.dom.clustered.bed', 'r') as old, \
    open(f'{tmp}/{prefix}_{domain_name}.dom.clustered.adj.bed', 'w') as new:
        for line in old:
            line = line.split('\t')
            st = int(line[1])
            end = int(line[2])
            new_st = st - ref_to_adj
            new_end = end + ref_to_adj
            if new_st <= 0:
                new_st = 1
            if new_end > dict_genome_lens[line[0]]:
                new_end = dict_genome_lens[line[0]]
                
            new_line = f'{line[0]}\t{new_st}\t{new_end}\t{line[3]}\t0\t{line[5]}'
            new.write(new_line)
            
            cluster = line[3].split('::')[1]
            
            if cluster not in clusters:
                clusters[cluster] = []
            clusters[cluster].append(line[3])
    
    os.system(f'bedtools getfasta -fi {genome} -bed {tmp}/{prefix}_{domain_name}.dom.clustered.adj.bed -fo {tmp}/{prefix}_{domain_name}.dom.clustered.adj.fasta -s -nameOnly')
    
    dict_seqs = {}
    
    for seq in SeqIO.parse(f'{tmp}/{prefix}_{domain_name}.dom.clustered.adj.fasta', 'fasta'):
        new_seq_id = seq.id.split('(')[0]
        dict_seqs[new_seq_id] = seq.seq
        
    element_boundaries = {}

    for cluster in clusters:
        seq_records = []

        # Collect sequences as (id, seq)
        for seq_id in clusters[cluster]:
            seq_seq = dict_seqs[seq_id]
            seq_records.append((seq_id, seq_seq))

        all_results = []
    
        if len(seq_records) <= 100:
            fasta_file = f'{tmp}/{prefix}_{domain_name}.dom.clustered.adj.{cluster}.fasta'

            write_fasta(seq_records, fasta_file)
            
            res = run_blastn_merge_boundaries(fasta_file, ref_to_adj, processors)
            all_results.append(res)

            os.system(f'rm {fasta_file}')
            os.system(f'rm {fasta_file}*')

        else:
            for i, chunk in enumerate(chunk_records(seq_records, 50)):
                chunk_file = f'{tmp}/{prefix}_{domain_name}.dom.clustered.adj.{cluster}.chunk{i}.fasta'
                
                write_fasta(chunk, chunk_file)
                                
                res = run_blastn_merge_boundaries(chunk_file, ref_to_adj, processors)
                all_results.append(res)

                os.system(f'rm {chunk_file}')
                os.system(f'rm {chunk_file}*')
    
        # Merge results
        merged = merge_boundary_dicts(all_results)
        element_boundaries.update(merged)
        
    with open(f'{tmp}/{prefix}_{domain_name}.dom.clustered.bed', 'r') as old, \
    open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'w') as new:
        for line in old:
            line = line.strip().split('\t')
            id_ = line[3]
            if id_ in element_boundaries:
                adj_left = element_boundaries[id_][0]
                adj_right = element_boundaries[id_][1]
                if line[5] == '+':
                    new_st = int(line[1]) - int(adj_left)
                    if new_st <= 0:
                        new_st = 1
                    new_end = int(line[2]) + int(adj_right)
                    if new_end > dict_genome_lens[line[0]]:
                        new_end = dict_genome_lens[line[0]]
                    new_line = f'{line[0]}\t{new_st}\t{new_end}\t{id_}\t0\t{line[5]}\n'
                if line[5] == '-':
                    new_st = int(line[1]) - int(adj_right)
                    if new_st <= 0:
                        new_st = 1
                    new_end = int(line[2]) + int(adj_left)
                    if new_end > dict_genome_lens[line[0]]:
                        new_end = dict_genome_lens[line[0]]
                    new_line = f'{line[0]}\t{new_st}\t{new_end}\t{id_}\t0\t{line[5]}\n'
                new.write(new_line)

def recluster_nested(tmp, prefix, domain_name, cluster, processors):
    os.system(f'mmseqs easy-cluster {tmp}/{prefix}_{domain_name}.element.nested.{cluster}.fasta {tmp}/{prefix}_{domain_name}.element.nested.{cluster}.fasta.cluster {tmp}/tmp_mmseqs_{domain_name} --cov-mode 5 --min-seq-id 0.8 -c 0.95 --threads {processors} > /dev/null 2>&1')
    
    dict_seqs = {}
    for seq in SeqIO.parse(f'{tmp}/{prefix}_{domain_name}.element.nested.{cluster}.fasta', 'fasta'):
        seq.id = seq.id.split('(')[0]
        dict_seqs[seq.id] = seq.seq
    
    cluster_reassign_ = {}
    to_save = []
    
    with open(f'{tmp}/{prefix}_{domain_name}.element.nested.{cluster}.fasta.cluster_cluster.tsv', 'r') as old:
        for line in old:
            line = line.strip().split('\t')
            line[0] = line[0].split('(')[0]
            line[1] = line[1].split('(')[0]
            if line[0] not in cluster_reassign_:
                cluster_reassign_[line[0]] = []
            cluster_reassign_[line[0]].append(line[1])
        
        for cluster in cluster_reassign_:
            if len(cluster_reassign_[cluster]) < 2:
                continue
            cluster_representative_seq = dict_seqs[cluster]
            
            first_100 = cluster_representative_seq[:100]
            last_100 = Seq(str(cluster_representative_seq[-100:])).reverse_complement()
            
            result_tir = local_align(first_100, last_100, 1, -1, -2)
            if result_tir < 5:
                continue
            if 'MULE' in cluster and result_tir < 30:
                continue
            
            for element in cluster_reassign_[cluster]:
                to_save.append(element)
    return to_save
          
    
def local_align(s1, s2, match, mismatch, gap):
    """Simple Smith-Waterman local alignment (returns length too)"""
    n, m = len(s1), len(s2)
    dp = [[0]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_pos = (0, 0)

    # Fill DP matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            score = match if s1[i-1] == s2[j-1] else mismatch

            dp[i][j] = max(
                0,
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_pos = (i, j)

    # Traceback
    i, j = max_pos
    aln1, aln2 = [], []
    aln_len = 0  # <-- track alignment length

    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch):
            aln1.append(s1[i-1])
            aln2.append(s2[j-1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i-1][j] + gap:
            aln1.append(s1[i-1])
            aln2.append('-')
            i -= 1
        else:
            aln1.append('-')
            aln2.append(s2[j-1])
            j -= 1

        aln_len += 1  # count every aligned position (including gaps)

    if aln_len < 5:
        return 0
    else:
        return aln_len
    
def recluster_get_representative(genome, tmp, prefix, domain_name, processors):
    
    elements_passed = []
    
    os.system(f'bedtools getfasta -fi {genome} \
    -bed {tmp}/{prefix}_{domain_name}.element.nested.bed \
    -fo {tmp}/{prefix}_{domain_name}.element.nested.fasta -s -nameOnly')

    clusters = {}
    for seq in SeqIO.parse(f'{tmp}/{prefix}_{domain_name}.element.nested.fasta', 'fasta'):
        cluster = seq.id.split('::')[1]
        if cluster not in clusters:
            clusters[cluster] = {}
        clusters[cluster][seq.id] = seq.seq

    for cluster in clusters:
        with open(f'{tmp}/{prefix}_{domain_name}.element.nested.{cluster}.fasta', 'w') as new:
            for seq in clusters[cluster]:
                seq_seq = clusters[cluster][seq]
                new_line = f'>{seq}\n{seq_seq}\n'
                new.write(new_line)
        
        to_save = recluster_nested(tmp, prefix, domain_name, cluster, processors)
        elements_passed += to_save
    
    return elements_passed

def generate_final_annot_dna(genome, outpath, tmp, prefix, domain_name, te_name, name_to_report, class_, processors):
    
    elements_passed = recluster_get_representative(genome, tmp, prefix, domain_name, processors)
    c = 0
    count_to_domain = {}
    
    with open(f'{tmp}/{prefix}_{domain_name}.dom.clustered.bed', 'r') as domain, \
    open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'r') as elements, \
    open(f'{outpath}/{prefix}.{name_to_report}.gff3', 'w') as new:
        for line in elements:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in elements_passed:
                c += 1
                count_to_domain[tmp_id] = c
                new_id = f'{te_name}_{c}'
                alias = tmp_id.split('::')[2]
                new_line = f'{line[0]}\tJumpORF\ttransposable_element\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Family={alias}\n'
                new.write(new_line)
        for line in domain:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in count_to_domain:
                c = count_to_domain[tmp_id]
                new_id = f'{te_name}_{c}::{domain_name}'
                parent_id = f'{te_name}_{c}'
                alias = tmp_id.split('::')[2]
                new_line = f'{line[0]}\tJumpORF\tprotein_domain\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Family={alias};Parent={parent_id}\n'
                new.write(new_line)
    print(f'INFO: intact DNA-TEs found: {c}')
            
def find_dual_domain_elements(tmp, prefix, domain_name_element, domain_name_second):
    dual_elements = []
    dual_second_pair = {}
    os.system(f'bedtools intersect -wo -a {tmp}/{prefix}_{domain_name_element}.element.nested.bed -b {tmp}/{prefix}_{domain_name_second}.dom.clustered.bed > {tmp}/{prefix}_{domain_name_element}.{domain_name_second}.element.nested.bed')
    with open(f'{tmp}/{prefix}_{domain_name_element}.{domain_name_second}.element.nested.bed', 'r') as old:
        for line in old:
            line = line.strip().split('\t')
            element = line[3]
            dual_elements.append(element)
            dual_second_pair[line[9]] = element
            
    return dual_elements, dual_second_pair

def find_mult_domain_elements(tmp, prefix, domain_name_element, list_domain_second):
    mult_elements = []
    mult_second_pair = {}
    mult_second_pair_ = {}
    mult_second_pair_tmp = {}
    
    for domain in list_domain_second:
        os.system(f'bedtools intersect -wo -a {tmp}/{prefix}_{domain_name_element}.element.nested.bed -b {tmp}/{prefix}_{domain}.dom.clustered.bed > {tmp}/{prefix}_{domain_name_element}.{domain}.element.nested.bed')
        with open(f'{tmp}/{prefix}_{domain_name_element}.{domain}.element.nested.bed', 'r') as old:
            for line in old:
                line = line.strip().split('\t')
                element = line[3]
                if element not in mult_second_pair_tmp:
                    mult_second_pair_tmp[element] = []
                if domain not in mult_second_pair_tmp[element]:
                    mult_second_pair_tmp[element].append(domain)
                mult_second_pair_[line[9]] = element
                
    for element in mult_second_pair_tmp:
        if len(mult_second_pair_tmp[element]) == len(list_domain_second):
            mult_elements.append(element)
            for second in mult_second_pair_:
                if mult_second_pair_[second] == element:
                    mult_second_pair[second] = element
            
    return mult_elements, mult_second_pair

def generate_final_annot_dual(outpath, tmp, prefix, domain_name_element, domain_name_second, te_name, name_to_report, class_, processors):
    dual_elements, dual_second_pair = find_dual_domain_elements(tmp, prefix, domain_name_element, domain_name_second)
    count_to_domain = {}
    c = 0
    with open(f'{tmp}/{prefix}_{domain_name_element}.element.nested.bed', 'r') as element, \
    open(f'{tmp}/{prefix}_{domain_name_element}.dom.clustered.bed', 'r') as first, \
    open(f'{tmp}/{prefix}_{domain_name_second}.dom.clustered.bed', 'r') as second, \
    open(f'{outpath}/{prefix}.{name_to_report}.gff3', 'w') as new:
        for line in element:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in dual_elements:
                c += 1
                count_to_domain[tmp_id] = c
                new_id = f'{te_name}_{c}'
                new_line = f'{line[0]}\tJumpORF\ttransposable_element\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_}\n'
                new.write(new_line)
        for line in first:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in count_to_domain:
                c = count_to_domain[tmp_id]
                new_id = f'{te_name}_{c}::{domain_name_element}'
                parent_id = f'{te_name}_{c}'
                new_line = f'{line[0]}\tJumpORF\tprotein_domain\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Parent={parent_id}\n'
                new.write(new_line)
        for line in second:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in dual_second_pair:
                parent = dual_second_pair[tmp_id]
                c = count_to_domain[parent]
                new_id = f'{te_name}_{c}::{domain_name_second}'
                parent_id = f'{te_name}_{c}'
                new_line = f'{line[0]}\tJumpORF\tprotein_domain\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Parent={parent_id}\n'
                new.write(new_line)
    
    return c

def generate_final_annot_multiple(outpath, tmp, prefix, domain_name_element, list_domain_second, te_name, name_to_report, class_, processors):

    mult_elements, mult_second_pair = find_mult_domain_elements(tmp, prefix, domain_name_element, list_domain_second)
    
    count_to_domain = {}
    c = 0
    with open(f'{tmp}/{prefix}_{domain_name_element}.element.nested.bed', 'r') as element, \
    open(f'{tmp}/{prefix}_{domain_name_element}.dom.clustered.bed', 'r') as first, \
    open(f'{outpath}/{prefix}.{name_to_report}.gff3', 'w') as new:
        for line in element:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in mult_elements:
                c += 1
                count_to_domain[tmp_id] = c
                new_id = f'{te_name}_{c}'
                new_line = f'{line[0]}\tJumpORF\ttransposable_element\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_}\n'
                new.write(new_line)
        for line in first:
            line = line.strip().split('\t')
            tmp_id = line[3]
            if tmp_id in count_to_domain:
                c = count_to_domain[tmp_id]
                new_id = f'{te_name}_{c}::{domain_name_element}'
                parent_id = f'{te_name}_{c}'
                new_line = f'{line[0]}\tJumpORF\tprotein_domain\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Parent={parent_id}\n'
                new.write(new_line)
        for domain in list_domain_second:
            with open(f'{tmp}/{prefix}_{domain}.dom.clustered.bed', 'r') as second:
                for line in second:
                    line = line.strip().split('\t')
                    tmp_id = line[3]
                    if tmp_id in mult_second_pair:
                        parent = mult_second_pair[tmp_id]
                        c = count_to_domain[parent]
                        new_id = f'{te_name}_{c}::{domain}'
                        parent_id = f'{te_name}_{c}'
                        new_line = f'{line[0]}\tJumpORF\tprotein_domain\t{line[1]}\t{line[2]}\t.\t{line[5]}\t.\tID={new_id};Name={new_id};Class={class_};Parent={parent_id}\n'
                        new.write(new_line)
    
    return c

def is_file_non_empty(path):
    import os
    return os.path.getsize(path) > 0

def check_if_annotation(list_annotated, outpath, prefix):
    list_have_annotation = []
    for te_type in list_annotated:
        state = is_file_non_empty(f'{outpath}/{prefix}.{te_type}.gff3')
        if state == True:
            list_have_annotation.append(te_type)
            
    return list_have_annotation

def extract_TEs_and_basal_proteins(list_have_annotation, genome, outpath, tmp, prefix):
    for te_type in list_have_annotation:
        os.system(f'grep "transposable_element" {outpath}/{prefix}.{te_type}.gff3 > {tmp}/{prefix}.{te_type}.transposable_element.gff3')
        os.system(f'grep "protein_domain" {outpath}/{prefix}.{te_type}.gff3 > {tmp}/{prefix}.{te_type}.protein_domain.gff3')
    
    with open(f'{tmp}/{prefix}.all.transposable_element.bed', 'w') as te, \
    open(f'{tmp}/{prefix}.all.protein_domain.gff3', 'w') as dom:
        for te_type in list_have_annotation:
            with open(f'{tmp}/{prefix}.{te_type}.transposable_element.gff3', 'r') as old:
                for line in old:
                    line = line.split('\t')
                    id_ = line[8].split('ID=')[1].split(';')[0]
                    new_line = f'{line[0]}\t{line[3]}\t{line[4]}\t{id_}\t0\t{line[6]}\n'
                    te.write(new_line)
            with open(f'{tmp}/{prefix}.{te_type}.protein_domain.gff3', 'r') as old1:
                for line in old1:
                    dom.write(line)
    os.system(f'bedtools getfasta -name -fi {genome} -bed {tmp}/{prefix}.all.transposable_element.bed -fo {tmp}/{prefix}.all.transposable_element.fasta')

def run_augustus(fasta_file_prefix, gene_model):
    os.system('augustus --species={1} --softmasking=0 --genemodel=complete --singlestrand=true --alternatives-from-evidence=false --alternatives-from-sampling=false --noInFrameStop=true --minexonintronprob=0.8 --minmeanexonintronprob=0.8 --gff3=on --UTR=off --AUGUSTUS_CONFIG_PATH=/module/augustus_config {0}.fasta > {0}.predicted.aORF.gff3'.format(fasta_file_prefix, gene_model))
    os.system('rm {0}.fasta'.format(fasta_file_prefix))
    
def process_batch(batch, gene_model, output_dir):
    prefix = os.path.join(output_dir, f"batch_{uuid.uuid4().hex}")
    fasta_file = f"{prefix}.fasta"
    SeqIO.write(batch, fasta_file, "fasta")
    run_augustus(prefix, gene_model)
    
def run_augustus_on_batches(tmp, prefix, batch_size, gene_model, processors):
    
    processors = int(processors)
    
    outdir = Path(f'{tmp}/{prefix}_augustus')
    outdir.mkdir(parents=True, exist_ok=True)

    sequences = list(SeqIO.parse(f'{tmp}/{prefix}.all.transposable_element.fasta', "fasta"))
    batches = [sequences[i:i + batch_size] for i in range(0, len(sequences), batch_size)]

    wrapped_batch_runner = partial(process_batch, gene_model=gene_model, output_dir=outdir)

    with ProcessPoolExecutor(max_workers=processors) as executor:
        futures = [executor.submit(wrapped_batch_runner, batch) for batch in batches]
        for future in futures:
            future.result()

def process_augustus_out(tmp, prefix):
    files = [f'{tmp}/{prefix}_augustus/{f}' for f in os.listdir(f'{tmp}/{prefix}_augustus')]

    with open(f'{tmp}/{prefix}.aORF.nested.gff3', 'w') as new:
        dict_augustus_prediction_on_batches = {}
        dict_augustus_prediction_on_batches_exons = {}
        for file in files:
            with open(file, 'r') as old:
                for line in old:
                    if '#' in line:
                        continue
                    line = line.strip().split('\t')

                    if line[2] == 'start_codon' or line[2] == 'intron' or line[2] == 'stop_codon':
                        continue
                    if line[2] == 'CDS':
                        line[2] = 'exon'
                    
                    parent = line[0].split('::')[0]
                    chr = line[0].split('::')[1].split(':')[0]
                    start = line[0].split('::')[1].split(':')[1].split('-')[0]
                    start_adj = line[3]
                    end_adj = line[4]
                    new_start = int(start) + int(start_adj)
                    new_end = int(start) + int(end_adj)
                    type_ = line[2]

                    if parent not in dict_augustus_prediction_on_batches:
                        dict_augustus_prediction_on_batches[parent] = 0
                    if parent not in dict_augustus_prediction_on_batches_exons:
                        dict_augustus_prediction_on_batches_exons[parent] = 0
                    
                    if type_ == 'gene':
                        dict_augustus_prediction_on_batches[parent] += 1
                        orf_num = dict_augustus_prediction_on_batches[parent]
                        new_line = f'{chr}\tJumpORF\tgene\t{new_start}\t{new_end}\t.\t{line[6]}\t.\tID={parent}:aORF_{orf_num};Name={parent}:aORF_{orf_num}\n'
                        new.write(new_line)
                    if type_ == 'exon':
                        dict_augustus_prediction_on_batches_exons[parent] += 1
                        orf_num = dict_augustus_prediction_on_batches[parent]
                        exon_number = dict_augustus_prediction_on_batches_exons[parent]
                        new_line = f'{chr}\tJumpORF\texon\t{new_start}\t{new_end}\t.\t{line[6]}\t.\tID={parent}:aORF_{orf_num}:exon_{exon_number};Name={parent}:aORF_{orf_num}:exon_{exon_number};Parent={parent}:aORF_{orf_num}\n'
                        new.write(new_line)

def remove_basal_protein_from_pred(outpath, tmp, prefix):
    os.system(f'grep "JumpORF	gene" {tmp}/{prefix}.aORF.nested.gff3 > {tmp}/{prefix}.aORF.nested.tmp.gff3')
    os.system(f'bedtools intersect -v -a {tmp}/{prefix}.aORF.nested.tmp.gff3 -b {tmp}/{prefix}.all.protein_domain.gff3 > {tmp}/{prefix}.aORF.tmp.gff3')
    with open(f'{tmp}/{prefix}.aORF.tmp.gff3', 'r') as old, \
    open(f'{tmp}/{prefix}.aORF.nested.gff3', 'r') as old1, \
    open(f'{outpath}/{prefix}.aORF.gff3', 'w') as new:
        dict_saved_exons = {}
        ids = []
        lines = {}
        for line in old1:
            line = line.strip().split('\t')
            type_ = line[2]
            if type_ == 'gene':
                name = line[8].split('ID=')[1].split(';')[0]
                if name not in dict_saved_exons:
                    dict_saved_exons[name] = []
                lines[name] = '\t'.join(line)
            if type_ == "exon":
                name = line[8].split('ID=')[1].split(';')[0]
                parent = line[8].split('Parent=')[1]
                dict_saved_exons[parent].append(name)
                lines[name] = '\t'.join(line)
        for line in old:
            line = line.split('\t')
            id = line[8].split('ID=')[1].split(';')[0]
            ids.append(id)
        for gene in dict_saved_exons:
            if len(dict_saved_exons[gene]) != 0 and gene in ids:
                new_line = f'{lines[gene]}\n'
                new.write(new_line)
                for i in dict_saved_exons[gene]:
                    new_line = f'{lines[i]}\n'
                    new.write(new_line)

def extract_transcripts_from_gff3(outpath, tmp, prefix, genome_):
    genome = SeqIO.to_dict(SeqIO.parse(genome_, "fasta"))
    genes = {}
    exons_by_gene = defaultdict(list)

    with open(f'{outpath}/{prefix}.aORF.gff3', 'r') as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            attr_dict = dict(item.strip().split("=") for item in attributes.strip().split(";") if "=" in item)
            if feature_type == "gene":
                gene_id = attr_dict.get("ID")
                if gene_id:
                    genes[gene_id] = {"seqid": seqid, "strand": strand}
            elif feature_type == "exon":
                parent_id = attr_dict.get("Parent")
                if parent_id:
                    exons_by_gene[parent_id].append((start, end))

    with open(f'{tmp}/{prefix}.aORF.transcripts.fasta', "w") as out_f:
        for gene_id, gene_info in genes.items():
            seqid = gene_info["seqid"]
            strand = gene_info["strand"]
            if gene_id not in exons_by_gene:
                continue
            exon_coords = sorted(exons_by_gene[gene_id], key=lambda x: x[0])
            exon_seqs = []

            for start, end in exon_coords:
                exon_seq = genome[seqid].seq[start - 1:end]
                exon_seqs.append(exon_seq)

            transcript_seq = Seq("").join(exon_seqs)
            if strand == "-":
                transcript_seq = transcript_seq.reverse_complement()

            out_f.write(f'>{gene_id}\n{transcript_seq}\n')
        
def run_trans_decoder(tmp, prefix):
    os.system(f'/module/TransDecoder/TransDecoder.LongOrfs -t {tmp}/{prefix}.aORF.transcripts.fasta -m 100 -S --complete_orfs_only --output_dir {tmp} > /dev/null 2>&1')
    with open(f'{tmp}/{prefix}.aORF.transcripts.tmp.faa', 'w') as new:
        for seq in SeqIO.parse(f'{tmp}/{prefix}.aORF.transcripts.fasta.transdecoder_dir/longest_orfs.pep', 'fasta'):
            if '.p1' in seq.id:
                new_seq = seq.seq.split('*')[0]
                new_line = f'>{prefix}::{seq.id}\n{new_seq}\n'
                new.write(new_line)

def generate_final_annotation(outpath, tmp, prefix, processors, min_identity_to_cluster, min_size_cluster, min_coverage):
    
    min_size_cluster = int(min_size_cluster)

    os.system(f'mmseqs easy-cluster -c {min_coverage} --cov-mode 1 --min-seq-id {min_identity_to_cluster} --cluster-mode 0 --threads {processors} {tmp}/{prefix}.aORF.transcripts.tmp.faa {tmp}/{prefix}.aORF.transcripts.tmp.faa {tmp}/tmp_mmseqs_aORF > /dev/null 2>&1')

    with open(f'{tmp}/{prefix}.aORF.transcripts.tmp.faa_cluster.tsv', 'r') as old, \
    open(f'{outpath}/{prefix}.aORF.gff3', 'r') as old1, \
    open(f'{outpath}/{prefix}.aORF.clustered.gff3', 'w') as new, \
    open(f'{outpath}/{prefix}.aORF.clusters', 'w') as new1, \
    open(f'{outpath}/{prefix}.aORF.clustered.faa', 'w') as new2:

        clusters = {}
        orig_seqs = {}
        clusters_final = {}

        for seq in SeqIO.parse(f'{tmp}/{prefix}.aORF.transcripts.tmp.faa', 'fasta'):
            seq_id = seq.id.split('::')[1].split('.p')[0]
            orig_seqs[seq_id] = seq.seq

        for line in old:
            line = line.strip().split('\t')
            if line[0] not in clusters:
                clusters[line[0]] = []
            member = line[1].split('::')[1].split('.p')[0]
            clusters[line[0]].append(member)
        
        count = 0
        for cluster in clusters:
            if len(clusters[cluster]) < min_size_cluster:
                continue
            count += 1
            clusters_final[f'cluster_{count}'] = clusters[cluster]
            for memeber in clusters[cluster]:
                new_line = f'cluster_{count}\t{memeber}\n'
                new1.write(new_line)

                member_seq = orig_seqs[memeber]
                new_line_1 = f'>{prefix}::{memeber}::cluster_{count}\n{member_seq}\n'
                new2.write(new_line_1)
        
        for line in old1:
            line1 = line.split('\t')
            if line1[2] == 'gene':
                id = line1[8].split('ID=')[1].split(';')[0]
                for cluster in clusters_final:
                    if id in clusters_final[cluster]:
                        new.write(line)
            if line1[2] == 'exon':
                id = line1[8].split('ID=')[1].split(';')[0].split(':exon')[0]
                for cluster in clusters_final:
                    if id in clusters_final[cluster]:
                        new.write(line)

def check_if_LTR_DIRS(fasta, outpath):
    real = False
    with open(f'{outpath}/first_200.fasta', 'w') as new:
        for seq in SeqIO.parse(fasta, 'fasta'):
            first_200 = seq.seq[:200]
            new.write(f'>first_200\n{first_200}\n')
    with open(f'{outpath}/last_100.fasta', 'w') as new:
        for seq in SeqIO.parse(fasta, 'fasta'):
            last_100 = seq.seq[-100:]
            new.write(f'>last_100\n{last_100}\n')
    os.system(f'makeblastdb -dbtype nucl -in {fasta} > /dev/null 2>&1')
    os.system(f'blastn -query {outpath}/first_200.fasta \
    -task blastn -evalue 0.0001 -outfmt "6 sstart send length"\
    -db {fasta} -dust no -soft_masking false -perc_identity 80 -num_threads 10 \
    -out {outpath}/first_200.m8')
    os.system(f'blastn -query {outpath}/last_100.fasta \
    -task blastn -evalue 0.0001 -outfmt "6 sstart send length"\
    -db {fasta} -dust no -soft_masking false -perc_identity 80 -num_threads 10 \
    -out {outpath}/last_100.m8')
    
    with open(f'{outpath}/first_200.m8', 'r') as first, \
    open(f'{outpath}/last_100.m8', 'r') as last:
        count_first = 0
        count_last = 0
        for line in first:
            line = line.split('\t')
            if int(line[0]) < int(line[1]) and int(line[2]) > 100:
                count_first += 1
        for line in last:
            line = line.split('\t')
            if int(line[0]) < int(line[1]) and int(line[2]) > 50:
                count_last += 1
        if count_first >= 2 and count_last >= 2:
            real = True
    return real

def check_DIRS_ends(genome, tmp, prefix, domain_name):
    
    true_DIRS = []
    
    os.system(f'bedtools getfasta -fi {genome} -bed {tmp}/{prefix}_{domain_name}.element.nested.bed -fo {tmp}/{prefix}_{domain_name}.element.nested.fasta -s -nameOnly')
    
    for seq in SeqIO.parse(f'{tmp}/{prefix}_{domain_name}.element.nested.fasta', 'fasta'):
        with open(f'{tmp}/{prefix}_{domain_name}.element.fasta', 'w') as new:
            new.write(f'>{seq.id}\n{seq.seq}\n')
            
        res = check_if_LTR_DIRS(f'{tmp}/{prefix}_{domain_name}.element.fasta', 
                          tmp)
        te_name = seq.id.split('(')[0]
        
        if res == True:
            true_DIRS.append(te_name)
    
    lines_to_save = []
    
    with open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'r') as old:
        for line in old:
            line1 = line.split('\t')
            if line1[3] in true_DIRS:
                lines_to_save.append(line)
    
    with open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'w') as new:
        for line in lines_to_save:
            new.write(line)

def find_polyA(sequence, min_length=5, max_mismatches=1):
    state = False
    seq = sequence.upper()
    for i in range(len(seq) - min_length + 1):
        window = seq[i:i + min_length]
        mismatches = sum(1 for base in window if base != "A")
        if mismatches <= max_mismatches:
            state = True
    return state

def check_polyA_tail(genome, tmp, prefix, domain_name):
    
    true_LINEs = []
    
    os.system(f'bedtools getfasta -fi {genome} -bed {tmp}/{prefix}_{domain_name}.element.nested.bed -fo {tmp}/{prefix}_{domain_name}.element.nested.fasta -s -nameOnly')
    
    for seq in SeqIO.parse(f'{tmp}/{prefix}_{domain_name}.element.nested.fasta', 'fasta'):
        res = find_polyA(seq.seq[-20:])
        te_name = seq.id.split('(')[0]
        if res == True:
            true_LINEs.append(te_name)
    
    lines_to_save = []
    
    with open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'r') as old:
        for line in old:
            line1 = line.split('\t')
            if line1[3] in true_LINEs:
                lines_to_save.append(line)
    
    with open(f'{tmp}/{prefix}_{domain_name}.element.nested.bed', 'w') as new:
        for line in lines_to_save:
            new.write(line)

def make_empty_annotation(outpath, prefix, name_to_report):
    os.system(f'touch {outpath}/{prefix}.{name_to_report}.gff3')
        