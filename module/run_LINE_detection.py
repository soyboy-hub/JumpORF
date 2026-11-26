#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from utils import mask_annotated_LTR_RTE_free_INT, mask_annotated_DNA_TE, mask_annotated_Helitron_TE, annotate_domain, cluster_domain_nucl, parse_annotation_clustering, adjast_align_clusters, check_polyA_tail, generate_final_annot_dual, make_empty_annotation

def get_args():
    desc = (
        """"""
    )
    epi = """"""
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("genome", action="store", help='path to genome .fasta')
    parser.add_argument("processors", action="store", help='processors')
    parser.add_argument("outpath", action="store", help='outpath')
    parser.add_argument("prefix", action="store", help='prefix')
    parser.add_argument("tmp", action="store", help='path to tmp dir')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def main():
    args = get_args()
    
    print('INFO: starting LINE prediction', flush=True)
    print('INFO: masking LTR-RTEs in reference genome', flush=True)
    path_to_masked_ltr = mask_annotated_LTR_RTE_free_INT(args.genome, args.tmp, args.outpath, args.prefix)
    print('INFO: masking DNA-TEs in reference genome', flush=True)
    path_to_masked_ltr_dna = mask_annotated_DNA_TE(path_to_masked_ltr, args.tmp, args.outpath, args.prefix)
    print('INFO: masking DNA-TEs in reference genome', flush=True)
    path_to_masked_ltr_dna_hel = mask_annotated_Helitron_TE(path_to_masked_ltr_dna, args.tmp, args.outpath, args.prefix)
    print('INFO: annotation of LINE RT-EN domains', flush=True)
    
    domains = ['RT', 'EN']
    
    for domain in domains:
        annotate_domain(f'/hmm/LINE_{domain}.hmm', domain, path_to_masked_ltr_dna_hel, args.processors, args.tmp, args.prefix)
        try:
            cluster_domain_nucl(domain, path_to_masked_ltr_dna_hel, args.tmp, args.prefix, args.processors, 0.5, 0.9)
            parse_annotation_clustering(domain, args.tmp, args.prefix)
        except:
            print(f'INFO: intact LINEs found: 0')
            make_empty_annotation(args.outpath, args.prefix, 'line')
            sys.exit(1)
    
    print('INFO: clustering and transposon boundary adjastment', flush=True)
    adjast_align_clusters(path_to_masked_ltr_dna_hel, args.tmp, args.prefix, 'RT', 10000, args.processors)
    print('INFO: checking polyA tail')
    check_polyA_tail(path_to_masked_ltr_dna_hel, args.tmp, args.prefix, 'RT')
    print('INFO: processing output', flush=True)
    element_count = generate_final_annot_dual(args.outpath, args.tmp, args.prefix, 'RT', 'EN', 'LINE', 'line', 'I', args.processors)
    print(f'INFO: intact LINEs found: {element_count}')

if __name__ == '__main__':
    main()