#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from utils import mask_annotated_LTR_RTE_free_INT, annotate_domain, cluster_domain_nucl, parse_annotation_clustering, adjast_align_clusters, generate_final_annot_dna, make_empty_annotation

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

    print('INFO: starting DNA-TE prediction', flush=True)
    print('INFO: masking LTR-RTEs in reference genome', flush=True)
    path_to_masked = mask_annotated_LTR_RTE_free_INT(args.genome, args.tmp, args.outpath, args.prefix)
    print('INFO: annotation of transposase domains', flush=True)
    annotate_domain('rexdb-pnas', 'TPase', path_to_masked, args.processors, args.tmp, args.prefix)
    try:
        cluster_domain_nucl('TPase', path_to_masked, args.tmp, args.prefix, args.processors, 0.5, 0.9)
        parse_annotation_clustering('TPase', args.tmp, args.prefix)
    except:
        print(f'INFO: intact DNA-TEs found: 0')
        make_empty_annotation(args.outpath, args.prefix, 'dna')
        sys.exit(1)
        
    print('INFO: clustering and transposon boundary adjastment', flush=True)
    adjast_align_clusters(path_to_masked, args.tmp, args.prefix, 'TPase', 10000, args.processors)
    print('INFO: processing output', flush=True)
    generate_final_annot_dna(path_to_masked, args.outpath, args.tmp, args.prefix, 'TPase', 'DNA_TE', 'dna', 'II', args.processors)

if __name__ == '__main__':
    main()