#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys

from utils import annotate_domain, cluster_domain_nucl, parse_annotation_clustering, adjast_align_clusters, check_DIRS_ends, generate_final_annot_multiple, make_empty_annotation

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
    
    print('INFO: starting DIRS prediction', flush=True)
    print('INFO: annotation of DIRS-1 type RT-RH-MT-YR domains', flush=True)
    
    domains = ['RT', 'RH', 'MT', 'YR']
    
    for domain in domains:
        annotate_domain(f'/hmm/DIRS_{domain}.hmm', domain, args.genome, args.processors, args.tmp, args.prefix)
        try:
            cluster_domain_nucl(domain, args.genome, args.tmp, args.prefix, args.processors, 0.5, 0.9)
            parse_annotation_clustering(domain, args.tmp, args.prefix)
        except:
            print(f'INFO: intact DIRSs found: 0')
            make_empty_annotation(args.outpath, args.prefix, 'dirs')
            sys.exit(1)
    
    print('INFO: clustering and transposon boundary adjastment', flush=True)
    adjast_align_clusters(args.genome, args.tmp, args.prefix, 'YR', 20000, args.processors)
    
    print('INFO: detecting terminal repeats', flush=True)
    check_DIRS_ends(args.genome, args.tmp, args.prefix, 'YR')

    print('INFO: processing output', flush=True)
    element_count = generate_final_annot_multiple(args.outpath, args.tmp, args.prefix, 'YR', ['RT', 'RH', 'MT'], 'DIRS', 'dirs', 'I', args.processors)
    print(f'INFO: intact DIRSs found: {element_count}')

if __name__ == '__main__':
    main()