#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from utils import annotate_domain, cluster_domain_nucl, parse_annotation_clustering, adjast_align_clusters, generate_final_annot_dual, make_empty_annotation

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
    
    print('INFO: starting Penelope-like prediction', flush=True)
    print('INFO: annotation of RT-EN domains', flush=True)
    
    domains = ['RT', 'EN']
    
    for domain in domains:
        annotate_domain(f'/hmm/PLE_{domain}.hmm', domain, args.genome, args.processors, args.tmp, args.prefix)
        try:
            cluster_domain_nucl(domain, args.genome, args.tmp, args.prefix, args.processors, 0.5, 0.9)
            parse_annotation_clustering(domain, args.tmp, args.prefix)
        except:
            print(f'INFO: intact Penelope-like elements found: 0')
            make_empty_annotation(args.outpath, args.prefix, 'ple')
            sys.exit(1)
        
    print('INFO: clustering and transposon boundary adjastment', flush=True)
    adjast_align_clusters(args.genome, args.tmp, args.prefix, 'RT', 10000, args.processors)
    print('INFO: processing output', flush=True)
    element_count = generate_final_annot_dual(args.outpath, args.tmp, args.prefix, 'RT', 'EN', 'Penelope', 'ple', 'I', args.processors)
    print(f'INFO: intact Penelope-like elements found: {element_count}')

if __name__ == '__main__':
    main()