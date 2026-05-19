#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from utils import parse_dante_ltr

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
    parser.add_argument("fallback_mode", action="store", help='fallback_mode')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def run_dante(genome, tmp, prefix, processors):
    os.system(f'dante -q {genome} -D Viridiplantae_v4.0 -o {tmp}/{prefix}.domain_dante.gff3 -M BL80 -c {processors} > /dev/null 2>&1')
    
def run_dante_ltr(genome, tmp, prefix, processors):
    os.system(f'JumpORF_container/dante_ltr/dante_ltr -g {tmp}/{prefix}.domain_dante.gff3 -s {genome} -o {tmp}/{prefix}.dante_ltr -c {processors} -M 1 > /dev/null 2>&1')

def run_dante_ltr_fallback_mode(genome, tmp, prefix, processors):
    os.system(f'JumpORF_container/dante_ltr/dante_ltr -g {tmp}/{prefix}.domain_dante.gff3 -s {genome} -o {tmp}/{prefix}.dante_ltr -c {processors} -M 1 --fallback_mode coarse2 > /dev/null 2>&1')

def main():
    args = get_args()
    
    print('INFO: starting LTR-RTE prediction', flush=True)
    
    run_dante(args.genome, args.tmp, args.prefix, args.processors)
    if args.fallback_mode == 'no':
        run_dante_ltr(args.genome, args.tmp, args.prefix, args.processors)
    if args.fallback_mode == 'yes':
        run_dante_ltr_fallback_mode(args.genome, args.tmp, args.prefix, args.processors)
        
    intact_elements = parse_dante_ltr(args.tmp, args.outpath, args.prefix)

    print(f'INFO: intact LTR-RTEs found: {intact_elements}', flush=True)
    
if __name__ == '__main__':
    main()