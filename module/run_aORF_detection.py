#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from utils import check_if_annotation, extract_TEs_and_basal_proteins, run_augustus_on_batches, process_augustus_out, remove_basal_protein_from_pred, extract_transcripts_from_gff3, run_trans_decoder, generate_final_annotation

def get_args():
    desc = (
        """"""
    )
    epi = """"""
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("genome", action="store", help='path to genome .fasta')
    parser.add_argument("specie", action="store", help='specie from augustus list')
    parser.add_argument("min_identity_to_cluster", action="store", help='min_identity_to_cluster')
    parser.add_argument("min_size_cluster", action="store", help='min_size_cluster')
    parser.add_argument("min_coverage", action="store", help='min_coverage')
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
    
    list_have_annotation = check_if_annotation(['ltr', 'dna', 'helitron', 'line', 'dirs', 'ple'], args.outpath, args.prefix)
    extract_TEs_and_basal_proteins(list_have_annotation, args.genome, args.outpath, args.tmp, args.prefix)
    run_augustus_on_batches(args.tmp, args.prefix, 20, args.specie, args.processors)
    process_augustus_out(args.tmp, args.prefix)
    remove_basal_protein_from_pred(args.outpath, args.tmp, args.prefix)
    extract_transcripts_from_gff3(args.outpath, args.tmp, args.prefix, args.genome)
    run_trans_decoder(args.tmp, args.prefix)
    generate_final_annotation(args.outpath, args.tmp, args.prefix, args.processors, args.min_identity_to_cluster, args.min_size_cluster, args.min_coverage)

if __name__ == '__main__':
    main()