#!/usr/bin/env python3

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
from pathlib import Path
import shutil
import time

def get_args():
    desc = (
        """"""
    )
    epi = """"""
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("genome", action="store", help='path to the genome (.fasta / .fna / .fa)')
    parser.add_argument("outpath", action="store", help='path to output directory')
    parser.add_argument("prefix", action="store", help='prefix, e.g. ArabThal')
    parser.add_argument("processors", action="store", help='processors number')
    parser.add_argument("specie", action="store", help='specie from augustus list')
    parser.add_argument("min_identity_to_cluster", action="store", help='min_identity_to_cluster')
    parser.add_argument("min_size_cluster", action="store", help='min_size_cluster')
    parser.add_argument("min_coverage", action="store", help='min_coverage')
    parser.add_argument("fallback_mode", action="store", help='add dante_ltr fallback mode')
    parser.add_argument("--overwrite", action="store_true", help='overwright existing results')
    parser.add_argument("--debug", action="store_true", help='do not clean up temporary directories for debugging purposes')
    
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def main():
    args = get_args()

    outdir = Path(args.outpath)
    outdir_tmp = Path(f'{args.outpath}/tmp')
    if outdir.exists():
        if not args.overwrite:
            print(f'WARNING: output directory already exists: {outdir}\n'
              f'Use --overwrite or choose a new path.')
            sys.exit(1)
        if args.overwrite:
            print(f'INFO: Overwriting existing directory: {outdir}')
            shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outdir_tmp.mkdir(parents=True, exist_ok=True)

    print(r'      _                        ___   _____  _____ ', flush=True)
    print(r'     | |                      / _ \ |  __ ||  ___)', flush=True)
    print(r'     | |_   _ ________   ___ | | | || |__) | |__                       _________', flush=True)
    print(r' _   | | | | |  _   _ \ / _ \| | | ||  _  /|  __)          _______    /     ____\____', flush=True)
    print(r"| |__| | | | | | | | | | (_) | |_| || | \ \| |     _____  /       \  /     |_________|->  ", flush=True)
    print(r" \____/ \___/|_| |_| |_|  __/ \___/ |_|  \_\_|    /     \/         \/", flush=True)
    print(r'                       | |', flush=True)
    print(r'                       |_|', flush=True)

    #run LTR-RTE detection module
    start = time.time()
    if args.fallback_mode == 'yes':
        os.system(f'/module/run_LTR_RTE_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp} yes')
    if args.fallback_mode == 'no':
        os.system(f'/module/run_LTR_RTE_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp} no')
    #run DNA-TE detection module
    os.system(f'/module/run_DNA_TE_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    #run Helitron detection module
    os.system(f'/module/run_Helitron_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    #run LINE detection module
    os.system(f'/module/run_LINE_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    #run DIRS detection module
    os.system(f'/module/run_DIRS_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    #run Penelope-like detection module
    os.system(f'/module/run_PLE_detection.py {args.genome} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    #run aORF detection module
    os.system(f'/module/run_aORF_detection.py {args.genome} {args.specie} {args.min_identity_to_cluster} {args.min_size_cluster} {args.min_coverage} {args.processors} {args.outpath} {args.prefix} {outdir_tmp}')
    end = time.time()
    elapsed = end - start
    print(f'INFO: Pipeline finished successfully. Total execution time: {elapsed:.2f} seconds"')
    
    if not args.debug:
        shutil.rmtree(outdir_tmp)

if __name__ == '__main__':
    main()
    