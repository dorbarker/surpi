#!/usr/bin/env python3

import argparse

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--input-type',
                        choices=('fasta', 'fastq'),
                        default='fastq',
                        help='Input filetype [fastq]')

    parser.add_argument('--quality-mode',
                        choices=('illumina', 'sanger'),
                        default='sanger',
                        help='Quality score mode [sanger]')

    parser.add_argument('--adapter-set',
                        choices=('truseq', 'nextera',
                                 'nexsolb', 'nexsoltruseq'),
                        default='truseq',
                        help='Adapter set to use [truseq]')

    # TODO: fix formatting
    parser.add_argument('--verify-fastq',
                        choices=(0, 1, 2, 3),
                        default=1,
                        help='''
                        0 = skip fastq validation
                        1 = validate fastqs; quit on fail
                        2 = validate fastqs; check for unique names; quit on fail
                        3 = validate fastqs; check for unique names; do not quit on failure
                        [1]
                        ''')

    parser.add_argument('--run-mode',
                        choices=('comprehensive', 'fast'),
                        default='comprehensive',
                        help='''
                        Comprehensive mode allows SNAP to NT -> denovo contig assembly -> RAPSearch to Viral proteins or NR
                        Fast mode allows SNAP to curated FAST databases
                        [comprehensive]
                        ''')

    parser.add_argument('--skip-preprocess',
                        action='store_true',
                        help='Skip preprocessing step [off]')

    parser.add_argument('--length-cutoff',
                        type=int,
                        default=50,
                        help='Discard any post-trimming sequence less than this length [50]')

    parser.add_argument('--crop-start',
                        type=int,
                        default=10,
                        help='Prior to SNAP alignment, cropping start position [10]')

    parser.add_argument('--crop-length',
                        type=int,
                        default=10,
                        help='Prior to SNP alignment, length to crop after start position [75]')

    parser.add_argument('--quality-cutoff',
                        type=int,
                        help='Quality cutoff passed to cutadapt\' -q parameter [18]')

    parser.add_argument('--edit-distance',
                        type=int,
                        default=12,
                        help='SNAP edit distance for computational subtraction of host genome [12]')


    # snap_nt interator option only has one choice: 'inline' and so is omitted here

    parser.add_argument('--rapsearch-method',
                        choices=('viral', 'nr'),
                        default='viral',
                        help='RAPSearch database method to use [viral]')

    parser.add_argument('--rapsearch-cutoff',
                        default='1',
                        help='E-value cutoff for RAPSearch. Will parse \
                              notation such as "1e+1" [1]')

    parser.add_argument('--fast-rapsearch',
                        action='store_true',
                        help='Yields a 10-30x speed increase at the cost of sensitivity [off]')

    parser.add_argument('--abyss-kmer',
                        type=int,
                        default=34,
                        help='kmer size for ABySS de novo assembly [34]')

    parser.add_argument('--ignore-barcodes',
                        action='store_true',
                        help='Assemble all barcodes together into a single assembly [off]')

    parser.add_argument('--blastn-evalue',
                        default='1e-15',
                        help='E-value cutoff for BLASTn/ Will parse notation \
                             such as "1e-15" [1e-15]')

    # Reference Data


    return parser.parse_args()

def main():

    args = arguments()


if __name__ == '__main__':
    main()
