#!/usr/bin/env python3

import sys
import argparse
from datetime import datetime

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help='blast_file/sam_file')

    parser.add_argument('file_type', choices = ('blast', 'sam'),
                        help='File format')

    parser.add_argument('seq_type', choices = ('nucl', 'prot'),
                        help='Sequences type')

    parser.add_argument('cores', type=int, help='CPU cores to use')

    parser.add_argument('tax_dir', help='Taxonomy reference directory')

    return parser.parse_args()

def main():

    args = arguments()

    begintime = datetime.now()

if __name__ == '__main__':
    main()
