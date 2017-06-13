#!/usr/bin/env python3
'''Retrieve taxonomy information for SAM and BLAST files'''

import argparse
import sqlite3 as sql
import csv

def arguments():
    '''Command line argument parser for gather arguments when called
    externally.
    '''

    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help='blast_file/sam_file')

    parser.add_argument('output_file', help='blast_file/sam_file')

    parser.add_argument('file_type', choices=('blast', 'sam'),
                        help='File format')

    parser.add_argument('seq_type', choices=('nucl', 'prot'),
                        help='Sequence type')

    parser.add_argument('tax_dir', help='Taxonomy reference directory')

    return parser.parse_args()

def get_accession_numbers(path, file_type):
    '''Extracts a non-redundant sorted list of accessions from
    the SAM or BLAST file located at path.
    '''

    acc_column = 2 if file_type == 'sam' else 1
    accs = set()

    with open(path, 'r') as seq_file:
        for line in seq_file:
            try:
                acc = line.split()[acc_column]
                accs.add(acc)
            except IndexError:
                pass

    return sorted(accs)

def query_dbs(db_dir, seq_type, accs):
    '''Queries NCBI Taxonomy databases by accession number.

    All databases must be co-located in db_dir.
    '''

    acc_db_name = 'acc_taxid_{}.db'.format(seq_type)
    acc_txid_db = sql.connect(str(db_dir / acc_db_name))
    names_db = sql.connect(str(db_dir / 'names_nodes_scientific.db'))

    taxid_query = 'SELECT taxid FROM acc_taxid WHERE acc = ?'
    names_query = 'SELECT name FROM names WHERE taxid = ?'
    nodes_query = 'SELECT * FROM nodes WHERE taxid = ?'

    taxonomy = {}

    with acc_txid_db, names_db:

        for acc in accs:

            taxonomy[acc] = {}

            acc_txid_cur = acc_txid_db.execute(taxid_query, (acc, ))

            taxid, *_ = acc_txid_cur.fetchone()  # one taxid per accession

            while taxid is not 1:

                names_cur = names_db.execute(names_query, (taxid, ))

                name, *_ = names_cur.fetchone()

                nodes_cur = names_db.execute(nodes_query, (taxid, ))

                _, taxid, rank = nodes_cur.fetchone()

                taxonomy[acc][rank] = name

    return taxonomy

def sam_append(in_sam_path, out_sam_path, taxonomy):
    '''Appends taxonomy info line-by-line to the SAM file and writes the
    result to `out_sam_path`.
    '''

    ranks = ('family', 'genus', 'species', 'lineage')

    with open(in_sam_path, 'r') as in_sam, open(out_sam_path, 'w') as out_sam:

        out = csv.writer(out_sam, delimiter='\t', quoting=csv.QUOTE_NONE)

        for line in in_sam:

            sam_line = line.strip().split()
            acc = line[2]

            ranks = ('{}--{}'.format(rnk, taxonomy[acc][rnk]) for rnk in ranks)

            out.writerow(sam_line.extend(ranks))


def taxonomy_lookup(in_sam, out_sam, db_dir, seq_type):
    '''Given a SAM file, look up taxonomy information for each sequence,
    append the taxonomy information to the end of each SAM line, and write this
    to a new file.
    '''

    # TODO: add support for BLAST files
    accs = get_accession_numbers(in_sam, 'sam')

    taxonomy = query_dbs(db_dir, seq_type, accs)

    sam_append(in_sam, out_sam, taxonomy)

def main():
    '''Main function for calling taxonomy_lookup() directly from the command
    line.

    However, this module was designed to be imported into surpi.py
    '''

    args = arguments()

    taxonomy_lookup(args.input_file, args.output_file,
                    args.tax_dir, args.seq_type)


if __name__ == '__main__':
    main()
