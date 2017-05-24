#!/usr/bin/env python3

import argparse
import sqlite3 as sql

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help='blast_file/sam_file')

    parser.add_argument('file_type', choices = ('blast', 'sam'),
                        help='File format')

    parser.add_argument('seq_type', choices = ('nucl', 'prot'),
                        help='Sequence type')

    parser.add_argument('cores', type=int, help='CPU cores to use')

    parser.add_argument('tax_dir', help='Taxonomy reference directory')

    return parser.parse_args()

def get_accession_numbers(path, file_type):
    '''Extracts a non-redundant sorted list of accessions from
    the SAM or BLAST file located at path.
    '''

    acc_column = 2 if file_type == 'sam' else 1
    accs = set()

    with open(path, 'r') as f:
        for line in f:
            try:
                accs.add( line.split()[acc_column] )
            except IndexError:
                pass

    return sorted(accs)

def query_dbs(db_dir, seq_type, accs):

    acc_db_name  = 'acc_taxid_{}.db'.format(seq_type)
    acc_txid_db = sql.connect(str(db_dir / acc_db_name))
    names_db = sql.connect(str(db_dir / 'names_nodes_scientific.db'))

    taxid_query = 'SELECT taxid FROM acc_taxid WHERE acc = ?'
    names_query = 'SELECT name FROM names WHERE taxid = ?'
    nodes_query = 'SELECT * FROM nodes WHERE taxid = ?'

    with acc_txid_db, names_db:

        for acc in accs:

            acc_txid_cur = acc_txid_db.execute(taxid_query, (acc, ))

            taxid, *_ = acc_txid_cur.fetchone()  # one taxid per accession

            print(taxid)

            names_cur = names_db.execute(names_query, (taxid, ))

            print(names_cur.fetchall())




def main():

    args = arguments()

if __name__ == '__main__':
    main()
