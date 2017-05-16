#!/usr/bin/env python3
#
#	create_taxonomy_db.py
#
#	This program creates the SQLite taxonomy database used by SURPI
#	Chiu Laboratory
#	University of California, San Francisco
#	January, 2014
#
# Copyright (C) 2014 Scot Federman - All Rights Reserved
# SURPI has been released under a modified BSD license.
# Please see license file for details.
# Last revised 7/2/2014
#
# Modified by Dillon Barker 2017-05-15

import argparse
import sqlite3
import sys
from pathlib import Path
import subprocess
import os

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('db_directory',
                        type=Path,
                        help='Directory containing NCBI data')

    return parser.parse_args()

def user_msg(*args):

    print(*args, file=sys.stderr)

def unzip_downloads(db_dir):

    tar = ('tar', 'xzf', str(db_dir / 'taxdump.tar.gz'), '-C', str(db_dir))
    subprocess.call(tar)

    for gz in db_dir.glob('*.gz'):
        cmd = ('pigz -dck {} > {}'.format(gz, db_dir / gz.with_suffix('').name))
        subprocess.call(cmd, shell=True)

def verify_files(db_dir):

    files = ('nucl_est.accession2taxid.gz', 'nucl_wgs.accession2taxid.gz',
             'nucl_gb.accession2taxid.gz', 'nucl_gss.accession2taxid.gz',
             'prot.accession2taxid', 'nt.gz', 'nr.gz', 'taxdump.tar.gz')

    if not all((db_dir / f).exists() for f in files):
        user_msg('Did not find all files. Quitting...', file=sys.stderr)
        sys.exit(1)

def trim_names(db_dir):

    names = db_dir / 'names.dmp'
    scinames = db_dir / 'names_scientificname.dmp'

    grep_scinames = ('grep', 'scientific name', str(names))

    scinames_lines = subprocess.check_output(grep_scinames,
                                             universal_newlines=True)

    with open(str(scinames), 'w') as f:
        f.write(scinames_lines)

def tidy(db_dir):

    for i in Path(db_dir).glob('*.dmp'):
        os.remove(i)

    os.remove('gc.prt')
    os.remove('readme.txt')

def create_dbs(db_dir):

    def acc_taxid(path, db_path):

        c = sqlite3.connect(str(db_path))

        with open(str(path), 'r') as f:
            for line in f:
                if line.startswith('accession'):
                    continue

                acc, acc_ver, taxid, gi = line.strip().split()

                cmd = 'INSERT INTO acc_taxid VALUES (?,?)'
                c.execute(cmd, (acc, taxid))
        c.commit()
        c.close()

    # Create names_nodes_scientific.db
    user_msg("Creating names_nodes_scientific.db...")
    conn = sqlite3.connect(str(db_dir / 'names_nodes_scientific.db'))
    c = conn.cursor()
    c.execute('''CREATE TABLE names (
                taxid INTEGER PRIMARY KEY,
                name TEXT)''')

    with open(str(db_dir / 'names_scientificname.dmp'), 'r') as map_file:
        for line in map_file:
            line = line.split("|")
            taxid = line[0].strip()
            name = line[1].strip()

            c.execute ("INSERT INTO names VALUES (?,?)", (taxid, name))

    d = conn.cursor()
    d.execute('''CREATE TABLE nodes (
                taxid INTEGER PRIMARY KEY,
                parent_taxid INTEGER,
                rank TEXT)''')

    with open(str(db_dir / 'nodes.dmp'), 'r') as map_file:
        for line in map_file:
            line = line.split("|")
            taxid = line[0].strip()
            parent_taxid = line[1].strip()
            rank = line[2].strip()

            d.execute ("INSERT INTO nodes VALUES (?,?,?)", (taxid, parent_taxid, rank))
    conn.commit()
    conn.close()

    # Create acc_taxid_nucl.db
    user_msg('Creating acc_taxid_nucl.db...')
    conn = sqlite3.connect(str(db_dir / 'acc_taxid_nucl.db'))
    c = conn.cursor()
    c.execute('''CREATE TABLE acc_taxid (
                 acc TEXT PRIMARY KEY,
                 taxid integer)''')

    # insert values to acc_taxid_nucl.db
    for i in ('est', 'wgs', 'gb', 'gss'):

        user_msg('Adding {} to acc_taxid_nucl.db'.format(i))

        path = db_dir / 'nucl_{}.accession2taxid'.format(i)

        acc_taxid(path, db_dir / 'acc_taxid_nucl.db')

    # Create gi_taxid_prot.db
    user_msg("Creating acc_taxid_prot.db...")
    conn = sqlite3.connect(str(db_dir / 'acc_taxid_prot.db'))
    c = conn.cursor()
    c.execute('''CREATE TABLE acc_taxid (
                acc TEXT PRIMARY KEY,
                taxid INTEGER)''')

    acc_taxid(db_dir / 'prot.accession2taxid', db_dir / 'acc_taxid_prot.db')

def main():

    args = arguments()

    verify_files(args.db_directory)

    #unzip_downloads(args.db_directory)

    trim_names(args.db_directory)

    create_dbs(args.db_directory)

    tidy(args.db_directory)

if __name__ == '__main__':
    main()
