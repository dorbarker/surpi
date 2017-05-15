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

    tar = ('tar', 'xzf', str(db_dir / 'taxdump.tar.gz'))
    subprocess.call(tar)

    for gz in db_dir.glob('*.gz'):
        cmd = ('pigz -dck {} > {}'.format(gz, gz.with_suffix('').name))
        subprocess.call(cmd, shell=True)

def verify_files(db_dir):

    # TODO add file names
    files = ()

    if not all((db_dir / f).exists() for f in files):
        user_msg('Did not find all files. Quitting...', file=sys.stderr)
        sys.exit(1)

def trim_names():

    subprocess.call('grep "scientific name" names.dmp > names_scientificname.dmp',
                    shell=True)

def tidy():

    for i in Path('.').glob('*.dmp'):
        os.remove(i)

    os.remove('gc.prt')
    os.remove('readme.txt')

def create_dbs():

    # Create names_nodes_scientific.db
    user_msg("Creating names_nodes_scientific.db...")
    conn = sqlite3.connect('names_nodes_scientific.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE names (
                taxid INTEGER PRIMARY KEY,
                name TEXT)''')

    with open('names_scientificname.dmp', 'r') as map_file:
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

    with open('nodes.dmp', 'r') as map_file:
        for line in map_file:
            line = line.split("|")
            taxid = line[0].strip()
            parent_taxid = line[1].strip()
            rank = line[2].strip()

            d.execute ("INSERT INTO nodes VALUES (?,?,?)", (taxid, parent_taxid, rank))
    conn.commit()
    conn.close()

    # Create gi_taxid_nucl.db
    user_msg("Creating gi_taxid_nucl.db...")
    conn = sqlite3.connect('gi_taxid_nucl.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE gi_taxid (
                gi INTEGER PRIMARY KEY,
                taxid INTEGER)''')

    with open('gi_taxid_nucl.dmp', 'r') as map_file:
        for line in map_file:
            fst, snd, *rest = line.split()
            c.execute('INSERT INTO gi_taxid VALUES ({},{})'.format(fst, snd))
    conn.commit()
    conn.close()

    # Create gi_taxid_prot.db
    user_msg("Creating gi_taxid_prot.db...")
    conn = sqlite3.connect('gi_taxid_prot.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE gi_taxid (
                gi INTEGER PRIMARY KEY,
                taxid INTEGER)''')

    with open('gi_taxid_prot.dmp', 'r') as map_file:
        for line in map_file:
            fst, snd, *rest  = line.split()
            c.execute('INSERT INTO gi_taxid VALUES ({},{})'.format(fst, snd))

    conn.commit()
    conn.close()

def main():

    args = arguments()

    verify_files(args.db_directory)

    unzip_downloads(args.db_directory)

    trim_names()

    create_dbs()

    tidy()

if __name__ == '__main__':
    main()
