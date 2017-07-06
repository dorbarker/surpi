#!/usr/bin/env python3
'''Retrieve taxonomy information for SAM and BLAST files'''

import argparse
import sqlite3 as sql
import csv
import subprocess

from pathlib import Path
from typing import List, Dict

def get_accession_numbers(path: Path, file_type: str) -> List[str]:
    '''Extracts a non-redundant sorted list of accessions from
    the SAM or BLAST file located at path.
    '''

    acc_column = 2 if file_type == 'sam' else 1
    accs = set()

    with path.open('r') as seq_file:
        for line in seq_file:
            try:
                acc = line.split()[acc_column]
                accs.add(acc)
            except IndexError:
                pass

    return sorted(accs)

def query_dbs(db_dir: Path, seq_type: str, accs: List[str]) -> Dict[str, Dict[str, str]]:
    '''Queries NCBI Taxonomy databases by accession number.

    All databases must be co-located in db_dir.
    '''

    acc_db_name = 'acc_taxid_{}.db'.format(seq_type)
    acc_txid_db = sql.connect(str(db_dir / acc_db_name))
    names_db = sql.connect(str(db_dir / 'names_nodes_scientific.db'))

    taxid_query = 'SELECT taxid FROM acc_taxid WHERE acc = ?'
    names_query = 'SELECT name FROM names WHERE taxid = ?'
    nodes_query = 'SELECT * FROM nodes WHERE taxid = ?'

    taxonomy = {}  # type: Dict[str, Dict]

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

def result_append(in_path: Path, out_path: Path,
                  taxonomy: Dict[str, Dict[str, str]], file_type: str) -> None:
    '''Appends taxonomy info line-by-line to the SAM file and writes the
    result to `out_sam_path`.
    '''

    acc_column = 2 if file_type == 'sam' else 1
    ranks = ('family', 'genus', 'species', 'lineage')

    with in_path.open('r') as infile, out_path.open('w') as outfile:

        out = csv.writer(outfile, delimiter='\t', quoting=csv.QUOTE_NONE)

        for line in infile:

            sam_line = line.strip().split()
            acc = line[acc_column]

            ranks = ['{}--{}'.format(rnk, taxonomy[acc][rnk]) for rnk in ranks]

            out.writerow(sam_line.extend(ranks))

def table_generator(annotated: Path, snap_rap: str, acc: str, species: str,
                    genus: str, family: str) -> None:
    '''Runs external script `table_generator.sh`'''

    assert snap_rap in ('SNAP', 'RAP')
    assert acc in ('Y', 'N')
    assert species in ('Y', 'N')
    assert genus in ('Y', 'N')
    assert family in ('Y', 'N')

    # TODO: MUST remove this diagnostic code
    cmd = ('bash', '/home/dbarker/Projects/surpi/table_generator.sh', str(annotated),
           snap_rap, acc, species, genus, family)

    subprocess.check_call(cmd)

def taxonomy_lookup(infile: Path, outfile: Path, db_dir: Path,
                    seq_type: str, file_type: str) -> None:
    '''Given a SAM file, look up taxonomy information for each sequence,
    append the taxonomy information to the end of each SAM line, and write this
    to a new file.
    '''

    accs = get_accession_numbers(infile, file_type)

    taxonomy = query_dbs(db_dir, seq_type, accs)

    result_append(infile, outfile, taxonomy, file_type)
