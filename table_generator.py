import re
from pathlib import Path
from typing import List

def create_tab_delimited_table(mode: str, annotated: Path) -> List[List[str]]:

    def extract_pattern(pattern: str, to_search: str) -> str:

        pattern_ = '{}--(.*)\t'.format(pattern)

        search_result = re.search(pattern_, to_search)

        return search_result.groups()[0] if search_result else ''

    table = []

    with annotated.open('r') as data:

        for line in data:

            columns = line.split('\t')

            searches = [extract_pattern(tax, line)
                        for tax in ('species', 'genus', 'family')]

            if mode == 'SNAP':

                table_line = columns[0] + columns[2] + searches

            else:

                table_line = columns[:2] + searches

            table.append(table_line)

    return table

def extract_barcodes(table: List[List[str]]) -> List[str]:

    tr = str.maketrans('#/', '\t\t')

    def extract_barcode(line: List[str]) -> str:

        header = line[0]

        barcode, *_ = header.translate(tr).split()

        return barcode

    preliminary_barcodes = set(extract_barcode(line) for line in table)

    # TODO: Confirm filtering on N from original source is correct
    barcodes = sorted(filter(lambda x: 'N' not in x, preliminary_barcodes))

    return barcodes

def generate_table(mode: str, barcodes: List[str], table: List[List[str]]):

    headers = ['Accession', 'Species', 'Genus', 'Family', '@=contigbarcode']

    idx = headers.index(mode)

    d = {}

    for barcode in barcodes:

        matching_rows = [row for row in table if barcode in row[0]]

        bar = 'bar{}'.format(barcode)
        output = [bar] + [matching_rows.count(row[idx]) for row in table]

        d[barcode] = {'output': output, 'matching': matching_rows}


