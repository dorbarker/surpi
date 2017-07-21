import re
import csv
from operator import itemgetter
from pathlib import Path
from typing import List

def create_tab_delimited_table(file_type: str, annotated: Path) -> List[List[str]]:

    def extract_pattern(pattern: str, to_search: str) -> str:

        match_slice = len(pattern) + 2
        split_line = to_search.split('\t')
        contains_pattern = [pattern in elem for elem in split_line]

        try:
            match = split_line[contains_pattern.index(True)]
            return match[match_slice:]
        except IndexError:
            return 'no rank'

    table = []

    with annotated.open('r') as data:

        for line in data:

            columns = line.split('\t')  # type: List[str]

            searches = [extract_pattern(tax, line)
                        for tax in ('species', 'genus', 'family')]

            if file_type == 'SNAP':

                table_line = [columns[0], columns[2]] + searches

            else:

                table_line = columns[:2] + searches

            table.append(table_line)

    return table

def extract_barcodes(table: List[List[str]]) -> List[str]:

    tr = str.maketrans('#/', '\t\t')

    def extract_barcode(line: List[str]) -> str:

        header = line[0]

        _, barcode, *_ = header.translate(tr).split()

        return barcode

    preliminary_barcodes = set(extract_barcode(line) for line in table)

    # TODO: Confirm filtering on N from original source is correct
    barcodes = sorted(filter(lambda x: 'N' not in x, preliminary_barcodes))

    return barcodes

def generate_table(mode: str, barcodes: List[str],
                   table: List[List[str]]) -> List[List[str]]:

    def count_bar(to_search: str, table: List[List[str]]) -> str:
        return str(sum(row.count(to_search) for row in table))

    headers = ['Accession', 'Species', 'Genus', 'Family', '@=contigbarcode']
    idx = headers.index(mode)


    output_table = sorted([row[idx:] for row in table], key=itemgetter(0))
    # uniq.column
    uniq_column = set(row[idx] for row in table)

    output_table.insert(0, headers[idx:])

    for barcode in barcodes:

        # bar.$f.$inputfile.tempsorted
        matching_rows = [row for row in table if barcode in row[0]]

        # bar$f.$inputfile.gi.output
        bar = 'bar{}'.format(barcode)
        output = {name: count_bar(name, matching_rows) for name in uniq_column}
        output[headers[idx]] = bar

        for row in output_table:
            key = row[0]
            row.append(output[key])

    return output_table

def write_table(table: List[List[str]], outpath: Path) -> None:

    joined_rows = ('\t'.join(row) for row in table)

    to_write = '\n'.join(joined_rows)

    outpath.write_text(to_write)

def table_generator(file_type: str, mode: str, annotated: Path) -> None:

    outpath = annotated.with_suffix('.{}_counttable'.format(mode))

    table = create_tab_delimited_table(file_type, annotated)

    barcodes = extract_barcodes(table)

    output_table = generate_table(mode, barcodes, table)

    write_table(output_table, outpath)
