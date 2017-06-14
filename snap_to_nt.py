#!/usr/bin/env python3

import os
import argparse
import subprocess
from functools import partial
from pathlib import Path
from multiprocessing import cpu_count
from utilities import run_shell, logtime, annotated_to_fastq
from Bio import SeqIO
from taxonomy_lookup import taxonomy_lookup, table_generator

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', type=Path,
                        help='Input fastq file')

    parser.add_argument('--workdir', type=Path,
                        help='Work directory for SURPI')

    parser.add_argument_group('--cores', type=int,
                              default=cpu_count(),
                              help='CPU cores to use [all]')

    return parser.parse_args()

@logtime('Separate matched and unmatched sequences from SAM file')
def separate_sam_lines(sam_file):
    '''Separates SAM-formatted file.

    Lines where the third column is '*' are appended to matched,
    and appended to unmatched otherwise.
    '''

    matched, unmatched = [], []

    with open(sam_file, 'r') as sam:

        for line in sam:
            if not line.startswith('@'):

                split_line = line.split('\t')

                if split_line[2] == '*':
                    matched.append(line)
                else:
                    unmatched.append(line)

    return matched, unmatched

def write_sam(lines, outpath):
    '''Writes SAM lines generated by separate_sam_lines()'''

    with open(outpath, 'w') as f:

        f.write('\n'.join(lines))

@logtime('SNAP unmatched sequences')
def snap_unmatched(basef, workdir, db_dir, unmatched, cache_reset,
                   edit_distance, cores):

    # edit_distance == d_NT_alignment
    snap_sam = basef.with_suffix('.NT.snap.sam')

    if not snap_sam.exists():

        snap_nt_cmd = ('snap_nt.sh', unmatched, db_dir, cores, cache_reset,
                       edit_distance, 'snap-aligner')

        subprocess.check_call(map(str,  snap_nt_cmd))

    os.rename(str(workdir / unmatched), str(workdir / snap_sam))

def separate_ranked_taxonomy(ranks, annotations, result_dir):
    '''Recursively searches through `ranks` and writes taxonomy results
    to files name after the given rank.

    Ranks must be ordered from broad to narrow,
    e.g. ['Eurkaryota', 'Chordata', 'Mammalia', 'Primates']
    '''

    def condition(line, to_find, to_exclude):
        return to_find in line and to_exclude not in line

    def tax_search_recursor(rank_list):

        to_find, to_exclude = rank_list[-2:]

        test = partial(condition, to_find=to_find, to_exclude=to_exclude)

        result_file = (result_dir / to_find).with_suffix('.annotated')

        sub_annot = (line for line in annotations if test(line))

        result_file.write_text('\n'.join(sub_annot))

        if len(rank_list) > 2:

            rank_list = rank_list[:-1]

            tax_search_recursor(rank_list)

        else:

            return result_file

    ranks_ = ranks.append('loremipsumdolorsitamet')  # ensure never found

    result_file = tax_search_recursor(ranks_)

    return result_file

def ribo_snap(name, mode, cores, ribo_dir):

    ribo_snap = ('ribo_snap_bac_euk.sh', name, mode, cores, ribo_dir)

    subprocess.check_call(map(str, ribo_snap))

def lookup_taxonomy(basef, workdir, tax_db_dir, ribo_dir, cores):

    full_length_fastq = basef.with_suffix('.NT.snap.matched.fulllength.fastq')
    full_length_sam = full_length_fastq.with_suffix('.sam')
    annotated = basef.with_suffix('.annotated')

    if not annotated.exists():
        extract_headers = (
            'extractHeaderFromFastq_ncores.sh', cores,
            workdir / basef.with_suffix('.cutadapt.fastq'),
            workdir / basef.with_suffix('.NT.snap.matched.sam'),
            workdir / basef.with_suffix('.NT.snap.matched.fulllength.fastq'),
            workdir / basef.with_suffix('.NT.snap.unmatched.sam'),
            workdir / basef.with_suffix('.NT.snap.unmatched.fulllength.fastq')
            )

        subprocess.check_call(map(str, extract_headers))
        # sort sam

        with open(workdir / basef.with_suffix('.NT.snap.matched.sam')) as sam:

            lines = sam.readlines()
            sorted_lines = sorted(lines, key=lambda x: x[0])

        # may need to double-check sequence sorting here
        with open(workdir / full_length_fastq) as fastq:
            seq_quals = []
            for rec in SeqIO.parse(fastq, 'fastq'):
                fastq_lines = rec.splitlines()
                seq_quals.append([fastq_lines[1], fastq_lines[3]])

        # paste
        with open(full_length_sam, 'w') as out_sam:

            for line, qual  in zip(sorted_lines, seq_quals):

                to_write = line[:9] + qual + line[11:] + ['\n']

                out_sam.write('\t'.join(to_write))

        taxonomy_lookup(infile=full_length_sam,
                        outfile=annotated,
                        db_dir=tax_db_dir,
                        file_type='sam',
                        seq_type='nucl')

        # check final annotation sort -k 13.7n

    viruses = separate_ranked_taxonomy(ranks=('Viruses',),
                                       annotations=annotated,
                                       result_dir=workdir)

    bacteria = separate_ranked_taxonomy(ranks=('Bacteria',),
                                        annotations=annotated,
                                        result_dir=workdir)

    euks = separate_ranked_taxonomy(ranks=('Eukaryota', 'Chordata',
                                           'Mammalia', 'Primates'),
                                    annotations=annotated,
                                    result_dir=workdir)

    ribo_snap(bacteria, 'BAC', cores, ribo_dir)
    ribo_snap(euks, 'EUK', cores, ribo_dir)

    return viruses, bacteria, euks

def extract_to_fast():

    pass

def snap(sample, workdir, snap_db_dir, tax_db_dir, ribo_dir, cores):

    basef = Path(sample.stem)


    snap_unmatched(basef, workdir, snap_db_dir, unmatched, cache_reset,
                   d_NT_alignment, cores)

    snap_sam = workdir / basef.with_suffix('.NT.snap.sam')

    matched, unmatched = separate_sam_lines(snap_sam)

    write_sam(matched, workdir / basef.with_suffix('.NT.snap.matched.sam'))
    write_sam(unmatched, workdir / basef.with_suffix('.NT.snap.unmatched.sam'))

    viruses, bacteria, euks = lookup_taxonomy(basef, workdir, tax_db_dir,
                                              ribo_dir, cores)

    table_generator(viruses, 'SNAP', 'Y', 'Y', 'Y', 'Y')

    # if comprehensive
    viruses.with_suffix('.fastq').write_text(annotated_to_fastq(viruses))
    # TODO: convert unmatched.fullength.fastq to fasta
        # sorted by sequence length
        # Run crop_reads.csh on sorted.fasta <- temporary
        # Run gt sequniq on sorted.cropped.fasta <- temporary
        # extractAlltoFast on sorted.cropped.uniq.fasta

def main():

    args = arguments()

    #snap(args.sample, args.workdir, args.db_dir, args.cores)

if __name__ == '__main__':
    main()
