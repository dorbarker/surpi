#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from multiprocessing import cpu_count
from typing import List, Tuple, Union
from utilities import logtime, annotated_to_fastq, fastq_to_fasta, sam_to_fasta, user_msg
from Bio import SeqIO
from taxonomy_lookup import taxonomy_lookup
from table_generator import table_generator
from preprocessing import crop_reads
import update_sam

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
def separate_sam_lines(sam_file: Path) -> Tuple[List[str], List[str]]:
    '''Separates SAM-formatted file.

    Lines where the third column is '*' are appended to matched,
    and appended to unmatched otherwise.
    '''

    matched, unmatched = [], []

    with sam_file.open('r') as sam:
        for line in sam:
            if not line.startswith('@') and not line.isspace():

                split_line = line.split('\t')

                if split_line[2] != '*':
                    matched.append(line)
                else:
                    unmatched.append(line)

    return matched, unmatched

def write_sam(lines: List[str], outpath: Path) -> None:
    '''Writes SAM lines generated by separate_sam_lines()'''

    with outpath.open('w') as out:

        out.write('\n'.join(lines))

@logtime('SNAP unmatched sequences')
def snap_unmatched(infile: Path, workdir: Path, tempdir: Path,
                   snap_db_dir: Path, edit_distance: int, cores: int) -> Path:

    first_pass_done = False
    snap_sam = infile.with_suffix('.NT.snap.sam')

    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:

        throw_away = Path(ephemeral)

        tmp_fastq = (throw_away / infile.stem).with_suffix('.tmp.fastq')
        tmp_sam = tmp_fastq.with_suffix('.sam')
        prev_sam = tmp_sam.with_suffix('.prev')

        for snap_db in snap_db_dir.glob('*'):

            if not first_pass_done:

                snap_cmd = ('snap-aligner', 'single', snap_db, infile,
                            '-o', tmp_sam, '-t', cores, '-x', '-f', '-h', 250,
                            '-d', edit_distance, '-n', 25)

            else:

                snap_cmd = ('snap-aligner', 'single', snap_db, tmp_fastq,
                            '-o', tmp_sam, '-t', cores, '-x', '-f', '-h', 250,
                            '-d', edit_distance, '-n', 25)

            subprocess.check_call([str(arg) for arg in snap_cmd if arg])

            update_sam.compare_sam(tmp_sam, prev_sam)
            tmp_fastq.write_text(annotated_to_fastq(prev_sam))

            first_pass_done = True

        # *.NT.sam doesn't appear anywhere else in the original SURPI - WTF?
            # Using *.NT.snap.sam instead


        update_sam.update_sam(prev_sam, snap_sam)

    return snap_sam

def separate_ranked_taxonomy(ranks: List[str], annotations: List[str],
                             result_dir: Path) -> Path:
    '''Recursively searches through `ranks` and writes taxonomy results
    to files name after the given rank.

    Ranks must be ordered from broad to narrow,
    e.g. ['Eurkaryota', 'Chordata', 'Mammalia', 'Primates']
    '''

    def condition(line: str, to_find: str, to_exclude: str) -> bool:
        '''Returns True if the `to_find` pattern is in `line` and
        `to_exclude` is not in `line`
        '''
        return to_find in line and to_exclude not in line

    def tax_search_recursor(rank_list) -> Path:

        to_find, to_exclude = rank_list[-2:]

        test = partial(condition, to_find=to_find, to_exclude=to_exclude)

        result_file = (result_dir / to_find).with_suffix('.annotated')

        sub_annot = (line for line in annotations if test(line))

        result_file.write_text('\n'.join(sub_annot))

        if len(rank_list) > 2:

            rank_list = rank_list[:-1]

            result = tax_search_recursor(rank_list)

        else:

            result = result_file

        return result

    ranks_ = ranks + ['loremipsumdolorsitamet']  # ensure never found

    result_file = tax_search_recursor(ranks_)

    return result_file

def ribo_snap(inputfile: Path, mode: str, cores: int, ribo_dir: Path, tempdir: Path):
    '''Runs SNAP to subtract `inputfile` reads from databases
    of bacterial and ribosomal reads
    '''

    def extract_sam_from_sam(infile, outfile, riboremoved):
        '''Writes to `outfile` the reads in `infile` whose headers
        are found in `riboremoved`
        '''

        with riboremoved.open('r') as noribo:
            headers = set(line.split()[0] for line in noribo)

        with infile.open('r') as inf:
            lines = [line.strip() for line in inf if line.split()[0] in headers]

        outfile.write_text('\n'.join(lines))

    noribo = inputfile.with_suffix('.noribo.annotated')

    if mode == 'BAC':

        snap_large = ribo_dir / '23s.fa.snap_index'
        snap_small = ribo_dir / 'rdp_typed_iso_goodq_9210seqs.fa.snap_index'

    else:

        snap_large = ribo_dir / \
                '28s_rRNA_gene_NOT_partial_18s_spacer_5.8s.fa.snap_index'
        snap_small = ribo_dir / '18s_rRNA_gene_not_partial.fa.snap_index'

    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:

        throw_away = Path(ephemeral)

        fasta = (throw_away / inputfile.stem).with_suffix('.fasta')
        fastq = fasta.with_suffix('.fastq')
        crop = fastq.with_suffix('.crop.fastq')
        large_out_sam = crop.with_suffix('.noLargeS.unmatched.sam')
        large_out_fq = large_out_sam.with_suffix('.fastq')
        small_out = crop.with_suffix('.noSmallS_LargeS.sam')

        sam_to_fasta(inputfile, fasta)

        fastq.write_text(subprocess.check_output(('fasta_to_fastq', str(fasta)),
                                                 universal_newlines=True))

        crop_reads(fastq, crop, 10, 75)

        snap1 = ('snap-aligner', 'single', snap_large, crop, '-o', large_out_sam,
                 '-t', cores, '-f', '-h', 250, '-d', 18, '-n', 200, '-F', 'u')

        subprocess.check_call([str(arg) for arg in snap1])

        # conversion of sam -> fastq
        large_out_fq.write_text(annotated_to_fastq(large_out_sam))

        snap2 = ('snap-aligner', 'single', snap_small, large_out_fq,
                 '-o', small_out, '-t', cores, '-h', 250, '-d', 18, '-n', 200,
                 '-F', 'u')

        subprocess.check_call([str(arg) for arg in snap2])

        extract_sam_from_sam(inputfile, noribo, small_out)


    table_generator('SNAP', 'Genus', noribo)

def extract_headers(parentfile: Path, queryfile: Path, output: Path):

    with queryfile.open('r') as query:
        headers = set()
        for line in query:
            try:
                headers.add(line.split()[0])
            except IndexError:
                pass  # Skip empty lines

    with parentfile.open('r') as parent, output.open('w') as out:
        for rec in SeqIO.parse(parent, 'fastq'):
            if rec.id in headers:
                SeqIO.write(rec, out, 'fastq')

def lookup_taxonomy(matched_fastq: Path, matched_sam: Path,  annotated: Path,
                    workdir: Path, tax_db_dir: Path) -> Tuple[Path, Path, Path]:
    '''Looks up the taxonomy of annotated reads'''

    def insert_seq_quals_to_line(line: List[str], seq_qual: List[str]) -> str:

        line_tuple = (*line[:9], *seq_qual, *line[11:])
        return '\t'.join(line_tuple)

    full_length_sam = matched_fastq.with_suffix('.sam')

    # sort sam
    with matched_sam.open('r') as sam:
        lines = []
        for line in sam:

            lines.append(line.strip().split())

        sorted_lines = sorted(filter(None, lines), key=lambda x: x[0])
        #lines = sam.readlines()
        #sorted_lines = sorted(lines, key=lambda x: x.split()[0])

    # may need to double-check sequence sorting here
    with matched_fastq.open('r') as fastq:

        fastq_ids, seq_quals = [], []
        for rec in SeqIO.parse(fastq, 'fastq'):

            fastq_ids.append(rec.id)

            fastq_lines = rec.format('fastq').splitlines()
            seq_quals.append((fastq_lines[1], fastq_lines[3]))

    # sort seq_quals by their fastq sequence name
    _,  seq_quals = zip(*sorted(zip(fastq_ids, seq_quals)))

    # paste
    lines_to_write = (insert_seq_quals_to_line(line, seq_qual)
                      for line, seq_qual in zip(sorted_lines, seq_quals))

    full_length_sam.write_text('\n'.join(lines_to_write))

    taxonomy_lookup(infile=full_length_sam,
                    outfile=annotated,
                    db_dir=tax_db_dir,
                    file_type='sam',
                    seq_type='nucl')

    # TODO: check final annotation sort -k 13.7n

    annotations = annotated.read_text().splitlines()

    viruses = separate_ranked_taxonomy(ranks=['Viruses'],
                                       annotations=annotations,
                                       result_dir=workdir)

    bacteria = separate_ranked_taxonomy(ranks=['Bacteria'],
                                        annotations=annotations,
                                        result_dir=workdir)

    euks = separate_ranked_taxonomy(ranks=['Eukaryota', 'Chordata',
                                           'Mammalia', 'Primates'],
                                    annotations=annotations,
                                    result_dir=workdir)

    return viruses, bacteria, euks

def extract_to_fast(fastq: Path, fasta: Path, output: Path, tempdir: Path):

    def subseq(infile: Path, parent: Path, outfile: Path):

        header_file = infile.with_suffix('.headers')
        headers = [line[1:] for line in infile.read_text().splitlines()
                   if line.startswith('>')]

        header_file.write_text('\n'.join(headers))

        cmd = ('seqtk', 'subseq', str(parent), str(header_file))

        output = subprocess.check_output(cmd, universal_newlines=True)

        outfile.write_text(output)

        header_file.unlink()

    def sequniq(infile: Path, outfile: Path):
        '''Filters redundant sequence from `infile`'''

        cmd = ('gt', 'sequniq', '-seqit', '-force', '-o', outfile, infile)

        subprocess.check_call([str(arg) for arg in cmd])

    def sort_fasta_by_length(unsorted_fasta: Path, sorted_fasta: Path):
        '''Sorts a FASTA file by descending sequence length and writes it to
        a new FASTA file
        '''

        with unsorted_fasta.open('r') as fasta, sorted_fasta.open('w') as out:

            records = SeqIO.parse(fasta, 'fasta')

            sorted_records = sorted(records,
                                    key=lambda x: len(x.seq),
                                    reverse=True)

            SeqIO.write(sorted_records, out, 'fasta')

    fastq_to_fasta(fastq, fasta)

    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:

        temp_dir_path = Path(ephemeral)

        sorted_fasta = (temp_dir_path / fasta).with_suffix('.sorted')

        cropped = sorted_fasta.with_suffix('.cropped')

        uniq = cropped.with_suffix('.uniq')

        sort_fasta_by_length(fasta, sorted_fasta)

        crop_reads(sorted_fasta, cropped, 25, 50)

        sequniq(cropped, uniq)

        subseq(uniq, fasta, output)

def snap(subtracted: Path, workdir: Path, snap_db_dir: Path, tax_db_dir: Path,
         ribo_dir: Path, cores: int, edit_distance: int, cache_reset,
         comprehensive: bool, temp_dir: Path):

    basef = Path(subtracted.stem)
    annotated = workdir / basef.with_suffix('.annotated')

    matched = workdir / basef.with_suffix('.NT.snap.matched.sam')
    unmatched = workdir / basef.with_suffix('.NT.snap.unmatched.sam')

    fulllength_unmatched_fastq = unmatched.with_suffix('.fulllength.fastq')
    fulllength_matched_fastq = matched.with_suffix('.fulllength.fastq')

    cutadapt_fastq = workdir / basef.with_suffix('.cutadapt.fastq')

    snap_sam = snap_unmatched(subtracted, workdir, temp_dir,
                              snap_db_dir, edit_distance, cores)

    matched_lines, unmatched_lines = separate_sam_lines(snap_sam)

    write_sam(matched_lines, matched)
    write_sam(unmatched_lines, unmatched)

    with ProcessPoolExecutor(max_workers=2) as ppe1:

        ppe1.submit(extract_headers, cutadapt_fastq, matched, fulllength_matched_fastq)

        ppe1.submit(extract_headers, cutadapt_fastq, unmatched, fulllength_unmatched_fastq)

    viruses, bacteria, euks = lookup_taxonomy(fulllength_matched_fastq, matched,
                                              annotated, workdir, tax_db_dir)

    viruses_fastq = viruses.with_suffix('.fastq')

    try:
        ribo_snap(bacteria, 'BAC', cores, ribo_dir, temp_dir)
    except subprocess.CalledProcessError:
        user_msg('No bacterial reads to align')

    try:
        ribo_snap(euks, 'EUK', cores, ribo_dir, temp_dir)
    except subprocess.CalledProcessError:
        user_msg('No eukaryotic reads to align')

    # Long running, but low-memory, so parallel process
    with ProcessPoolExecutor(max_workers=4) as ppe2:
        for outtype in ('Accession', 'Genus', 'Species', 'Family'):
            ppe2.submit(table_generator, 'SNAP', outtype, viruses)

    if comprehensive:

        viruses_fastq.write_text(annotated_to_fastq(viruses))

        # output of this given to denovo assembly and RAPSearch
        extract_to_fast(fulllength_unmatched_fastq,
                        fulllength_unmatched_fastq.with_suffix('.fasta'),
                        fulllength_unmatched_fastq.with_suffix('.uniq.f1.fasta'),
                        temp_dir)

    return viruses, viruses_fastq, fulllength_unmatched_fastq.with_suffix('.uniq.f1.fasta')
