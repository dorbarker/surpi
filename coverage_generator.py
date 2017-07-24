import requests
import tempfile
from pathlib import Path
from operator import itemgetter
import subprocess
from typing import List, Set
from Bio import SeqIO

from table_generator import create_tab_delimited_table, extract_barcodes
from utilities import fastq_to_fasta
from coveragePlot import coverage_plot

def coverage_generator(snap_file: Path, rap_file: Path, evalue: str,
                       cores: int, top_acc: int, top_coverage_plots: int):

    snap = snap_file.read_text().splitlines()
    rap = rap_file.read_text().splitlines()

    table = create_tab_delimited_table('SNAP', snap_file)

    barcodes = extract_barcodes(table)

    for barcode in barcodes:

        snap_pattern = '#{}/'.format(barcode)
        rap_pattern = '#{}'.format(barcode)

        snap_matches = [line for line in snap if snap_pattern in line]
        rap_matches = [line for line in rap if rap_pattern in line]

        barcode_genera = set(line[3] for line in table
                             if snap_pattern in line[0])

        for genus in barcode_genera:

            snap_mixed = [line for line in snap_matches if genus in line]

            snap_fastas = ['>{}\n{}'.format(line[0], line[9])
                           for line in snap_mixed]

            rap_fastas = ['>{}\n{}'.format(line[0], line[12])
                          for line in rap_matches if genus in line]

            snap_rap = snap_fastas + rap_fastas

            numreads = len(snap_rap)

            accs = [line[2] for line in snap_mixed][:200]

            fastas = get_genbank_fasta('fasta', accs)

            for acc in accs:

                plot_reads_to_acc(snap_rap, acc, genus, workdir, evalue, cores)

def get_genbank_fasta(rettype, accessions):
    '''Fetches requested accessions and returns them from NCBI'''

    acc_string = ' '.join(str(acc) for acc in accessions)

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
           db=sequences&\
           id={accs}&\
           rettype={rettype}&\
           retmode=text'.replace(' ','').format(accs=acc_string,
                                                rettype=rettype)

    res = requests.get(url)

    return res.text

def extract_accession_from_genbank(fasta_string):

    def extract_accession(header: str) -> str:

        acc = header.split()[0].lstrip('>')

        return acc[:acc.index('.')] if '.' in acc else acc

    headers = [line for line in fasta_string.splitlines()
               if line.startswith('>')]

    complete_genomes = (header for header in headers
                        if 'omplete' in header and 'enome' in header)

    accessions = [extract_accession(header) for header in complete_genomes]

    if not accessions:  # no complete genomes

        complete_sequences = (header for header in headers
                              if 'omplete' in header and 'equence' in header)

        accessions = [extract_accession(header)
                      for header in complete_sequences]

        if not accessions:  # no complete sequences

            accessions = [extract_accession(header) for header in headers]

def extract_to_fast(infile: Path, parent: Path, outfile: Path):

    header_file = infile.with_suffix('.headers')
    headers = [line.split()[0] for line in infile.read_text().splitlines()]

    header_file.write_text('\n'.join(headers))

    cmd = ('seqtk', 'subseq', str(parent), str(header_file))

    output = subprocess.check_output(cmd, universal_newlines=True)

    outfile.write_text(output)

    header_file.unlink()

def uniq_blast_results(infile: Path, outfile: Path) -> None:

    with infile.open('w') as inf:
        table = [line.split('\t') for line in inf]

    table.sort(key=itemgetter(10))
    table.sort(key=itemgetter(0))

    table2 = []
    done = set()  # type: Set[str]

    for line in table:
        if line[0] not in done:
            table2.append('\t'.join(line))
            done.add(line[0])

    outfile.write_text('\n'.join(table2))

def map_perfect_blast_to_genome(blast_out: Path, output: Path, genome_length):

    genome = []  # type: List[int]

    with blast_out.open('r') as blast, output.open('w') as out:
        for line in blast:

            data = line.split()

            oStart = int(data[6])
            oEnd = int(data[7])
            oLength = oEnd - oStart + 1
            gStart = min(int(data[8]), int(data[9]))
            gEnd = max(int(data[8]), int(data[9]))
            hitStart = max(0, (gStart - oStart + 1))
            hitEnd = min(genome_length, (gEnd + (oLength - oEnd)))

            for hit in range(hitStart, hitEnd + 1):
               genome[hit - 1] += 1

        for idx, val in enumerate(genome, 1):
            to_write = '{}\t{}\n'.format(idx, val)
            out.write(to_write)

    coverage_bp = len(genome) - genome.count(0)
    coverage_pc  = coverage_bp / genome_length
    depth = sum(genome) / genome_length
    return coverage_bp, coverage_pc, depth

def format_report(query, subject, length: int, coverage_bp: int, coverage_pc: float,
                  depth: float, num_reads: int, outpath: Path) -> None:

    out = '''Mapping: {}
    Against: {}
    Reference sequence length in bp: {}
    Coverage in bp: {}
    % Coverage: {}
    Average depth of coverage: {}
    Number of reads contributing to assembly: {}
    '''.format(query, subject, length, coverage_bp,
               coverage_pc, depth, num_reads)

    outpath.write_text(out)

def plot_reads_to_acc(snap_rap: Path, acc: str, genus: str, workdir: Path, evalue: str, cores: int):

    fasta_text = get_genbank_fasta('fasta', acc)

    fastafile = (workdir / acc).with_suffix('.fasta')
    blast_out = fastafile.with_suffix('.blastout')
    blast_uniq = blast_out.with_suffix('.uniq.blastout')
    blast_map = blast_out.with_suffix('.map')
    extracted = blast_out.with_suffix('.uniq.ex.fasta')

    report_path = snap_rap.parent / '.'.join([snap_rap.stem, acc, genus, 'report'])

    fastafile.write_text(fasta_text)

    makeblastdb = ['makeblastdb', '-in', str(fastafile), '-dbtype', 'nucl']

    blastn = ['blastn', '-task', 'blastn', '-outfmt', '8', '-evalue', evalue,
              '-num_threads', str(cores), '-num_alignments', '1', str(fastafile),
              '-out', str(blast_out)]

    subprocess.check_call(makeblastdb)
    subprocess.check_call(blastn)

    uniq_blast_results(blast_out, blast_uniq)
    extract_to_fast(blast_uniq, snap_rap, extracted)

    with extracted.open('r') as fasta:

        record = SeqIO.read(fasta, 'fasta')
        genome_length = len(record.seq)

    coverage_bp, coverage_pc, depth = map_perfect_blast_to_genome(blast_out,
                                                                  blast_map,
                                                                  genome_length)

    num_reads = len(blast_uniq.read_text().splitlines())

    format_report(snap_rap, acc, genome_length, coverage_bp, coverage_pc,
                  depth, num_reads, report_path)

    coverage_plot(str(blast_map), '{}_{}'.format(genus, acc))
