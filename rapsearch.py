'''Manages running RAPSearch and processing RAPSearch output for SURPI'''

import subprocess
import re
import tempfile
from pathlib import Path
from typing import List

from Bio import SeqIO
from utilities import concatenate
from taxonomy_lookup import taxonomy_lookup
from table_generator import table_generator
def run_rapsearch(query: Path, output: Path, database: Path, cores: int,
                  cutoff: int, fast: bool, log: Path) -> None:
    '''Runs RAPSearch with some default arguments overriden.'''

    fast_mode = 'T' if fast else 'F'

    rapsearch_cmd = ('rapsearch', '-q', query, '-o', output, '-d', database,
                     '-z', cores, '-e', cutoff, '-v', '1', '-b', '1',
                     '-t', 'N', '-a', fast_mode)

    rapsearch_stdout = subprocess.check_output([str(x) for x in rapsearch_cmd],
                                               stderr=subprocess.STDOUT,
                                               universal_newlines=True)

    try:
        log.write_text(rapsearch_stdout)

    except TypeError:
        pass

def remove_octothorpe_lines(filepath: Path) -> None:
    '''Deletes lines starting with the "#" character and
    writes back into the same file.
    '''

    with filepath.open('r') as infile:
        lines = [l for l in infile.readlines() if l.startswith('#')]

    filepath.write_text('\n'.join(lines))

def subseq(parent: Path, query: Path, output: Path) -> None:
    '''Runs external process seqtk subseq.

    No return value, but creates a file specified by `output`
    '''

    subseq_cmd = [str(x) for x in ('seqtk', 'subseq', parent, query)]

    result = subprocess.check_output(subseq_cmd, universal_newlines=True)

    output.write_text(result)

def fasta_to_seq(fastafile: Path) -> List[str]:
    '''Reads a fasta file, and returns a list of headerless sequences'''

    with fastafile.open('r') as fasta:

        recs = SeqIO.parse(fasta, 'fasta')

        seqs = [str(rec.seq) for rec in recs]

    return seqs

def paste(a_list: List[str], b_list: List[str]) -> List[str]:
    '''Pastes two [str] together, line-by-line,
    tab separated and newline terminated
    '''

    def merge(a_elem: str, b_elem: str):
        '''Merges paired line elements'''

        return '{}\t{}\n'.format(a_elem.strip(), b_elem.strip())

    return [merge(a_elem, b_elem) for a_elem, b_elem in zip(a_list, b_list)]

def filter_taxonomy(pattern: str, annotated: Path) -> List[str]:
    '''Filters the taxonomy annotation file based on a regular expression.

    Returns a list of matching patterns.
    '''
    comppat = re.compile(pattern)

    with annotated.open('r') as annot:

        return [line for line in annot if re.search(comppat, line)]

def coverage_map(snap_matched: Path, annotated: Path,
                 evalue: str, cores: int) -> None:
    '''Removes "@" from taxonomy annotation
    and runs coverage_generator_bp.sh
    '''

    name = annotated.stem.split('.')[0]

    annot = annotated.read_text().replace('@', '')

    with tempfile.NamedTemporaryFile() as inc_bar:
        inc_bar.write(annot)
        inc_bar.seek(0)
        inc_bar_name = Path(inc_bar.name)

        cov_gen = [str(x) for x in ('coverage_generator_bp.sh', snap_matched,
                                    inc_bar_name, evalue, cores, 10, 1, name)]

        subprocess.check_call(cov_gen)

    inc_bar_name.unlink()

def rapsearch_shared(rapsearch_query: Path, rapsearch_output: Path,
                     rap_database: Path, tax_db_dir: Path, cores: int,
                     cutoff: int, fast: bool, log: Path) -> None:
    '''Functions shared by all RAPSearch strategies.'''

    annotated = rapsearch_output.with_suffix('.annotated')
    addseq = rapsearch(rapsearch_query, rapsearch_output, rap_database,
                       cores, cutoff, fast, log)

    taxonomy_lookup(infile=addseq, outfile=annotated, db_dir=tax_db_dir,
                    seq_type='nucl', file_type='blast')

def rapsearch(query: Path, output: Path, database: Path, cores: int,
              cutoff: int, fast: bool, log: Path) -> Path:
    '''Run RAPSearch and handle immediate processing of its output'''

    m8 = output.with_suffix(output.suffix + '.m8')
    m8fasta = m8.with_suffix(m8.suffix + '.fasta')
    addseq = m8.with_suffix('.addseq.m8')


    run_rapsearch(query, output, database, cores, cutoff, fast, log)

    remove_octothorpe_lines(m8)

    m8_lines = m8.read_text().splitlines()

    subseq(query, m8, m8fasta)

    seqs = fasta_to_seq(m8fasta)

    pasted = paste(m8_lines, seqs)

    addseq.write_text('\n'.join(pasted))

    return addseq

def vir_nr_difference(vir_annot: List[str], nr_annot: List[str]) -> List[str]:
    '''Returns the viral annotations not found in NR'''

    def get_headers(annot: List[str]) -> List[str]:
        '''Extracts the headers from a taxomomy annotation'''
        return [line.split('\t')[0] for line in annot]

    vir_headers = get_headers(vir_annot)
    nr_headers = get_headers(nr_annot)

    headers_not_in_nr = set(vir_headers) - set(nr_headers)

    not_in_nr = [annot for header, annot in zip(vir_headers, vir_annot)
                 if header in headers_not_in_nr]

    return not_in_nr

def rapsearch_viral(query: Path, workdir: Path, vir_database: Path,
                    nr_database: Path, vir_cutoff: int, nr_cutoff: int,
                    fast: bool, abyss_output: Path, tax_db_dir: Path,
                    snap_match_annot: Path, evalue: str, cores: int) -> None:
    '''Performs a RAPSearch first against a specialized viral database,
    and then against the NCBI NR database.
    '''

    vir_output = (workdir / query.stem).with_suffix('.rapsearch.vir')
    nr_output = (workdir / query.stem).with_suffix('.rapsearch.nr')

    m8 = vir_output.with_suffix(vir_output.suffix  + '.m8')
    log = vir_output.with_suffix('.virlog')
    virus_tax = nr_output.with_suffix('.viruses.annotated')
    contig_tax = nr_output.with_suffix('.NRcontigs.annotated')
    contigs_nt_unmatched_fasta = m8.with_suffix('Contigs.NTunmatched.fasta')
    not_in_nr_annot = m8.with_suffix('.notNR.annotated')

    rapsearch_shared(query, vir_output, vir_database, tax_db_dir,
                     cores, vir_cutoff, fast, log)

    table_generator('RAP', 'Genus', vir_output.with_suffix('.annotated'))

    concatenate(m8, abyss_output, output=contigs_nt_unmatched_fasta)

    assert contigs_nt_unmatched_fasta.exists(), "Contigs.NTunmatched.fasta \
                                                 doesn't exist"

    rapsearch_shared(contigs_nt_unmatched_fasta, nr_output, nr_database,
                     tax_db_dir, cores, nr_cutoff, fast,
                     log.with_suffix('.nrlog'))

    viruses = filter_taxonomy('Viruses', nr_output)
    virus_tax.write_text('\n'.join(viruses))


    contigs = filter_taxonomy('^contig', nr_output)
    contig_tax.write_text('\n'.join(contigs))

    for outtype in ('Accession', 'Species', 'Genus', 'Family'):

        table_generator('RAP', outtype, contig_tax)

        table_generator('RAP', outtype, virus_tax)

    coverage_map(snap_match_annot, virus_tax, evalue, cores)

    not_in_nr = vir_nr_difference(viruses, contigs)
    not_in_nr_annot.write_text('\n'.join(not_in_nr))

    table_generator('RAP', 'Genus', not_in_nr_annot)

def rapsearch_nr(snap_unmatched: Path, workdir: Path, abyss_output: Path,
                 rap_database: Path, tax_db_dir: Path, cutoff: int, fast: bool,
                 cores: int) -> None:
    '''Performs a RAPSearch against the NCBI NR database'''

    output = (workdir / snap_unmatched).with_suffix('.rapsearch.nr')
    contigs_nt = abyss_output.parent / 'contigs_nt_snap_unmatched.fasta'
    concatenate(snap_unmatched, abyss_output, output=contigs_nt)
    log = output.with_suffix('.nrlog')
    virus_tax = output.with_suffix('.viruses.annotated')
    novir_tax = output.with_suffix('.novir.annotated')

    rapsearch_shared(contigs_nt, output, rap_database, tax_db_dir, cores,
                     cutoff, fast, log)

    virus = filter_taxonomy('Viruses', output)
    novir = [tax for tax in filter_taxonomy('^contig', output)
             if 'Viruses' not in tax]

    virus_tax.write_text('\n'.join(virus))
    novir_tax.write_text('\n'.join(novir))


    for outtype in ('Accession', 'Species', 'Genus', 'Family'):

        table_generator('RAP', outtype, virus_tax)

    table_generator('RAP', 'Genus', novir_tax)
