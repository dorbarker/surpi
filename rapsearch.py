'''Manages running RAPSearch and processing RAPSearch output for SURPI'''

import subprocess
import re
import tempfile
from Bio import SeqIO
from utilities import concatenate
from taxonomy_lookup import taxonomy_lookup
from pathlib import Path

def run_rapsearch(query, output, database, cores, cutoff, fast, log):
    '''Runs RAPSearch with some default arguments overriden.'''

    fast_mode = 'T' if fast else 'F'

    rapsearch_cmd = ('rapsearch', '-q', query, '-o', output, '-d', database,
                     '-z', str(cores), '-e', cutoff, '-v', '1', '-b', '1',
                     '-t', 'N', '-a', fast_mode)

    rapsearch_stdout = subprocess.check_output(rapsearch_cmd,
                                               stderr=subprocess.STDOUT,
                                               universal_newlines=True)

    try:
        with open(log, 'w') as logfile:
            logfile.write(rapsearch_stdout)

    except TypeError:
        pass

def remove_octothorpe_lines(filepath):
    '''Deletes lines starting with the "#" character and
    writes back into the same file.
    '''

    with open(filepath, 'r') as infile:
        lines = [l for l in infile.readlines() if l.startswith('#')]

    with open(filepath, 'w') as outfile:
        outfile.writelines(lines)

def subseq(parent, query, output):
    '''Runs external process seqtk subseq.

    No return value, but creates a file specified by `output`
    '''

    subseq_cmd = ('seqtk', 'subseq', parent, query)

    result = subprocess.check_output(subseq_cmd, universal_newlines=True)

    with open(output, 'w') as outfile:
        outfile.write(result)

def fasta_to_seq(fastafile):
    '''Reads a fasta file, and returns a list of headerless sequences'''

    with open(fastafile, 'r') as fasta:

        recs = SeqIO.parse(fasta, 'fasta')

        seqs = [str(rec.seq) for rec in recs]

    return seqs

def paste(a_list, b_list):
    '''Pastes two [str] together, line-by-line,
    tab separated and newline terminated
    '''

    def merge(a_elem, b_elem):
        '''Merges paired line elements'''

        return '{}\t{}\n'.format(a_elem.strip(), b_elem.strip())

    return [merge(a_elem, b_elem) for a_elem, b_elem in zip(a_list, b_list)]

def filter_taxonomy(pattern, annotated):
    '''Filters the taxonomy annotation file based on a regular expression.

    Returns a list of matching patterns.
    '''
    comppat = re.compile(pattern)

    with open(annotated, 'r') as annot:

        return [line for line in annot if re.search(comppat, line)]

def table_generator(annotated, snap_rap, acc, species, genus, family):
    '''Runs external script `table_generator.sh`'''

    assert snap_rap in ('SNAP', 'RAP')
    assert acc in ('Y', 'N')
    assert species in ('Y', 'N')
    assert genus in ('Y', 'N')
    assert family in ('Y', 'N')

    cmd = ('table_generator.sh', str(annotated),
           snap_rap, acc, species, genus, family)

    subprocess.check_call(cmd)

def coverage_map(snap_matched, annotated, evalue, cores):
    '''Removes "@" from taxonomy annotation
    and runs coverage_generator_bp.sh
    '''

    name = annotated.stem.split('.')[0]

    annot = annotated.read_text().replace('@', '')

    with tempfile.NamedTemporaryFile() as inc_bar:
        inc_bar.write(annot)
        inc_bar.seek(0)
        inc_bar_name = Path(inc_bar.name)

        cov_gen = ('coverage_generator_bp.sh', snap_matched, inc_bar_name,
                   evalue, cores, 10, 1, name)

        subprocess.check_call(map(str, cov_gen))

    inc_bar_name.unlink()

def rapsearch_shared(rapsearch_query, rapsearch_output, rap_database,
                     tax_db_dir, cores, cutoff, fast, log):
    '''Functions shared by all RAPSearch strategies.'''

    annotated = rapsearch_output.with_suffix('.annotated')
    addseq = rapsearch(rapsearch_query, rapsearch_output, rap_database,
                       cores, cutoff, fast, log)

    taxonomy_lookup(infile=addseq, outfile=annotated, db_dir=tax_db_dir,
                    seq_type='nucl', file_type='blast')

def rapsearch(query, output, database, cores, cutoff, fast, log):
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

def rapsearch_viral(query, vir_output, nr_output, vir_database, nr_database,
                    cores, vir_cutoff, nr_cutoff, fast, abyss_output,
                    contigs_nt_unmatched_fasta, tax_db_dir, snap_match_annot,
                    evalue):

    # TODO: make contigs_nt_unmatched generated internally rather than
    # passing in by parameter

    m8 = vir_output.with_suffix(vir_output.suffix  + '.m8')
    log = vir_output.with_suffix('.virlog')
    virus_tax = nr_output.with_suffix('.viruses.annotated')
    contig_tax = nr_output.with_suffix('.NRcontigs.annotated')

    rapsearch_shared(query, vir_output, vir_database, tax_db_dir,
                     cores, vir_cutoff, fast, log)

    table_generator(vir_output.with_suffix('.annotated'),
                    'RAP', 'N', 'Y', 'N', 'N')

    concatenate(m8, abyss_output, output=contigs_nt_unmatched_fasta)

    assert contigs_nt_unmatched_fasta.exists(), "Contigs.NTunmatched.fasta \
                                                 doesn't exist"

    rapsearch_shared(contigs_nt_unmatched_fasta, nr_output, nr_database,
                     tax_db_dir, cores, nr_cutoff, fast,
                     log.with_suffix('.nrlog'))

    viruses = filter_taxonomy('Viruses', nr_output)
    virus_tax.write_text('\n'.join(viruses))

    table_generator(virus_tax, 'RAP', 'Y', 'Y', 'Y', 'Y')

    contigs = filter_taxonomy('^contig', nr_output)
    contig_tax.write_text('\n'.join(contigs))

    table_generator(contig_tax, 'RAP', 'Y', 'Y', 'Y', 'Y')


    coverage_map(snap_match_annot, virus_tax, evalue, cores)

    # TODO: Get viral geaders no longer found in NR rapsearch

    # TODO: table_generator goes here

def rapsearch_nr(snap_unmatched, abyss_output, output, rap_database,
                 tax_db_dir, cores, cutoff, fast, log):

    contigs_nt = abyss_output.parent / 'contigs_nt_snap_unmatched.fasta'
    concatenate(snap_unmatched, abyss_output, output=contigs_nt)
    log = output.with_suffix('.log')
    virus_tax = output.with_suffix('.viruses.annotated')
    novir_tax = output.with_suffix('.novir.annotated')

    rapsearch_shared(contigs_nt, output, rap_database, tax_db_dir, cores,
                     cutoff, fast, log)

    virus = filter_taxonomy('Viruses', output)
    novir = [tax for tax in filter_taxonomy('^contig', output)
               if 'Viruses' not in tax]

    virus_tax.write_text('\n'.join(virus))
    novir_tax.write_text('\n'.join(novir))

    table_generator(virus_tax, 'RAP', 'Y', 'Y', 'Y', 'Y')
    table_generator(novir_tax, 'RAP', 'N', 'Y', 'N', 'N')

    # TODO: coverage map

