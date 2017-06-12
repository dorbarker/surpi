import subprocess
import re
from Bio import SeqIO
from utilities import concatenate

def run_rapsearch(query, output, database, cores, cutoff, fast, log):
    '''Runs RAPSearch with some default arguments overriden.'''

    a = 'T' if fast else 'F'

    rapsearch_cmd = ('rapsearch', '-q', query, '-o', output, '-d', database,
                     '-z', str(cores), '-e', cutoff, '-v', '1', '-b', '1',
                     '-t', 'N', '-a', a)

    rapsearch_stdout = subprocess.check_output(rapsearch_cmd,
                                               stderr=subprocess.STDOUT,
                                               universal_newlines=True)

    try:
        with open(log, 'w') as f:
            f.write(rapsearch_stdout)

    except TypeError:
        pass

def remove_octothorpe_lines(filepath):
    '''Deletes lines starting with the "#" character and
    writes back into the same file.
    '''

    with open(filepath, 'r') as i:
        lines = [l for l in i.readlines() if l.startswith('#')]

    with open(filepath, 'w') as o:
        o.writelines(lines)

def subseq(parent, query, output):
    '''Runs external process seqtk subseq.

    No return value, but creates a file specified by `output`
    '''

    subseq_cmd = ('seqtk', 'subseq', parent, query)

    result = subprocess.check_output(subseq_cmd, universal_newlines=True)

    with open(output, 'w') as f:
        f.write(result)

def fasta_to_seq(fastafile):
    '''Reads a fasta file, and returns a list of headerless sequences'''

    with open(fastafile, 'r') as f:

        recs = SeqIO.parse(f, 'fasta')

        seqs = [str(rec.seq) for rec in recs]

    return seqs

def paste(a_list, b_list):
    '''Pastes two [str] together, line-by-line,
    tab separated and newline terminated
    '''

    def merge(a, b):
        return '{}\t{}\n'.format(a.strip(), b.strip())

    return [merge(a, b) for a, b in zip(a_list, b_list)]

def filter_taxonomy(pattern, annotated):
    '''Filters the taxonomy annotation file based on a regular expression.

    Returns a list of matching patterns.
    '''
    p = re.compile(pattern)

    with open(annotated, 'r') as f:

        return [line for line in f if re.search(p, line)]

def rapsearch(rapsearch_query, rapsearch_output, database, cores,
              cutoff, fast, log):
    '''Functions shared by all RAPSearch strategies.'''

    m8 = rapsearch_output.with_suffix( rapsearch_output.suffix  + '.m8')
    m8fasta = m8.with_suffix( m8.suffix + '.fasta' )
    addseq = m8.with_suffix('.addseq.m8')

    run_rapsearch(rapsearch_query, rapsearch_output,
                  database, cores, cutoff, fast, log)


    remove_octothorpe_lines(m8)

    m8_lines = m8.read_text().splitlines()

    subseq(unmatched, m8, m8fasta)

    seqs = fasta_to_seq(m8fasta)

    pasted = paste(m8_lines, seqs)

    addseq.write_text('\n'.join(pasted))

    # TODO: taxonomy lookup goes here

def rapsearch_viral(query, vir_output, nr_output, vir_database, nr_database,
                    cores, vir_cutoff, nr_cutoff, fast, abyss_output,
                    contigs_nt_unmatched_fasta):

    # TODO: make contigs_nt_unmatched generated internally rather than
    # passing in by parameter

    m8 = output.with_suffix( output.suffix  + '.m8')
    log = output.with_suffix('.virlog')

    rapsearch(query, vir_output, vir_database, cores, cutoff, fast, log)

    # TODO: table_generator goes here

    concatenate(m8, abyss_output, contigs_nt_unmatched_fasta)
    assert contigs_nt_unmatched_fasta.exists(), "Contigs.NTunmatched.fasta \
                                                 doesn't exist"

    rapsearch(contigs_nt_unmatched_fasta, nr_output, nr_database, cores,
              nr_cutoff, fast, log.with_suffix('.nrlog'))

    viruses = filter_taxonomy('Viruses', nr_output)
    contigs = filter_taxonomy('^contig', nr_output)

    # TODO: table_generator goes here

    # TODO: coverage mapping

    # TODO: Get viral geaders no longer found in NR rapsearch

    # TODO: table_generator goes here


