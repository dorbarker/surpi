from Bio import SeqIO
from utilities import concatenate
import argparse
import subprocess

def arguments():

    parser = argparse.ArgumentParser()

    return parser.parse_args()

def rapsearch(query, output, database, cores, cutoff, fast, log):

    a = 'T' if fast else 'F'

    rapsearch_cmd = ('rapsearch', '-q', query, '-o', output, '-d', database,
                     '-z', str(cores), '-e', cutoff, '-v', '1', '-b', '1',
                     '-t', 'N', '-a', a)

    rapsearch_stdout = subprocess.check_output(rapsearch_cmd,
                                               stderr=subprocess.STDOUT,
                                               universal_newlines=TRUE)

    try:
        with open(log, 'w') as f:
            f.write(stdout)

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


