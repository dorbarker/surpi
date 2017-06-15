'''Utilities and shared functions for SURPI'''

from datetime import datetime
import sys
import functools
import fileinput
import pathlib
import subprocess

from Bio import SeqIO

def user_msg(*args):
    '''Wrapper for print() that prints to stderr'''

    print(*args, file=sys.stderr)

def logtime(name):
    '''Function decorator that print to stderr the runtime of the
    decorated function.
    '''

    def decorator(func):
        '''Interface between wrapper and the outer logtime() function'''

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            '''Wraps func and prints to STDERR the runtime of func'''

            msg = 'Elapsed time for {}: {}'
            before = datetime.now()

            result = func(*args, *kwargs)

            after = datetime.now()

            user_msg(msg.format(name, after - before))
            return result

        return wrapper
    return decorator

def run_shell(cmd, **kwargs):
    '''Runs a shell command using the subprocess module and returns the
    command's stdout.'''

    out = subprocess.check_output(cmd, shell=True,
                                  universal_newlines=True, **kwargs)
    return out.strip()

def concatenate(*files, output):
    '''Concatenate files together and write them to output'''

    with open(output, 'w') as out, fileinput.input(files) as infiles:

        for line in infiles:
            out.write(line)

def annotated_to_fastq(annotated):
    '''Converts and annotation Path to FASTQ format as a str.'''

    def fastq(line):
        '''Reformats a SAM line as FASTQ'''

        split_line = line.split()

        replacements = {'fst': split_line[0],
                        'tenth': split_line[9],
                        'eleventh': split_line[10]}

        if line[2] == '*':

            rec = '@{fst}\n{tenth}\n+{fst}\n{eleventh}\n'

            return rec.format(**replacements)

    lines = annotated.read_text().splitlines()

    fastq_lines = (fastq(line) for line in lines if not line.startswith('@'))

    return '\n'.join(fastq_lines)

def fastq_to_fasta(fastq_path, fasta_path):
    '''Reads a FASTQ file located at fastq_path and writes
    its FASTA conversion to fasta_path
    '''
    msg = 'Parameter {} must be of type pathlib.Path'
    assert isinstance(fastq_path, pathlib.Path), msg.format('fastq_path')
    assert isinstance(fasta_path, pathlib.Path), msg.format('fasta_path')


    with fastq_path.open('r') as fastq, fasta_path.open('w') as fasta:
        for record in SeqIO.parse(fastq, 'fastq'):
            SeqIO.write(record, fasta, 'fasta')
