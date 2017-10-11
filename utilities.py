'''Utilities and shared functions for SURPI'''

from datetime import datetime
import sys
import functools
import fileinput
from itertools import chain
from pathlib import Path
import subprocess
import operator

from Bio import SeqIO

def user_msg(*args) -> None:
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

def run_shell(cmd, **kwargs) -> str:
    '''Runs a shell command using the subprocess module and returns the
    command's stdout.'''

    out = subprocess.check_output(cmd, shell=True,
                                  universal_newlines=True, **kwargs)
    return out.strip()

def concatenate(*files, output):
    '''Concatenate files together and write them to output'''

    files_ = [str(filename) for filename in files]
    with output.open('w') as out, fileinput.input(files_) as infiles:

        for line in infiles:
            out.write(line)

def annotated_to_fastq(annotated: Path, matches: bool) -> str:
    '''Converts and annotation Path to FASTQ format as a str.'''

    def fastq(line: str):
        '''Reformats a SAM line as FASTQ'''

        compare_match = operator.ne if matches else operator.eq
        split_line = line.split()

        replacements = {'fst': split_line[0],
                        'tenth': split_line[9],
                        'eleventh': split_line[10]}

        if compare_match(split_line[2], '*'):

            rec = '@{fst}\n{tenth}\n+{fst}\n{eleventh}'

            return rec.format(**replacements)

    lines = annotated.read_text().splitlines()

    fastq_lines = (fastq(line) for line in lines if not line.startswith('@'))

    return filter(None, fastq_lines)

def write_fastq_generators(filepath: Path, *fastqs):

    with filepath.open('w') as out:
        lines = ('{}\n'.format(line) for line in chain.from_iterable(fastqs))
        out.writelines(lines)

def fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> None:
    '''Reads a FASTQ file located at fastq_path and writes
    its FASTA conversion to fasta_path
    '''
    msg = 'Parameter {} must be of type pathlib.Path'
    assert isinstance(fastq_path, Path), msg.format('fastq_path')
    assert isinstance(fasta_path, Path), msg.format('fasta_path')


    with fastq_path.open('r') as fastq, fasta_path.open('w') as fasta:
        for record in SeqIO.parse(fastq, 'fastq'):
            SeqIO.write(record, fasta, 'fasta')

def sam_to_fasta(samfile: Path, fastafile: Path) -> None:
    '''Coverts a SAM-formatted file to FASTA format'''

    with samfile.open('r') as sam, fastafile.open('w') as fasta:

        for line in sam:

            split_line = line.strip().split()

            fasta_line = '>{}\n{}\n'.format(split_line[0], split_line[9])
            fasta.write(fasta_line)

def pigz(src: Path, dst: Path) -> None:
    '''Decompresses src with pigz, preserving the original file'''

    pigz_cmd = 'pigz -d -c -k {} > {}'.format(src, dst)

    ret = subprocess.run(pigz_cmd, shell=True, check=True)
