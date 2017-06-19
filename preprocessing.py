import tempfile
import subprocess
from itertools import chain
import re

ADAPTERS = {
    'truseq': ('GTTTCCCACTGGAGGATA',
               'TATCCTCCAGTGGGAAAC',
               'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
               'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC',
               'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
               'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC'),

    'nextera': ('CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
                'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
                'CTGTCTCTTATACACATCT')
}

ADAPTERS['nexsoltruseq'] = ADAPTERS['truseq'] + ADAPTERS['nextera']
ADAPTERS['nexsolb'] = ADAPTERS['truseq'][:2] + ADAPTERS['nextera']

def cutadapt(infile, outfile, adapter_set, qual, keep_short, quality_cutoff,
             length_cutoff):
    '''Runs cutadapt on an input FASTQ file'''

    def format_adapters(flags, adapters):
        '''Pairs each adapter with the cutadapt flags'''
        return tuple(chain(*zip(flags, adapters)))

    def post_proc(filepath):

        filepath.read_text()

    # Not including keep short reads because in the original SURPI script,
    # short reads are *never* kept

    adapter_info = infile.with_suffix('.adapterinfo')
    summary_log = infile.with_suffix('.summary.log')

    args = ('cutadapt', '-n', 15, '-O', 10, '-q', quality_cutoff,
            '-m', length_cutoff, '--quality-base={}'.format(qual),
            '-o', outfile, '--info-file={}'.format(adapter_info))

    if adapter_set == 'truseq':

        adapter_args = format_adapters(('-g', '-a', '-a', '-g', '-a' '-a'),
                                       ADAPTERS['truseq'])

        cmd = args + adapter_args + (infile,)

    elif adapter_set == 'nextera':

        adapter_args = format_adapters(('-a', '-a', '-a'), ADAPTERS['nextera'])

        cmd = args + adapter_args + (infile, )

    elif adapter_set == 'nexsoltruseq':

        flags = ('-g', '-a', '-a', '-g', '-a', '-a', '-a', '-a', '-a')
        adapter_args = format_adapters(flags, ADAPTERS['nexsoltruseq'])

        cmd = args + adapter_args + (infile, )

    elif adapter_set == 'nexsolb':

        adapter_args = format_adapters(('-g', '-a', '-a', '-a', '-a'),
                                       ADAPTERS['nexsolb'])

        cmd = args + adapter_args + (infile, )

    else:
        # TODO: Error handling
        '''Something has gone awry'''

    log = subprocess.check_output(map(str, cmd), universal_newlines=True)

    summary_log.write_text(log)

    # TODO: Confirm that outfile.cutadapt.fastq is created and the
    # sed 's/^$/N/g' call from the original SURPI is useless



def crop():
    pass

def dust():
    pass

def preprocess(infile, outfile, tempdir):

    # TODO: confirm fastq filter step is never run
    # cutadapt
    # cutadapt postprocessing
    pass
