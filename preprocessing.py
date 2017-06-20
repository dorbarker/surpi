'''Preprocesses fastq file for surpi, by running cutadapt,
trimming read lengths, and DUST masking sequence.
'''
import tempfile
import subprocess
from itertools import chain
from pathlib import Path

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

def cutadapt(infile, outfile, adapter_set, qual, quality_cutoff, length_cutoff):
    '''Runs cutadapt on an input FASTQ file'''

    def format_adapters(flags, adapters):
        '''Pairs each adapter with the cutadapt flags'''
        return tuple(chain(*zip(flags, adapters)))

    # Not including keep short reads because in the original SURPI script,
    # short reads are *never* kept

    adapter_info = infile.with_suffix('.adapterinfo')
    summary_log = infile.with_suffix('.summary.log')

    args = ('cutadapt', '-n', 15, '-O', 10, '-q', quality_cutoff,
            '-m', length_cutoff, '--quality-base={}'.format(qual),
            '-o', outfile, '--info-file={}'.format(adapter_info))

    if adapter_set == 'truseq':

        adapter_args = format_adapters(('-g', '-a', '-a', '-g', '-a', '-a'),
                                       ADAPTERS['truseq'])

    elif adapter_set == 'nextera':

        adapter_args = format_adapters(('-a', '-a', '-a'), ADAPTERS['nextera'])

    elif adapter_set == 'nexsoltruseq':

        flags = ('-g', '-a', '-a', '-g', '-a', '-a', '-a', '-a', '-a')
        adapter_args = format_adapters(flags, ADAPTERS['nexsoltruseq'])

        cmd = args + adapter_args + (infile, )

    elif adapter_set == 'nexsolb':

        adapter_args = format_adapters(('-g', '-a', '-a', '-a', '-a'),
                                       ADAPTERS['nexsolb'])

    else:
        # TODO: Error handling
        assert False

    cmd = args + adapter_args + (infile, )

    log = subprocess.check_output(map(str, cmd), universal_newlines=True)

    summary_log.write_text(log)

    # TODO: Confirm that outfile.cutadapt.fastq is created and the
    # sed 's/^$/N/g' call from the original SURPI is useless


def crop_reads(infile, outfile, start, length):
    '''Crops reads found at infile for `length` reads beginning
    at the 1-indexed position `start`

    Writes results to `outfile`, which may be the same Path as `infile`
    '''

    strt = start - 1
    end = strt + length

    with infile.open('r') as reads:

        lines = [line.strip() if idx % 2 is 1 else line.strip()[strt:end]
                 for idx, line in enumerate(reads.readlines(), 1)]

    outfile.write_text('\n'.join(lines))

def dust(infile, dusted):
    '''DUST mask infile and write it to `dusted`.fastq'''

    cmd = ('prinseq-lite.pl', '-fastq', infile, '-out_format', 3,
           '-out_good', dusted,
           '-out_bad', infile.with_suffix('.cutadapt.cropped.dusted.bad'),
           '-log', '-lc_method', 'dust', '-lc_threshold', 7)

    subprocess.check_call(map(str, cmd))

def preprocess(infile, workdir, tempdir, adapter_set, fastq_type,
               quality_cutoff, length_cutoff, crop_start, crop_length):
    '''Entry point for preprocessing.py

    Preprocesses FASTQ file for use in SURPI
    '''

    qual = 33 if fastq_type == 'S' else 64

    out_base = workdir / infile.stem

    cutadapt_output = out_base.with_suffix('.cutadapt.fastq')
    preproc = out_base.with_suffix('.preprocessed.fastq')

    # TODO: handle space removing
    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:

        throw_away = Path(ephemeral) / infile.stem
        cropped = throw_away.with_suffix('.cropped.fastq')
        dusted = throw_away.with_suffix('.dusted')

        cutadapt(infile, cutadapt_output, adapter_set, qual, quality_cutoff,
                 length_cutoff)

        crop_reads(cutadapt_output, cropped, crop_start, crop_length)

        dust(cropped, dusted)

        dusted.with_suffix('.dusted.fastq').rename(preproc)
