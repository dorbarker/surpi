from multiprocessing import cpu_count
import argparse
import subprocess
from utilities import concatenate, fastq_to_fasta

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-k', '--kmer', type=int, default=34,
                        help='k-mer size for ABySS [34]')

    # TODO handle this automatically
    parser.add_argument('--length', type=int, help='Length of line 1 in FASTQ')

    parser.add_argument('--ignore-barcodes',
                        choices=('Y', 'N'),
                        default='N',
                        help="If 'Y', then assemble all barcodes together in \
                              a single assembly, otherwise assemble each \
                              barcode separately [N]")

    parser.add_argument('--cores', default=cpu_count(), type=int,
                        help='CPU cores to use [all]')

    parser.add_argument('matched', help='Matched viral FASTQ')

    parser.add_argument('unmatched', help='Unmatched FASTQ')

    parser.add_argument('unmatched-add-vir',
                        help='Path to create FASTA of unmatched reads \
                              combined with viral reads')

    return parser.parse_args()

def sequniq(matched, uniqvir, uniqunmatched, unmatched_add_vir):
    '''Extracts the unique sequences from `matched`, concatenates them with
    the unique unmatched sequences, and writes it all to `unmatched_add_vir`
    '''

    sequniq_cmd = ('gt', 'sequniq', '-seqit', '-force', '-o', uniqvir, matched)

    subprocess.check_call(sequniq_cmd)

    concatenate(uniqvir, uniqunmatched, output=unmatched_add_vir)

def abyss(unmatched_add_vir, length, cores, kmer, ignore_barcodes):
    '''Executes abyss_minimus.sh'''

    contig_cutoff = int(1.75 * length)

    abyss_cmd = map(str, ('abyss_minimus.sh', unmatched_add_vir, length,
                          contig_cutoff, cores, kmer, ignore_barcodes))

    subprocess.call(abyss_cmd)

def assemble(matched_vir_fastq, matched_vir_fasta, matched_vir_uniq,
             uniqunmatched, addvir, sample_fastq_length, cores, kmer,
             ignore_barcodes):

    fastq_to_fasta(matched_vir_fastq, matched_vir_fasta)

    sequniq(matched_vir_fasta, matched_vir_uniq, uniqunmatched, addvir)

    abyss(addvir, sample_fastq_length, cores, kmer, ignore_barcodes)

def main():

    args = arguments()


if __name__ == '__main__':
    main()
