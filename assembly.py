from Bio import SeqIO
from multiprocessing import cpu_count
import argparse
import subprocess

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

def fastq_to_fasta(infile, outfile):

    with open(infile, 'r') as fastq, open(outfile, 'w') as fasta:

        for rec in SeqIO.parse(fastq, 'fastq'):
            rec.name = 'Vir{}'.format( rec.name )
            print(rec.name)
            SeqIO.write(rec, fasta, 'fasta')

def sequniq(matched, uniqvir, uniqunmatched, unmatched_add_vir):

    sequniq_cmd = ('gt', 'sequniq', '-seqit', '-force', '-o', uniqvir, matched)

    subprocess.call(sequniq_cmd)

    # cat uniqvir and uniqunmatched
    with open(uniqvir, 'r') as m, open(uniqunmatched, 'r') as u:

        with open(unmatched_add_vir, 'w') as o:
            o.write('{}\n{}'.format(m.read(), u.read()))

def abyss(unmatched_add_vir, length, cores, kmer, ignore_barcodes):

    contig_cutoff = int(1.75 * length)

    abyss_cmd = map(str, ('abyss_minimus.sh', unmatched_add_vir, length,
                          contig_cutoff, cores, kmer, ignore_barcodes))

    subprocess.call(abyss_cmd)

def process(matched_vir_fq, matched_vir_fa, matched_vir_uniq, uniqunmatched,
            addvir, length, cores, kmer, ignore_barcodes):

    sequniq(matched_vir_fa, matched_vir_uniq, uniqunmatched, addvir)

    abyss(addvir, length, cores, kmer, ignore_barcodes)

def main():

    args = arguments()

    process()

if __name__ == '__main__':
    main()
