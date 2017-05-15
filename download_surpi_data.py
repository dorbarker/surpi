#!/usr/bin/env python3

import sys
import argparse
import requests
import subprocess
from datetime import date
from pathlib import Path

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('destination',
                        default='NCBI_{}'.format(date.today()),
                        type=Path,
                        help='Specify directory for downloaded data \
                              [NCBI_[current date]]')

    parser.add_argument('curated',
                        default='curated_{}'.format(date.today()),
                        type=Path,
                        help='Specify directory for curated data \
                              [curated_current date]]')

    parser.add_argument('--overwrite',
                        action='store_true',
                        help='Overwrite existing data [off]')

    return parser.parse_args()

def ensure_dir(d):

    if not Path(d).is_dir():
        Path(d).mkdir()

def download_file(src, dest, overwrite=False):

    if overwrite or not dest.exists():

        cmd = ('curl', '-o', str(dest), str(src))

        subprocess.call(cmd)

def download_ncbi(dest, overwrite):

    ensure_dir(dest)

    ncbi  = Path('ftp.ncbi.nih.gov')
    fasta = ncbi / 'blast/db/FASTA'
    taxny = ncbi / 'pub/taxonomy'

    download_file(fasta / 'nt.gz', dest / 'nt.gz', overwrite)
    download_file(fasta / 'nt.gz.md5', dest / 'nt.gz.md5', overwrite)

    download_file(fasta / 'nr.gz', dest / 'nr.gz', overwrite)
    download_file(fasta / 'nr.gz.md5', dest / 'nr.gz.md5', overwrite)

    download_file(taxny / 'taxdump.tar.gz', dest / 'taxdump.tar.gz', overwrite)
    download_file(taxny / 'taxdump.tar.gz.md5', dest / 'taxdump.tar.gz.md5', overwrite)

    for db in ('est', 'gb', 'gss', 'wgs'):

        nucl = 'nucl_{}.accession2taxid.gz'.format(db)
        nuclmd5 = str(nucl) + '.md5'
        url  = taxny / 'accession2taxid' / nucl
        urlmd5 = str(url) + '.md5'

        download_file(url, dest / nucl, overwrite)
        download_file(urlmd5, dest / nuclmd5, overwrite)

    prot = 'prot.accession2taxid.gz'

    download_file(taxny / 'accession2taxid' / prot , dest / prot, overwrite)
    download_file(taxny / 'accession2taxid' / prot + '.md5', dest / prot + '.md5', overwrite)

def download_curated(dest):

    ensure_dir(dest)

    chiu = Path('http://chiulab.ucsf.edu/SURPI/databases')

    download_list=('Bacterial_Refseq_05172012.CLEAN.LenFiltered.uniq.fa.gz',
                   'hg19_rRNA_mito_Hsapiens_rna.fa.gz',
                   'rapsearch_viral_aa_130628_db_v2.12.fasta.gz',
                   'viruses-5-2012_trimmedgi-MOD_addedgi.fa.gz',
                   '18s_rRNA_gene_not_partial.fa.gz' '23s.fa.gz',
                   '28s_rRNA_gene_NOT_partial_18s_spacer_5.8s.fa.gz',
                   'rdp_typed_iso_goodq_9210seqs.fa.gz')

    for dl in download_list:
        download_file(chiu / dl, dest / dl, overwrite=False)
        download_file(chiu / dl + '.md5', dest / dl + '.md5', overwrite=False)

def md5check(directory):

    md5s = Path(directory).glob('*.md5')

    for md5 in md5s:
        cmd = ('md5sum', '-c', '--status', md5)

        # if hashes don't match
        if subprocess.call(cmd):
            print('Checksum of {} failed'.format(md5), file=sys.stderr)

def main():

    args = arguments()

    download_test(Path('/home/dbarker/Desktop/test'), True)
    # NCBI data
    download_ncbi(args.destination, args.overwrite)

    md5check(args.destination)

    # Chiu lab data
    download_curated(args.curated)
    md5check(args.curated)

if __name__ == '__main__':
    main()
