#!/usr/bin/env python3

import sys
import argparse
import subprocess
import hashlib
from math import ceil
from concurrent.futures import ProcessPoolExecutor
from datetime import date
from pathlib import Path
from typing import Tuple, Any
from Bio import SeqIO
from Bio.Seq import Seq
from utilities import user_msg, pigz
from create_taxonomy_db import create_taxonomy_database

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('ncbi',
                        default='NCBI_{}'.format(date.today()),
                        type=Path,
                        help='Specify directory for downloaded data \
                              [NCBI_[current date]]')

    parser.add_argument('curated',
                        default='curated_{}'.format(date.today()),
                        type=Path,
                        help='Specify directory for curated data \
                              [curated_current date]]')

    parser.add_argument('reference',
                        type=Path,
                        help='Path for SURPI reference data')

    parser.add_argument('--overwrite',
                        action='store_true',
                        help='Overwrite existing data [off]')

    return parser.parse_args()

def ensure_dir(directory: Path) -> None:

    if not Path(directory).is_dir():
        Path(directory).mkdir()

def download_file(src, dest, overwrite=False):

    if overwrite or not dest.exists():

        db = ('curl', '-o', str(dest), str(src))
        md5 = ('curl', '-o', str(dest) + '.md5', str(src) + '.md5')

        subprocess.call(md5)
        subprocess.call(db)

def download_ncbi(dest: Path, overwrite: bool):

    ensure_dir(dest)

    ncbi  = Path('ftp.ncbi.nih.gov')
    fasta = ncbi / 'blast/db/FASTA'
    taxny = ncbi / 'pub/taxonomy'

    download_file(fasta / 'nt.gz', dest / 'nt.gz', overwrite)


    download_file(fasta / 'nr.gz', dest / 'nr.gz', overwrite)


    download_file(taxny / 'taxdump.tar.gz', dest / 'taxdump.tar.gz', overwrite)


    for db in ('est', 'gb', 'gss', 'wgs'):

        nucl = 'nucl_{}.accession2taxid.gz'.format(db)

        url  = taxny / 'accession2taxid' / nucl

        download_file(url, dest / nucl, overwrite)


    prot = 'prot.accession2taxid.gz'

    download_file(taxny / 'accession2taxid' / prot , dest / prot, overwrite)


def download_curated(dest: Path):

    ensure_dir(dest)

    chiu = Path('chiulab.ucsf.edu/SURPI/databases')

    download_list=('Bacterial_Refseq_05172012.CLEAN.LenFiltered.uniq.fa.gz',
                   'hg19_rRNA_mito_Hsapiens_rna.fa.gz',
                   'rapsearch_viral_aa_130628_db_v2.12.fasta.gz',
                   'viruses-5-2012_trimmedgi-MOD_addedgi.fa.gz',
                   '18s_rRNA_gene_not_partial.fa.gz',
                   '23s.fa.gz',
                   '28s_rRNA_gene_NOT_partial_18s_spacer_5.8s.fa.gz',
                   'rdp_typed_iso_goodq_9210seqs.fa.gz')

    for dl in download_list:
        download_file(chiu / dl, dest / dl, overwrite=False)


def md5check(directory: Path):
    '''Checks that each downloaded file matches its expected MD5 sum'''

    def chunkfile(size, datafile):
        '''Lazily read a file and yield `size` bytes at a time'''

        while True:
            data = datafile.read(size)
            if not data:
                break
            yield data

    success = True
    md5s = directory.glob('*.md5')

    for md5 in md5s:

        expected = md5.read_text().strip().split()[0]

        with md5.with_suffix('').open('rb') as datafile:

            md5hash = hashlib.md5()

            for chunk in chunkfile(4096, datafile):
                md5hash.update(chunk)

            if not md5hash.hexdigest() == expected:
                user_msg('MD5 sum of {} failed'.format(md5))
                success = False

    # Exit only after all MD5 sums have been checked, so that the user
    # can see and fix all errors before re-running
    if not success:
        sys.exit(1)

def prerapsearch(database: Path, name: str, output: Path) -> None:

    cmd = ('prerapsearch', '-d', str(database), '-n', str(name))
    subprocess.check_call(cmd)

    info = Path(name).with_suffix('.info')
    (output / name).symlink_to(name)
    (output / info).symlink_to(info)

def organize_data(ncbi: Path, curated: Path, reference: Path) -> None:


    def snap_index(fasta: Path, dest_dir: Path, extra_args: Tuple[Any, ...]) -> None:

        tempout = fasta.with_suffix(fasta.suffix + '.snap_index')
        dest = dest_dir / tempout.name

        cmd = ['snap-aligner', 'index', fasta, tempout]

        cmd.extend(extra_args)

        subprocess.check_call([str(arg) for arg in cmd if arg])

        try:
            dest.symlink_to(tempout)
        except FileExistsError:
            dest.unlink()
            dest.symlink_to(tempout)

    def setup_reference_dirs(reference: Path):
        subdirs = [reference / sub for sub in ('taxonomy', 'rapsearch',
                                               'fast_snap', 'comp_snap',
                                               'riboclean_snap', 'host_snap')]

        for subdir in subdirs:
            if not subdir.exists():
                subdir.mkdir(parents=True)

        return subdirs

    def snap_index_nt(nt: Path, destdir: Path, prefix: str, n_chunks: int):

        # mask ambiguous positions with N, else SNAP fails
        tr = str.maketrans('WSMKRYBDHV', 'NNNNNNNNNN')

        output = nt.with_name('{}_{}'.format(prefix, nt.name))

        splitfasta = ('gt', 'splitfasta', '-force', '-numfiles', n_chunks, nt)

        subprocess.check_call([str(arg) for arg in splitfasta])

        for chunk in nt.parent.glob('*.[0-9]*'):

            dest = destdir / (prefix + nt.name + chunk.suffix)

      #      with chunk.open('r') as infile:
       #         records = list(SeqIO.parse(infile, 'fasta'))

#            with chunk.open('w') as outfile:

 #               for rec in records:
  #                  rec.description = ''
   #                 rec.id = rec.id[:rec.id.rfind('.')]
    #                rec.seq = Seq(str(rec.seq).translate(tr))
     #               SeqIO.write(rec, outfile, 'fasta')

            snap_index(chunk, destdir, ('-s', 22, '-locationSize', 5))

    today = '{}-{}-{}'.format(*date.today().timetuple()[:3])

    tax, rap, fast, comp, ribo, host = setup_reference_dirs(reference)

    nr = ncbi / 'nr'
    nt = ncbi / 'nt'

    if not nr.exists():
        pigz(nr.with_suffix('.gz'), nr)

    if not nt.exists():
        pigz(nt.with_suffix('.gz'), nt)

    # Create taxonomy databases
    #create_taxonomy_database(ncbi)

    #for tax_src in ncbi.glob('*.db'):
    #    tax_dst = reference / tax_src.name
    #    tax_dst.symlink_to(tax_src)

    # SNAP indices
    for gz_file in curated.glob('*.gz'):
        break  # remove
        decomp = gz_file.with_suffix('')

        if not decomp.exists():
            pigz(gz_file, decomp)

        if decomp.name.startswith('rapsearch'):

            viral = decomp  # processing happens below

        elif decomp.name.startswith('hg19'):

            snap_index(decomp, reference, ('-s', 22, '-hg19'))

        elif decomp.name.startswith('Bacterial_Refseq'):

            # Peak RAM usage ~47 GB
            snap_index(decomp, fast, ('-s', 20, '-locationSize', 5))

        elif decomp.name.startswith('viruses'):

            snap_index(decomp, fast, ('',))

        else:
            snap_index(decomp, ribo, ('',))

    # Run prerapsearch procs concurrently,
    # as they're long-running and relatively low-memory
    #with ProcessPoolExecutor(max_workers=2) as executor:

     #   executor.submit(prerapsearch, nr, 'rapsearch_nr', rap)
      #  executor.submit(prerapsearch, viral, 'rapsearch_viral', rap)

    snap_index_nt(nt, comp, today, 120)

def main():

    args = arguments()

    # NCBI data
    #download_ncbi(args.ncbi, args.overwrite)
    #md5check(args.ncbi)

    # Chiu lab data
    #download_curated(args.curated)
    #md5check(args.curated)

    organize_data(args.ncbi, args.curated, args.reference)

if __name__ == '__main__':
    main()
