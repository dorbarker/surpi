#!/usr/bin/env python3

# system imports
import glob
import os
import argparse
import subprocess
from multiprocessing import cpu_count
import sys
from pathlib import Path
from typing import Dict

# surpi imports
from preprocessing import preprocess
from utilities import run_shell, user_msg, logtime
from snap_to_nt import snap
from map_host import host_subtract
from assembly import assemble
from rapsearch import rapsearch_viral, rapsearch_nr

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--input-type',
                        choices=('fasta', 'fastq'),
                        default='fastq',
                        help='Input filetype [fastq]')

    parser.add_argument('--quality-mode',
                        choices=('illumina', 'sanger'),
                        default='sanger',
                        help='Quality score mode [sanger]')

    parser.add_argument('--adapter-set',
                        choices=('truseq', 'nextera',
                                 'nexsolb', 'nexsoltruseq'),
                        default='truseq',
                        help='Adapter set to use [truseq]')

    # TODO: fix formatting
    parser.add_argument('--verify-fastq',
                        choices=(0, 1, 2, 3),
                        default=1,
                        help='''
                        0 = skip fastq validation
                        1 = validate fastqs; quit on fail
                        2 = validate fastqs; check for unique names; quit on fail
                        3 = validate fastqs; check for unique names; do not quit on failure
                        [1]
                        ''')

    parser.add_argument('--run-mode',
                        choices=('comprehensive', 'fast'),
                        default='comprehensive',
                        help='''
                        Comprehensive mode allows SNAP to NT -> denovo contig assembly -> RAPSearch to Viral proteins or NR
                        Fast mode allows SNAP to curated FAST databases
                        [comprehensive]
                        ''')

    parser.add_argument('--skip-preprocess',
                        action='store_true',
                        help='Skip preprocessing step [off]')

    parser.add_argument('--length-cutoff',
                        type=int,
                        default=50,
                        help='Discard any post-trimming sequence less than this length [50]')

    parser.add_argument('--crop-start',
                        type=int,
                        default=10,
                        help='Prior to SNAP alignment, cropping start position [10]')

    parser.add_argument('--crop-length',
                        type=int,
                        default=10,
                        help='Prior to SNP alignment, length to crop after start position [75]')

    parser.add_argument('--quality-cutoff',
                        type=int,
                        help='Quality cutoff passed to cutadapt\' -q parameter [18]')

    parser.add_argument('--edit-distance',
                        type=int,
                        default=12,
                        help='SNAP edit distance for computational subtraction of host genome [12]')


    # snap_nt interator option only has one choice: 'inline' and so is omitted here

    parser.add_argument('--rapsearch-method',
                        choices=('viral', 'nr'),
                        default='viral',
                        help='RAPSearch database method to use [viral]')

    parser.add_argument('--rapsearch-cutoff',
                        default='1',
                        help='E-value cutoff for RAPSearch. Will parse \
                              notation such as "1e+1" [1]')

    parser.add_argument('--fast-rapsearch',
                        action='store_true',
                        help='Yields a 10-30x speed increase at the cost of sensitivity [off]')

    parser.add_argument('--abyss-kmer',
                        type=int,
                        default=34,
                        help='kmer size for ABySS de novo assembly [34]')

    parser.add_argument('--ignore-barcodes',
                        action='store_true',
                        help='Assemble all barcodes together into a single assembly [off]')

    parser.add_argument('--blastn-evalue',
                        default='1e-15',
                        help='E-value cutoff for BLASTn/ Will parse notation \
                             such as "1e-15" [1e-15]')

    # Reference Data

    parser.add_argument('--snap-subtract-dir',
                        default='reference/hg19',
                        help='Directory containing the SNAP-indexed \
                             database of the host genome for the \
                             subtraction phase. SURPI will subtract \
                             all SNAP databases found in this directory. \
                             [reference/hg19]')

    parser.add_argument('--snap-comprehensive-dir',
                        default='reference/COMP_SNAP',
                        help='Directory for SNAP-indexed databases of \
                             NCBI NT. Directory must contain *only* SNAP \
                             indexed databases [reference/COMP_SNAP]')

    parser.add_argument('--snap-fast-dir',
                        default='reference/FAST_SNAP',
                        help='Directory for FAST-mode SNAP databases \
                             [reference/FAST_SNAP]')

    # update if create_taxonomy_db.sh changes
    parser.add_argument('--tax-ref-dir',
                        default='reference/taxonomy',
                        help='Directory containing 3 SQLite databases \
                             created by "create_taxonomy_db.sh" \
                             [reference/taxonomy]')

    parser.add_argument('--rapsearch-viral',
                        default='reference/RAPSearch/rapsearch_viral_db',
                        help='Directory containing the RAPSearch \
                             viral database \
                             [reference/RAPSearch/rapsearch_viral_db]')

    parser.add_argument('--rapsearch-nr',
                        default='reference/RAPSearch/rapsearch_nr_db',
                        help='Directory containing the RAPSearch NCBI NR \
                             database [reference/RAPSearch/rapsearch_nr_db')

    parser.add_argument('--snap-riboclean',
                        default='reference/RiboClean_SNAP',
                        help='SNAP RiboClean directory \
                             [reference/RiboClean_SNAP]')

    parser.add_argument('--cores',
                        type=int,
                        default=cpu_count(),
                        help='CPU cores to use [{}]'.format(cpu_count()))

    parser.add_argument('--temp',
                        default='/tmp/',
                        help='Temporary working directory [/tmp/]')

    parser.add_argument('--cache-reset',
                        type=int,
                        default=0,
                        help='dropcache if free RAM (GB) falls below this \
                             value. If 0, then the cache is never reset [0]')

    return parser.parse_args()

def ensure_fastq(inputfile: Path, workdir: Path, mode: str) -> Path:
    '''If mode is "fasta", convert inputfile to fastq.

    In either case, guarantee that file extension is ".fastq"
    '''

    assert mode in ('fastq', 'fasta'), 'Mode must be "fasta" or "fastq"'

    fq_name = workdir / inputfile.with_suffix('.fastq').name

    if mode == 'fasta':

        cmd = ('fasta_to_fastq', str(inputfile))

        fq_name.write_text(subprocess.check_output(cmd,
                                                   universal_newlines=True))

    else:

        fq_name.symlink_to(inputfile)

    return fq_name

def free_memory(cache_reset: int):

    freemem = float(run_shell("free -g | awk 'NR==2{print $4}'"))

    if freemem < cache_reset:
        run_shell('dropcache')

def get_memory_limit() -> int:
    '''Return a usable memory limit in GB based on system RAM'''

    cmd = ('grep', 'MemTotal', '/proc/meminfo')

    meminfo = subprocess.check_output(cmd, universal_newlines=True)

    total_memory = int(meminfo.split()[1])

    if total_memory >= 500000000:

        memory_limit = 200

    elif total_memory >= 200000000:

        memory_limit = 150

    else:

        memory_limit = 50

    return memory_limit

def check_program_versions() -> Dict[str, str]:
    '''Return the versions of each external dependency as a dict

    {program: version}
    '''

    version_cmds = {
        'gt':        "gt -version | head -1 | awk '{print $3}'",
        'seqtk':     "seqtk 2>&1 | head -3 | tail -1 | awk '{print $2}'",
        'cutadapt':  "cutadapt --version",
        'prinseqlite':  "prinseq-lite.pl --version 2>&1 | awk '{print $2}'",
        'snap-aligner': "snap 2>&1 | grep version | awk '{print $5}'",
        'rapsearch': "rapsearch 2>&1 | head -2 | tail -1 | awk '{print $2}'",
        'abyss_pe':  "abyss-pe version | head -2 | tail -1 | awk '{print $3}'",
        'abyss_p':   "ABYSS-P  --version | head -1 | awk '{print $3}'",
        'minimo':    "Minimo -h | tail -2 | awk '{print $2}'"
        }

    versions = {prog: run_shell(cmd) for prog, cmd in version_cmds.items()}

    malformed = [k for k, v in versions if v == '']

    if malformed:

        user_msg('The following dependencies are missing:')

        for mal in malformed:
            user_msg(mal)

        sys.exit(65)

    return versions

@logtime('SNAP database validation')
def verify_snap_databases(*db_dirs):
    '''Verifies SNAP databases.

    If databases are malformed, print which are broken and exit with code 65.
    '''

    def verify_database(db_dir):

        db_dir_path = os.path.abspath(db_dir)
        genomes_glob = os.path.join(db_dir_path, '*', 'Genome')

        return all(os.access(f, os.F_OK) for f in glob.glob(genomes_glob))

    malformed = [db_dir for db_dir in db_dirs if not verify_database(db_dir)]

    if malformed:

        user_msg('The following databases are malformed:')

        for mal in malformed:
            user_msg(mal)

        sys.exit(65)

    else:
        user_msg('SNAP Databases are OK.')

@logtime('RAPSearch database validation')
def verify_rapsearch_databases(viral: Path, nr: Path) -> None:
    '''Verifies RAPSearch databases.

    If databases are malformed, print which are broken at exit with code 65.
    '''

    def verify_rapsearch_database(db: Path) -> bool:
        '''Tests that rapsearch database exists'''
        return db.is_file() and db.with_suffix('.info').is_file()

    malformed = [db for db in (viral, nr) if not verify_rapsearch_database(db)]

    if malformed:
        user_msg('The following databases are malformed:')

        for mal in malformed:
            user_msg(mal)

        sys.exit(65)

    else:
        user_msg('RAPSearch databases are OK.')

@logtime('Fastq validation')
def validate_fastqs(fastq_file, logfile, mode):

    assert mode in range(4), 'Invalid fastq validation mode'
    modes = (
        '',
        'fastQValidator --file {} --printBaseComp \
                --avgQual --disableSeqIDCheck > {}',
        'fastQValidator --file {} --printBaseComp --avgQual > {}',
        'fastQValidator --file {} --printBaseComp --avgQual > {}'
        )

    # TODO
    if mode:

        try:
            cmd = modes[mode]
            result = run_shell(cmd.format(fastq_file, logfile))

        except subprocess.CalledProcessError:
            msg = '{} appears to be invalid. Check {} for details'
            user_msg(msg.format(fastq_file, logfile))

            if mode == 3:

                user_msg('Continuing anyway...')

            else:
                sys.exit(65)

@logtime('Sample and DB Validation')
def validation(inputfile: Path, workdir: Path, mode: str, fastq_log: Path,
               validation_mode: int, snap_db_dir: Path,
               rapsearch_vir_db: Path, rapsearch_nr_db: Path):

    sample = ensure_fastq(inputfile, workdir, mode)

    validate_fastqs(sample, fastq_log, validation_mode)

    verify_snap_databases(*snap_db_dir.glob('*'))

    verify_rapsearch_databases(rapsearch_vir_db, rapsearch_nr_db)

    return sample

@logtime('Total Runtime')
def surpi(sample: Path, workdir: Path, temp_dir: Path, fastq_type: str,
          quality_cutoff: int, length_cutoff: int, adapter_set: str,
          crop_start: int, crop_length: int, edit_distance: int,
          snap_db_dir: Path, tax_db_dir: Path, ribo_dir: Path, cache_reset,
          contig_cutoff:int, abyss_kmer: int, ignore_barcodes: bool,
          comprehensive: bool, rapsearch_mode: str, rapsearch_vir_db: Path,
          rapsearch_nr_db: Path, vir_cutoff: int, nr_cutoff: int,
          fast: bool, evalue: str, cores: int):
    '''Master function for SURPI pipeline'''

    preprocessed = preprocess(sample, workdir, temp_dir, adapter_set,
                              fastq_type, quality_cutoff, length_cutoff,
                              crop_start, crop_length)

    subtracted_fastq = host_subtract(preprocessed, snap_db_dir, edit_distance,
                                     temp_dir, cores)

    viruses, viruses_fastq, uniqunmatched = snap(sample, workdir, snap_db_dir,
                                                 tax_db_dir, ribo_dir, cores,
                                                 edit_distance, cache_reset,
                                                 comprehensive, temp_dir)

    assembled = assemble(viruses_fastq, uniqunmatched, workdir, temp_dir,
                         sample, contig_cutoff, abyss_kmer, ignore_barcodes,
                         cores)

    if rapsearch_mode == 'viral':
        rapsearch_viral(uniqunmatched, workdir, rapsearch_vir_db,
                        rapsearch_nr_db, vir_cutoff, nr_cutoff, fast,
                        assembled, tax_db_dir, viruses, evalue, cores)
    else:
        rapsearch_nr(uniqunmatched, workdir, assembled, rapsearch_nr_db,
                     tax_db_dir, nr_cutoff, fast, cores)

def main():
    pass
    #args = arguments()
    # TODO add taxonomy DB verification
    # TODO verification mode


if __name__ == '__main__':
    main()
