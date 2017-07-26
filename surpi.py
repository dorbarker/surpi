#!/usr/bin/env python3

# system imports
import json
import subprocess
import sys
import textwrap
from argparse import ArgumentParser, RawTextHelpFormatter, Namespace
from multiprocessing import cpu_count
from pathlib import Path
from typing import Dict

# surpi imports
from preprocessing import preprocess
from utilities import run_shell, user_msg, logtime, pigz
from snap_to_nt import snap
from map_host import host_subtract
from assembly import assemble
from rapsearch import rapsearch_viral, rapsearch_nr

def arguments():
    '''Reads commandline arguments and parses them.

    If --config is set, arguments will be read from the given config file
    '''

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)

    parser.add_argument('--fastq-type',
                        choices=('illumina', 'sanger'),
                        default='sanger',
                        help='Quality score mode [sanger]')

    parser.add_argument('--adapter-set',
                        choices=('truseq', 'nextera',
                                 'nexsolb', 'nexsoltruseq'),
                        default='truseq',
                        help='Adapter set to use [truseq]')


    parser.add_argument('--verify-fastq',
                        choices=(0, 1, 2, 3),
                        default=1,
                        help=textwrap.dedent('''\
                        0 = skip fastq validation
                        1 = validate fastqs; quit on fail
                        2 = validate fastqs; ensure unique names; quit on fail
                        3 = validate fastqs; ensure unique name; warn on fail
                        [1]
                        '''))

    parser.add_argument('--fast',
                        action='store_false',
                        dest='comprehensive',
                        help=textwrap.dedent('''
                        Run in fast mode instead of comprehensive mode.
                        Comprehensive mode runs:
                            1. SNAP to NT
                            2. de novo contig assembly
                            3. RAPSearch to viral proteins or NCBI NR

                        Fast mode will SNAP align to curated bacterial
                        and viral databases
                        [comprehensive]
                        '''))

    parser.add_argument('--skip-preprocess',
                        action='store_true',
                        help='Skip preprocessing step [off]')

    parser.add_argument('--length-cutoff',
                        type=int,
                        default=50,
                        help=textwrap.dedent('''\
                                Discard any post-trimming sequence
                                less than this length [50]'''))

    parser.add_argument('--crop-start',
                        type=int,
                        default=10,
                        help='Start position for read cropping [10]')

    parser.add_argument('--crop-length',
                        type=int,
                        default=75,
                        help='Length to crop after start position [75]')

    parser.add_argument('--quality-cutoff',
                        type=int,
                        default=18,
                        help='Quality cutoff for cutadapt [18]')

    parser.add_argument('--edit-distance',
                        type=int,
                        default=12,
                        help='SNAP edit distance for host genome [12]')

    parser.add_argument('--rapsearch-method',
                        choices=('viral', 'nr'),
                        default='viral',
                        help='RAPSearch database method to use [viral]')

    parser.add_argument('--rapsearch-cutoff',
                        type=float,
                        default='1',
                        help='E-value cutoff for RAPSearch [1]')

    parser.add_argument('--fast-rapsearch',
                        action='store_true',
                        help=textwrap.dedent('''\
                                Run RAPSearch much faster (10 to 30-fold) at
                                the cost of sensitivity [off]'''))

    parser.add_argument('--vir-cutoff',
                        type=float,
                        default=1.0,
                        help='Virus cutoff for RAPSearch [1]')

    parser.add_argument('--nr-cutoff',
                        type=float,
                        default=1.0,
                        help='NR cutoff for RAPSearch [1]')

    parser.add_argument('--abyss-kmer',
                        type=int,
                        default=34,
                        help='kmer size for ABySS de novo assembly [34]')

    parser.add_argument('--ignore-barcodes',
                        action='store_true',
                        help='Assemble all barcodes into a single assembly [off]')

    parser.add_argument('--evalue',
                        default='1e-15',
                        help=textwrap.dedent('''\
                              E-value cutoff for BLASTn.
                              Will parse notation such as "1e-15" [1e-15]'''))

    parser.add_argument('--cores',
                        type=int,
                        default=cpu_count(),
                        help='CPU cores to use [{}]'.format(cpu_count()))

    parser.add_argument('--temp',
                        type=Path,
                        default='/tmp/',
                        help='Temporary working directory [/tmp/]')

    parser.add_argument('--config',
                        type=Path,
                        help='Config file [optional]')

    parser.add_argument('reference',
                        type=Path,
                        help='Reference data directory')

    parser.add_argument('workdir',
                        type=Path,
                        help='Location for working files and results')

    parser.add_argument('sample',
                        type=Path,
                        help='Sample to analyze')

    args = parser.parse_args()

    ignore_args = set(('reference', 'workdir', 'sample', 'config'))
    if args.config:

        new_args = Namespace()
        with args.config.open('r') as config_file:

            data = json.load(config_file)

            for key, value in data.items():

                if key not in ignore_args:

                    if value == 'True':
                        new_value = True

                    elif value == 'False':
                        new_value = False

                    else:
                        new_value = type(getattr(args, key))(value)

                    setattr(new_args, key, new_value)

        new_args.reference = args.reference
        new_args.workdir = args.workdir
        new_args.sample = args.sample

        args = new_args

    new_config = (args.workdir / args.sample.name).with_suffix('.config')
    with new_config.open('w') as config_file:

        to_dump = {k: str(v) for k, v in args._get_kwargs()
                   if k not in ignore_args}

        json.dump(to_dump, config_file, indent=4)

    return args

def ensure_fastq(inputfile: Path, workdir: Path) -> Path:
    '''If mode is "fasta", convert inputfile to fastq.

    In either case, guarantee that file extension is ".fastq"
    '''

    fq_name = workdir / inputfile.with_suffix('.fastq').name

    # FASTA
    if inputfile.suffix in ('.fasta', '.fas'):


        cmd = ('fasta_to_fastq', str(inputfile))

        fq_name.write_text(subprocess.check_output(cmd,
                                                   universal_newlines=True))
    # Compressed FASTQ
    elif inputfile.suffix == '.gz':

        fq_name = workdir / inputfile.with_suffix('').name

        pigz(inputfile, fq_name)

    # Uncompressed FASTQ
    else:

        try:
            fq_name.symlink_to(inputfile)

        except FileExistsError:
            fq_name.unlink()
            fq_name.symlink_to(inputfile)

    return fq_name

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

#        db_dir_path = os.path.abspath(db_dir)
#        genomes_glob = os.path.join(db_dir_path, '*', 'Genome')
        db_dir_path = db_dir.resolve()
        return db_dir_path.glob('*/Genome')
        #return all(os.access(f, os.F_OK) for f in glob.glob(genomes_glob))

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
    '''Using external program fastQValidator, ensure input fastq files are
    valid.
    '''

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
            run_shell(cmd.format(fastq_file, logfile))

        except subprocess.CalledProcessError:
            msg = '{} appears to be invalid. Check {} for details'
            user_msg(msg.format(fastq_file, logfile))

            if mode == 3:

                user_msg('Continuing anyway...')

            else:
                sys.exit(65)

@logtime('Sample and DB Validation')
def validation(inputfile: Path, workdir: Path, validation_mode: int,
               reference: Path, comprehensive: bool) -> Path:

    snap_db = 'comp_snap' if comprehensive else 'fast_snap'
    snap_db_dir = reference / snap_db
    rapsearch_nr_db = reference / 'rapsearch' / 'rapsearch_nr'
    rapsearch_vir_db = reference / 'rapsearch' / 'rapsearch_viral'
    fastq_log = workdir / 'fastq_validation.log'

    sample = ensure_fastq(inputfile, workdir)

    validate_fastqs(sample, fastq_log, validation_mode)

    verify_snap_databases(*snap_db_dir.glob('*'))

    verify_rapsearch_databases(rapsearch_vir_db, rapsearch_nr_db)

    return sample

@logtime('Total Runtime')
def surpi(sample: Path, workdir: Path, temp_dir: Path, edit_distance: int,
          reference: Path, abyss_kmer: int, ignore_barcodes: bool,
          comprehensive: bool, rapsearch_mode: str, vir_cutoff: int,
          nr_cutoff: int, rap_fast: bool, evalue: str, cores: int) -> None:
    '''Master function for SURPI pipeline'''

    snap_db = 'comp_snap' if comprehensive else 'fast_snap'

    host_snap_dir = reference / 'host_snap'
    rapsearch_nr_db = reference / 'rapsearch' / 'rapsearch_nr'
    rapsearch_vir_db = reference / 'rapsearch' / 'rapsearch_viral'
    ribo_dir = reference / 'riboclean_snap'
    tax_db_dir = reference / 'taxonomy'
    snap_db_dir = reference / snap_db

    subtracted_fastq = host_subtract(preprocessed, host_snap_dir,
                                     edit_distance, temp_dir, cores)

    viruses, viruses_fastq, uniqunmatched = snap(subtracted_fastq, workdir,
                                                 snap_db_dir, tax_db_dir,
                                                 ribo_dir, cores,
                                                 edit_distance, cache_reset,
                                                 comprehensive, temp_dir)

    if comprehensive:

        assembled = assemble(viruses_fastq, uniqunmatched, workdir, temp_dir,
                             sample, abyss_kmer, ignore_barcodes, cores)

        if rapsearch_mode == 'viral':
            rapsearch_viral(uniqunmatched, workdir, rapsearch_vir_db,
                            rapsearch_nr_db, vir_cutoff, nr_cutoff,
                            assembled, tax_db_dir, viruses, evalue, cores)
        else:
            rapsearch_nr(uniqunmatched, workdir, assembled, rapsearch_nr_db,
                         tax_db_dir, nr_cutoff, rap_fast, cores)

def main():

    args = arguments()

    sample = validation(args.sample, args.workdir, args.verify_fastq,
               args.reference, args.comprehensive)

    if not args.skip_preprocess:
        sample = preprocess(sample, args.workdir, args.temp, args.adapter_set,
                            args.fastq_type, args.quality_cutoff,
                            args.length_cutoff, args.crop_start,
                            args.crop_length )

    surpi(sample, args.workdir, args.temp, args.edit_distance, args.reference,
          args.abyss_kmer, args.ignore_barcodes, args.comprehensive,
          args.rapsearch_method, args.vir_cutoff, args.nr_cutoff,
          args.fast_rapsearch, args.evalue, args.cores)

if __name__ == '__main__':
    main()
