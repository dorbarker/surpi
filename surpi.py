#!/usr/bin/env python3

import glob
import os
import argparse
import subprocess
from multiprocessing import cpu_count
import sys
from pathlib import Path

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

    # TODO: Add AWS-related values

    return parser.parse_args()

def user_msg(*msg):

    print(*msg, file=sys.stderr)

def ensure_fastq(inputfile, mode):
    '''If mode is "fasta", convert inputfile to fastq.

    In either case, guarantee that file extension is ".fastq"
    '''

    assert mode in ('fastq', 'fasta'), 'Mode must be "fasta" or "fastq"'

    name, ext = os.path.splitext(inputfile)

    fq_name = name + '.fastq'

    if mode == 'fasta':

        cmd = ('fasta_to_fastq', inputfile)

        fastq = subprocess.check_output(cmd, universal_newlines=True)

        with open(fq_name, 'w') as f:
            f.write(fastq)

    else:

        os.symlink(inputfile, fq_name)

def free_memory(cache_reset):

    freemem = float(run_shell("free -g | awk 'NR==2{print $4}'"))

    if freemem < cache_reset:
        run_shell('dropcache')

def get_memory_limit():
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

def run_shell(cmd, **kwargs):
    out = subprocess.check_output(cmd, shell=True, universal_newlines=True, **kwargs)
    return out.strip()

def check_program_versions():
    '''Return the versions of each external dependency as a dict

    {program: version}
    '''

    version_cmds = {
        'gt':        "gt -version | head -1 | awk '{print $3}'",
        'seqtk':     "seqtk 2>&1 | head -3 | tail -1 | awk '{print $2}'",
        'cutadapt':  "cutadapt --version",
        'prinseqlite':  "prinseq-lite.pl --version 2>&1 | awk '{print $2}'",
        'snap':      "snap 2>&1 | grep version | awk '{print $5}'",
        'snap-dev':  "snap-dev 2>&1 | grep version | awk '{print $5}'",
        'rapsearch': "rapsearch 2>&1 | head -2 | tail -1 | awk '{print $2}'",
        'abyss_pe':  "abyss-pe version | head -2 | tail -1 | awk '{print $3}'",
        'abyss_p':   "ABYSS-P  --version | head -1 | awk '{print $3}'",
        'minimo':    "Minimo -h | tail -2 | awk '{print $2}'"
        }

    versions = {prog: run_shell(cmd) for prog, cmd in version_cmds.items()}

    malformed = [k for k, v in versions if v == '']

    if malformed:

        print('The following dependencies are missing:', file=sys.stderr)

        for mal in malfored:
            print(mal, file=sys.stderr)

        sys.exit(65)

    return versions

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

        print('The following databases are malformed:', file=sys.stderr)

        for mal in malformed:
            print(mal, file=sys.stderr)

        sys.exit(65)

    else:
        print('SNAP Databases are OK.', file=sys.stderr)

def verify_rapsearch_databases(viral, nr):
    '''Verifies RAPSearch databases.

    If databases are malformed, print which are broken at exit with code 65.
    '''

    def verify_rapsearch_database(db):
        return os.access(db, os.F_OK) and os.access(db + '.info', os.F_OK)

    malformed = [db for db in (viral, nr) if not verify_rapsearch_database(db)]

    if malformed:
        print('The following databases are malfored:', file=sys.stderr)

        for mal in malformed:
            print(mal, file=sys.stderr)

        sys.exit(65)

    else:
        print('RAPSearch databases are OK.', file=sys.stderr)

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
            result = run_shell(modes[mode].format(fastq_file, logfile))

        except CalledProcessError:
            msg = '{} appears to be invalid. Check {} for details'.format(fastq_file, logfile)
            print(msg, file=sys.stderr)

            if mode == 3:

                print('Continuing anyway...', file=sys.stderr)

            else:
                sys.exit(65)

def preprocess(name, quality, length_cutoff, cores, adapter_set, start,
               crop_length, temp_dir):
    # TODO user messages

    cmd = 'preprocess_ncores.sh {name}.fastq {qual} N {len_cut} {cores} Y N \
            {adapters} {start} {crop_len} {temp_dir} >& {name}.preprocess.log'

    run_shell(cmd.format(name=name, qual=quality, len_cut=length_cutoff,
                         cores=cores, adapters=adapter_set, start=start,
                         crop_len=crop_length, temp_dir=temp_dir))
    try:
        os.path.getsize(name + 'preprocess.log')
        os.path.getsize(name + '.cutadapt.fastq')
    except OSError:
        user_msg('Preproessing {} appears to have failed'.format(name))
        sys.exit(65)

# TODO may need extensive modification to generalize from human to other hosts
def human_mapping(fastq, d_human, snap_subtraction_dir, cores):

    human_basefile = Path(fastq).with_suffix('.preprocessed.s20.h250n25d{}xfu'.format(d_human))

    to_subtract = Path(fastq).with_suffix('.preprocessed.fastq')

    subtracted_output = Path(fastq).with_suffix('human.snap.unmatched.sam')

    subtraction_counter = 0

    vmtouch= "vmtouch -m500G -f {} | grep 'Resident Pages' | awk '{print $5}'"

    snapdev = 'snap-dev single {sub_db} {to_sub} -o -sam {out} -t {cores} \
              -x -f -h 250 -d {d_human} -n 25 -F u {snap_cache_opt}'

    for subtract_db in Path(snap_subtraction_dir).glob('*'):

        outname = str(subtracted_output) + '.' + subtraction_counter + '.sam'

        subtraction_counter += 1

        snap_db_cached = run_shell(vmtouch.format(db))

        snap_cache_opt = ' -map ' if snap_db_cached == '100%' else ' -pre -map '

        # wtf? fix later
        snapdev_cmd = snapdev.format(sub_db=subtract_db, to_sub=to_subtract,
                                     out=outname,
                                     cores=cores, d_human=d_human, snap_cache_opt=snap_cache_opt)

        to_subtract = outname

    # dear god; must replace
    write_unmatched = "grep -Ev \"^@\" {} | awk '{if($3 == \"*\") print \"@\"$1\"\n\"$10\"\n\"\"+\"$1\"\n\"$11}' > \
            $echo {} | sed 's/\(.*\)\..*/\1/').fastq".format(str(human_basefile) + '.human.snap.unmatched.sam')
    run_shell(write_unmatched)

def snap_to_nt(fastq, run_mode, host_basefile, snap_edit_distance, snap_comp_db,
               snap_fast_db, cores):

    # SNAP unmatched sequences to NT
    if not Path(fastq).with_suffix('.NT.snap.sam').exists():

        snap_nt = 'snap_nt.sh {host_unmatched} {snap_db} {cores} \
                       {cache_reset} {snap_distance}'

        if run_mode == 'comprehensive':


            run_shell(snap_nt.format(host_unmatched=str(host_basefile) + '.human.snap.unmatched.fastq',
                                     snap_db=snap_comp_db, cores=cores, cache_reset=cache_reset,
                                     snap_distance = snap_edit_distance))

        else:  # run_mode == 'fast'

            run_shell(snap_nt.format(host_unmatched=str(host_basefile) + '.human.snap.unmatched.fastq',
                                     snap_db=snap_fast_db, cores=cores, cache_reset=cache_reset,
                                     snap_distance = snap_edit_distance))




    matches = run_shell('grep -Ev "^@" {}'.format(Path(fastq).with_suffix('')))


    matches_sam = Path(fastq).with_suffix('.NT.snap.matched.sam')
    unmatches_sam = Path(fastq).with_suffix('.NT.snap.unmatched.sam')

    with open(matches_sam, 'w') as m:

        matched = run_shell("awk '{if($3 != \"*\") print }'", input=matches)
        m.write(matched)

    with open(unmatches_sam, 'w') as u:

        unmatched = run_shell("awk '{if($3 == \"*\") print }'", input=matches)
        u.write(unmatched)

    if Path(fastq).with_suffix('.NT.snap.matched.all.annotated'):

        # convert to FASTQ and retrieve full-length sequences

        extract_fq_header = 'extractHeaderFromFastq_ncores.sh \
                {cores} {cutadapt_fastq} {nt_snap_match} \
                {nt_snap_match_fulllen} {nt_snap_unmatch} \
                {nt_snap_unmatch_fulllen}'.format(
                        cores=cores,
                        cutadapt_fastq=Path(fastq).with_suffix('.cutadapt.fastq'),
                        nt_snap_match=matches_sam,
                        nt_snap_match_fulllen=matches_sam.with_suffix('.fulllength.fastq'),
                        nt_snap_unmatch=unmatches_sam,
                        nt_snap_unmatch_fulllen=unmatches_sam.with_suffix('.fulllength.fastq')
                        )

        run_shell(extract_fq_header)

        # line 900 in original

def main():

    args = arguments()
    # TODO add taxonomy DB verification
    # TODO handle slave instances (line 633-637, 687-703 in original)
    # TODO verification mode


if __name__ == '__main__':
    main()
