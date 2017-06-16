'''Matches sample reads to known host sequence using SNAP'''

import tempfile
import subprocess
from pathlib import Path
from utilities import user_msg, logtime, annotated_to_fastq

def vmtouch(snap_db):
    '''Calls vmtouch to check if the SNAP database is loaded into memory'''

    cmd = ('vmtouch', '-m500G', snap_db)
    output = subprocess.check_output(cmd, universal_newlines=True).splitlines()

    for line in output:
        if 'Resident Pages' in line:

            return line.split()[4] == '100%'

def set_snap_cache_option(snap_db):
    '''Sets SNAP options based on whether the database is loaded into memory'''

    if vmtouch(snap_db):

        user_msg('{} is cached'.format(snap_db))

        snap_cache_option = ''

    else:

        user_msg('{} is not cached'.format(snap_db))

        snap_cache_option = '-pre'

    return snap_cache_option

def host_subtract(preprocessed, snap_db_dir, edit_distance, temp_dir, cores):
    '''Subtracts preprocessed reads from host SNAP databases.

    Writes the subtracted fastq and returns its file path.
    '''
    @logtime('Host subtraction')
    def subtract(snap_db, to_subtract, output, pre):
        '''Prepare and execute the SNAP search'''

        subtract_cmd = ('snap-aligner', 'single', snap_db, to_subtract,
                        '-o', '-sam', output,
                        '-t', cores, '-x', '-f', '-h', '250',
                        '-d', edit_distance, '-n', '25', '-F', 'u',
                        '-map', pre)

        subprocess.call(map(str, subtract_cmd))
        # end subtract()

    subtracted_fastq = preprocessed.with_suffix('.fastq')

    to_subtract = preprocessed

    with tempfile.TemporaryDirectory(dir=temp_dir) as temp_dir_:

        for i, snap_db in enumerate(snap_db_dir.glob('*')):

            snap_cache_option = set_snap_cache_option(snap_db)

            temp_output = Path(temp_dir_) / i

            subtract(snap_db, to_subtract, temp_output, snap_cache_option)

            to_subtract = temp_output

        subtracted_fastq.write_text(annotated_to_fastq(to_subtract))

    return subtracted_fastq
