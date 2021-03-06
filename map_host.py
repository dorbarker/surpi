'''Matches sample reads to known host sequence using SNAP'''

import tempfile
import subprocess
from pathlib import Path
from utilities import user_msg, logtime, write_fastq_generators, annotated_to_fastq

def vmtouch(snap_db: Path) -> bool:
    '''Calls vmtouch to check if the SNAP database is loaded into memory'''

    cmd = ('vmtouch', '-f', '-m500G', str(snap_db))
    output = subprocess.check_output(cmd, universal_newlines=True).splitlines()

    for line in output:
        if 'Resident Pages' in line and line.split()[4] == '100%':
            in_memory = True
            break
    else:
        in_memory = False

    return in_memory

def set_snap_cache_option(snap_db: Path) -> str:
    '''Sets SNAP options based on whether the database is loaded into memory'''

    if vmtouch(snap_db):

        user_msg('{} is cached'.format(snap_db))

        snap_cache_option = ''

    else:

        user_msg('{} is not cached'.format(snap_db))

        snap_cache_option = '-pre'

    return snap_cache_option

def host_subtract(preprocessed: Path, snap_db_dir: Path, edit_distance: int,
                  temp_dir: Path, cores: int) -> Path:
    '''Subtracts preprocessed reads from host SNAP databases.

    Writes the subtracted fastq and returns its file path.
    '''
    @logtime('Host subtraction')
    def subtract(snap_db: Path, to_subtract: Path,
                 output: Path, pre: str) -> None:
        '''Prepare and execute the SNAP search'''

        subtract_cmd = ('snap-aligner', 'single', snap_db, to_subtract,
                        '-o', '-sam', output,
                        '-t', cores, '-x', '-f', '-h', '250',
                        '-d', edit_distance, '-n', '25', '-F', 'u',
                        '-map', pre)

        subprocess.call([str(arg) for arg in subtract_cmd if arg])
        # end subtract()

    subtracted_fastq = preprocessed.with_suffix('')\
                                   .with_suffix('.unmatched.fastq')

    to_subtract = preprocessed

    with tempfile.TemporaryDirectory(dir=str(temp_dir)) as ephemeral:

        for i, snap_db in enumerate(snap_db_dir.glob('*')):

            temp_output = Path(ephemeral) / '{}.sam'.format(str(i))

            snap_cache_option = set_snap_cache_option(snap_db)

            subtract(snap_db, to_subtract, temp_output, snap_cache_option)

            to_subtract = temp_output

        fastq_lines = annotated_to_fastq(to_subtract, False)

        write_fastq_generators(subtracted_fastq, fastq_lines)

    return subtracted_fastq
