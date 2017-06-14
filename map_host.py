#!/usr/bin/env python3

import tempfile
import subprocess
from pathlib import Path
from utilities import user_msg, logtime
from shutil import copy
def vmtouch(db):

    cmd = ('vmtouch', '-m500G', db)
    output = subprocess.check_output(cmd, universal_newlines=True).splitlines()

    for line in output:
        if 'Resident Pages' in line:
            l = line.split()

            return l[4] == '100%'

def host_subtract(sample, output, snap_db_dir, cores, tempdir, edit_distance):

    @logtime('Host subtraction')
    def subtract(snap_db, to_subtract, output, pre):

        cmd = ('snap-dev', snap_db, to_subtract, '-o', 'sam', output, '-t',
                str(cores), '-x', '-f', '-h', '250', '-d', edit_distance,
                '-n', '25', '-F', 'u', '-map', pre)

        subprocess.call(cmd)

    to_subtract = sample.with_suffix('.preprocessed.fastq')

    with tempfile.TemporaryDirectory(dir=temp) as d:

        for i, snap_db in enumerate(snap_db_dir.glob('*')):

            if vmtouch(snap_db):

                user_msg('{} is cached'.format(snap_db))

                snap_cache_option = ''

            else:

                user_msg('{} is not cached'.format(snap_db))

                snap_cache_option = '-pre'

            temp_output = Path(d) / i

            subtract(snap_db, to_subtract, temp_output, snap_cache_option)

            to_subtract = temp_output
