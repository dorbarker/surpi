from datetime import datetime
import sys
import functools
import fileinput

def user_msg(*args):
    '''Wrapper for print() that prints to stderr'''

    print(*args, file=sys.stderr)

def logtime(name):
    '''Function decorator that print to stderr the runtime of the
    decorated function.
    '''

    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            msg = 'Elapsed time for {}: {}'
            before = datetime.now()

            result = func(*args, *kwargs)

            after  = datetime.now()

            user_msg(msg.format(name, after - before))
            return result

        return wrapper
    return decorator

def run_shell(cmd, **kwargs):
    '''Runs a shell command using the subprocess module and returns the
    command's stdout.'''

    out = subprocess.check_output(cmd, shell=True,
                                  universal_newlines=True, **kwargs)
    return out.strip()

def concatenate(*files, output):
    '''Concatenate files together and write them to output'''

    with open(output, 'w') as o, fileinput.input(files)  as i:
        for line in i:
            o.write(line)
