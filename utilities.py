from datetime import datetime
import sys
import functools

def user_msg(*args):
    print(*args, file=sys.stderr)

def logtime(name):
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
