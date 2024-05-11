#! /usr/bin/env python

from tqdm import trange
from subprocess import Popen, PIPE, STDOUT


def mrange(a, b=None, c=1, DEBUG=False):
    if DEBUG:
        if b is None:
            return trange(a, leave=False)
        else:
            return trange(a, b, c, leave=False)
    else:
        if b is None:
            return range(a)
        else:
            return range(a, b, c)


def mexec(command, DEBUG=False):
    if DEBUG:
        print(f'Executing {command}')
    if isinstance(command, str):
        command = command.split()
    captured_output = []
    with Popen(command, stdout=PIPE,
               stderr=STDOUT, text=True) as process:
        while True:
            output_line = process.stdout.readline()
            if not output_line and process.poll() is not None:
                break
            if DEBUG:
                print(output_line, end='')
            captured_output.append(output_line)
    return ''.join(captured_output).strip()
