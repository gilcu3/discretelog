#! /usr/bin/env python3 -u

import sys

from . import dlog

if __name__ == "__main__":
    assert len(sys.argv) == 4, "Expecting parameters g (base), v (value), n (modulo)"
    g, v, n = tuple(map(int, sys.argv[1:]))
    d, dt = dlog(g, v, n)
    print(d)
