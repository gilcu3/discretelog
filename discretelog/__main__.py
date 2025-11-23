#!/usr/bin/env python3

import argparse

from . import dlog


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="discretelog", description="Calculate discrete logarithm."
    )
    parser.add_argument("-g", type=int, help="Logarithm base")
    parser.add_argument("-n", type=int, help="Modulo")
    parser.add_argument("v", type=int, help="Element")

    args = parser.parse_args()

    d, _dt = dlog(args.g, args.v, args.n)
    print(d)


if __name__ == "__main__":
    main()
