#!/usr/bin/env python3

import argparse

from . import dlog


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="discretelog", description="Calculate discrete logarithm."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False, help="Verbose mode"
    )
    parser.add_argument("-g", "--base", type=int, help="Logarithm base")
    parser.add_argument("-n", "--modulo", type=int, help="Modulo")
    parser.add_argument("element", type=int, help="Element")

    args = parser.parse_args()

    d, _dt = dlog(args.base, args.element, args.modulo, args.verbose)
    print(d)


if __name__ == "__main__":
    main()
