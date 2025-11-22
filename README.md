# discretelog

A pure python package to compute discrete logs in the ring of integers
modulo n. It aims to provide relatively performant code with simple
implementations of well-known algorithms, inspired by the popular `primefac`
library.

## Usage

To find the discrete log of `v=641629670911834423534` modulo
`n=1540571422742786915303` with base `g=25`:

```bash
python -m discretelog 25 641629670911834423534 1540571422742786915303
```

### Installation

The library is available in `pypi`, therefore can be installed with:

```bash
python -m pip install discretelog
```

### Building from source

Using `poetry` to generate a wheel file:

```bash
git clone https://github.com/gilcu3/discretelog
cd discretelog
poetry build
```

## Status

This is a work in progress. Several issues will be worked on in the future:

- Looking at a similar library within `Pari/GP`, probably better performance
can be achieved
- Several parts can be trivially made parallel
- For the moment documentation is missing
- The library should either fail gracefully or continue computing indefinitely
in case of big inputs. At the moment it simply crashes

### Discrete log algorithms implemented

- Baby steps giant steps
- Pollard rho
- Pohlig Hellman
- Naive index calculus
- Linear sieve index calculus

## External tools

For large inputs the library uses external `c++` implementations of state
of the art algorithms. In order to use them the user needs to install
them separately. See `docs/external.md` for details
