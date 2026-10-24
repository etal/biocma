"""Command-line programs installed with biocma.

Each ``*_main`` function is the entry point of one console script. It takes
an optional argument list, as ``sys.argv[1:]`` would hold, and returns the
process exit status. A file argument of ``-`` means standard input or output.
"""

import argparse
import sys

from . import cma


def _input(path):
    """Return standard input for '-', otherwise the path itself."""
    return sys.stdin if path == '-' else path


def _fail(prog, exc):
    """Report an error on standard error and return a failing exit status."""
    print(f"{prog}: error: {exc}", file=sys.stderr)
    return 1


def blocks_main(argv=None):
    """cma-blocks: write the aligned columns of a CMA file as FASTA."""
    parser = argparse.ArgumentParser(
        prog='cma-blocks',
        description="Remove the unaligned insert columns from a CMA file and "
        "write the remaining aligned columns to standard output as a FASTA "
        "alignment. The sequences of all blocks in the file are written in "
        "order as one FASTA stream.")
    parser.add_argument(
        'infile', nargs='?', default='-',
        help="CMA file to read (default: standard input).")
    args = parser.parse_args(argv)

    try:
        for block in cma.parse(_input(args.infile)):
            for seq in block['sequences']:
                header = seq['id']
                if seq['description']:
                    header += ' ' + seq['description']
                aligned = ''.join(c for c in seq['seq'] if not c.islower())
                sys.stdout.write(f">{header}\n{aligned}\n")
    except OSError as exc:
        return _fail(parser.prog, exc)
    return 0
