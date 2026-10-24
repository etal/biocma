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


def convert_main(argv=None):
    """cma-convert: convert between CMA and other alignment formats."""
    # Imported here so that the other commands do not load Biopython.
    from . import biocma

    parser = argparse.ArgumentParser(
        prog='cma-convert',
        description="Convert a CMA file to another alignment format, or an "
        "alignment in another format to CMA, and write the result to standard "
        "output. FORMAT is a Bio.AlignIO format name such as fasta, clustal "
        "or stockholm.")
    direction = parser.add_mutually_exclusive_group(required=True)
    direction.add_argument(
        '-o', '--outformat', metavar='FORMAT',
        help="Convert CMA input to FORMAT.")
    direction.add_argument(
        '-i', '--informat', metavar='FORMAT',
        help="Convert FORMAT input to CMA.")
    parser.add_argument(
        'infile', nargs='?', default='-',
        help="Alignment file to read (default: standard input).")
    args = parser.parse_args(argv)

    if args.outformat is not None:
        in_fmt, out_fmt = 'cma', args.outformat
    else:
        in_fmt, out_fmt = args.informat, 'cma'
    try:
        biocma.convert(_input(args.infile), in_fmt, sys.stdout, out_fmt)
    except (OSError, ValueError) as exc:
        return _fail(parser.prog, exc)
    return 0
