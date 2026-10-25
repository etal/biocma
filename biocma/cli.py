"""Command-line programs installed with biocma.

Each ``*_main`` function is the entry point of one console script. It takes
an optional argument list, as ``sys.argv[1:]`` would hold, and returns the
process exit status. A file argument of ``-`` means standard input or output.
"""

import argparse
import contextlib
import functools
import os
import sys

from . import cma, utils


def _input(path):
    """Return standard input for '-', otherwise the path itself."""
    return sys.stdin if path == "-" else path


def _output(path):
    """Return standard output for '-', otherwise the path itself."""
    return sys.stdout if path == "-" else path


def _console(main):
    """Run a command and report its input and output errors in one line.

    The command function NAME_main takes the argument list and the program
    name, cma-NAME. OSError and ValueError, which include unreadable files,
    text that is not valid in the locale's encoding and malformed CMA input,
    are reported on standard error as 'cma-NAME: error: ...' with exit status
    1. When the reader of standard output exits early, as head does, the
    command exits with status 1 and prints nothing, following the note on
    SIGPIPE in the documentation of Python's signal module.
    """
    prog = "cma-" + main.__name__.removesuffix("_main")

    @functools.wraps(main)
    def wrapper(argv=None):
        try:
            status = main(argv, prog)
            # Inside the try, so that a closed pipe is not first found by the
            # flush at interpreter exit
            sys.stdout.flush()
        except BrokenPipeError:
            _discard_stdout()
            return 1
        except (OSError, ValueError) as exc:
            print(f"{prog}: error: {exc}", file=sys.stderr)
            return 1
        return status

    # The wrapper takes only the argument list; without this, introspection
    # would follow __wrapped__ and report the command function's signature
    del wrapper.__wrapped__
    return wrapper


def _discard_stdout():
    """Send standard output to the null device, as the SIGPIPE note advises.

    Output still buffered then cannot raise BrokenPipeError again when the
    interpreter flushes standard output at exit.
    """
    try:
        stdout = sys.stdout.fileno()
    except OSError:
        # A replacement without a file descriptor, as under pytest
        return
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, stdout)
    os.close(devnull)


def _add_infile(parser, what="CMA file"):
    """Add the optional input file argument that every command takes."""
    parser.add_argument(
        "infile",
        nargs="?",
        default="-",
        help=f"{what} to read, or '-' for standard input (the default).",
    )


@_console
def blocks_main(argv, prog):
    """cma-blocks: write the aligned columns of a CMA file as FASTA."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Remove the unaligned insert columns from a CMA file and "
        "write the remaining aligned columns to standard output as a FASTA "
        "alignment. The sequences of all blocks in the file are written in "
        "order as one FASTA stream.",
    )
    _add_infile(parser)
    args = parser.parse_args(argv)

    with contextlib.closing(cma.parse(_input(args.infile))) as blocks:
        for block in blocks:
            for seq in block["sequences"]:
                header = seq["id"]
                if seq["description"]:
                    header += " " + seq["description"]
                aligned = "".join(c for c in seq["seq"] if not c.islower())
                sys.stdout.write(f">{header}\n{aligned}\n")
    return 0


@_console
def convert_main(argv, prog):
    """cma-convert: convert between CMA and other alignment formats."""
    # Imported here so that the other commands do not load Biopython.
    from . import biocma

    parser = argparse.ArgumentParser(
        prog=prog,
        description="Convert a CMA file to another alignment format, or an "
        "alignment in another format to CMA, and write the result to standard "
        "output. FORMAT is a Bio.AlignIO format name such as fasta, clustal "
        "or stockholm.",
    )
    direction = parser.add_mutually_exclusive_group(required=True)
    direction.add_argument(
        "-o", "--outformat", metavar="FORMAT", help="Convert CMA input to FORMAT."
    )
    direction.add_argument(
        "-i", "--informat", metavar="FORMAT", help="Convert FORMAT input to CMA."
    )
    _add_infile(parser, "Alignment file")
    args = parser.parse_args(argv)

    if args.outformat is not None:
        in_fmt, out_fmt = "cma", args.outformat
    else:
        in_fmt, out_fmt = args.informat, "cma"
    biocma.convert(_input(args.infile), in_fmt, sys.stdout, out_fmt)
    return 0


@_console
def dedupe_main(argv, prog):
    """cma-dedupe: keep one record for each sequence ID in a CMA file."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Remove records with duplicate sequence IDs from each "
        "block of a CMA file and write the result to standard output. For "
        "each ID the record with the most aligned (uppercase) residues is "
        "kept, and of equally good records the last one. The number of "
        "records removed is reported on standard error.",
    )
    _add_infile(parser)
    args = parser.parse_args(argv)

    blocks = []
    dropped = 0
    with contextlib.closing(cma.parse(_input(args.infile))) as parsed:
        for block in parsed:
            kept = _best_by_id(block["sequences"])
            dropped += len(block["sequences"]) - len(kept)
            blocks.append(dict(block, sequences=kept))
    print(f"Dropped {dropped} duplicate sequences", file=sys.stderr)
    cma.write(blocks, sys.stdout)
    return 0


def _best_by_id(sequences):
    """Keep the best record for each sequence ID, in their original order.

    The best record has the most uppercase (aligned) residues; of equally good
    records, the last is kept.
    """
    best: dict[str, tuple[int, int]] = {}
    for idx, seq in enumerate(sequences):
        score = sum(c.isupper() for c in seq["seq"])
        if seq["id"] not in best or score >= best[seq["id"]][0]:
            best[seq["id"]] = (score, idx)
    keep = {idx for _score, idx in best.values()}
    return [seq for idx, seq in enumerate(sequences) if idx in keep]


@_console
def tweak_main(argv, prog):
    """cma-tweak: change block-level features of a CMA file."""
    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Change CMA-specific features of each block in a CMA file.",
        epilog=f"example:\n  {prog} foo.cma --iron --level 1 --name Testfam -o bar.cma",
    )
    _add_infile(parser)
    parser.add_argument(
        "-c",
        "--add-consensus",
        action="store_true",
        help="Insert a consensus sequence as the first record.",
    )
    parser.add_argument(
        "--iron",
        action="store_true",
        help="'Iron out' indels to remove deletion/insertion transitions.",
    )
    parser.add_argument("-l", "--level", type=int, help="Set the block level.")
    parser.add_argument("--name", type=_block_name, help="Set the block name.")
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="CMA file to write, or '-' for standard output (the default).",
    )
    args = parser.parse_args(argv)

    with contextlib.closing(cma.parse(_input(args.infile))) as parsed:
        blocks = [_tweak_block(block, args) for block in parsed]
    cma.write(blocks, _output(args.output))
    return 0


def _block_name(value):
    """Accept a block name that MAPGAPS and cma.parse can read back.

    MAPGAPS reads the name as the characters up to the first '(', and a line
    break would split the block header line.
    """
    if not value or any(c in value for c in "(\r\n"):
        raise argparse.ArgumentTypeError(
            "must be nonempty and contain no '(' or line break"
        )
    return value


def _tweak_block(block, args):
    """Apply the cma-tweak options in args to block in place and return it."""
    if args.add_consensus:
        block["sequences"].insert(0, utils.get_consensus(block))
    if args.level is not None:
        block["level"] = args.level
    if args.name is not None:
        block["name"] = args.name
    if args.iron:
        for sequence in block["sequences"]:
            sequence["seq"] = cma.iron(sequence["seq"])
    return block
