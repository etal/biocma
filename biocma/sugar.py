"""Helpers."""

import contextlib
import os


def make_reader(parser):
    def read(infile, *args, **kwargs):
        with contextlib.closing(parser(infile, *args, **kwargs)) as gen:
            try:
                first = next(gen)
            except StopIteration:
                raise ValueError("Input file contains no elements") from None
            try:
                next(gen)
            except StopIteration:
                return first
            else:
                raise ValueError(
                    "Input file contains multiple elements; use parse() instead."
                )

    return read


@contextlib.contextmanager
def maybe_open(infile, mode="r"):
    """Take a file name or a handle, and return a handle.

    Simplifies creating functions that accept either a file name or an already
    opened file handle. A file name is a str, bytes or os.PathLike path, as
    open() accepts; the file it names is closed on leaving the context, also
    on an exception. Anything else is taken to be a handle and left open.

    A binary mode given a text handle yields the byte stream under it, which
    Bio.AlignIO needs to read or write a format such as seqxml through
    standard input or output. The text layer is flushed first, so that text
    written through it precedes the bytes. For reading, nothing may have
    been read through the text layer, which can hold data read ahead from
    the byte stream that would then be skipped; this cannot be checked on a
    pipe. A handle with no byte stream, such as io.StringIO, is yielded as
    it is.
    """
    if isinstance(infile, (str, bytes, os.PathLike)):
        with open(infile, mode) as handle:
            yield handle
    elif "b" in mode:
        infile.flush()
        yield getattr(infile, "buffer", infile)
    else:
        yield infile


def unblank(stream):
    """Remove blank lines from a file being read iteratively."""
    return filter(None, (line.strip() for line in stream))
