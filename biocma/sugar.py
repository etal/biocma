"""Helpers."""

import contextlib
import os


def make_reader(parser):
    def read(infile, *args, **kwargs):
        with contextlib.closing(parser(infile, *args, **kwargs)) as gen:
            try:
                first = next(gen)
            except StopIteration:
                raise ValueError("Input file is empty") from None
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
    """
    if isinstance(infile, (str, bytes, os.PathLike)):
        with open(infile, mode) as handle:
            yield handle
    else:
        yield infile


def unblank(stream):
    """Remove blank lines from a file being read iteratively."""
    return filter(None, (line.strip() for line in stream))
