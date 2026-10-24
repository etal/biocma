"""Helpers."""

import contextlib


def make_reader(parser):
    def read(infile, *args, **kwargs):
        gen = parser(infile, *args, **kwargs)
        try:
            first = next(gen)
        except StopIteration:
            raise ValueError("Input file is empty") from None
        try:
            next(gen)
        except StopIteration:
            return first
        else:
            raise ValueError("Input file contains multiple elements; "
                    "use parse() instead.")
    return read


@contextlib.contextmanager
def maybe_open(infile, mode='r'):
    """Take a file name or a handle, and return a handle.

    Simplifies creating functions that automagically accept either a file name
    or an already opened file handle.
    """
    # ENH: Exception safety?
    if isinstance(infile, str):
        handle = open(infile, mode)
        do_close = True
    else:
        handle = infile
        do_close = False
    yield handle
    if do_close:
        handle.close()


def unblank(stream):
    """Remove blank lines from a file being read iteratively."""
    return filter(None, (line.strip() for line in stream))


