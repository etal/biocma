BioCMA
======

Consensus Multiple Alignment format, using Biopython alignments.

I/O support and relevant functionality for the Consensus Alignment Format (CMA).
This format represents protein sequence alignments. It is used by a few tools by
Dr. Andrew F. Neuwald, notably MAPGAPS_ and CHAIN_.

.. _MAPGAPS: https://mapgaps.igs.umaryland.edu/
.. _CHAIN: https://pubmed.ncbi.nlm.nih.gov/17962021/

Biopython objects and conventions are used where possible.


Installation
------------

BioCMA requires Python 3.12 or later and Biopython 1.86 or later. Install it
from a source checkout with::

    pip install .

This version is not yet released. The latest release on PyPI, 0.2.1, requires
Python 2.


Commands
--------

Installing the package also installs four commands. Each reads the CMA or
alignment file it is given, or standard input if none is given or the file is
named ``-``, and writes its result to standard output; ``cma-tweak`` can
instead write to the file named by ``--output``. Run a command with ``--help``
to see its options. A file that cannot be read, or input that is not valid
CMA, is reported in one line on standard error, and the command exits with
status 1.

``cma-blocks``
    Remove the unaligned insert columns from a CMA file and write the remaining
    aligned columns as a FASTA alignment.

``cma-convert``
    Convert a CMA file to another alignment format (``-o FORMAT``), or an
    alignment in another format to CMA (``-i FORMAT``). Formats are named as in
    Biopython's ``Bio.AlignIO``, for example ``fasta``, ``clustal`` or
    ``stockholm``.

``cma-dedupe``
    Keep one record for each sequence ID, the one with the most aligned
    residues, and report the number of records removed on standard error.

``cma-tweak``
    Set the level or name of each block, insert a consensus sequence, or iron
    out indels.


What can the CMA format do for me?
----------------------------------

- Like the A2M and Stockholm formats, alignments are shown with insertions as
  lowercase characters and deletions are dashes.
- Like the A3M format, alignments are pairwise versus a profile (or "consensus"
  sequence), which also dictates which sites are indels. By compressing the
  insert columns, a large alignment of many divergent (but related) sequences
  can be shown without filling it with mostly gap characters, as Stockholm
  would.
- Like Stockholm, but unlike A2M and A3M, more than one alignment can be
  contained in a single file.
- Typically, an ungapped consensus sequence will be included as the first
  sequence.
- A FASTA-like header contains additional, optional fields for the number of
  leading and trailing sites and NCBI taxonomy codes.


Who uses CMA?
-------------

The CMA format appears to have been invented by Dr. Andrew Neuwald at the
University of Maryland, and is used in these programs:

- MAPGAPS, described in `Neuwald (2009), Bioinformatics 25:1869-1875
  <https://pubmed.ncbi.nlm.nih.gov/19505947/>`_ and distributed from
  https://mapgaps.igs.umaryland.edu/
- CHAIN, described in `Neuwald (2007), Trends in Biochemical Sciences
  32:487-493 <https://pubmed.ncbi.nlm.nih.gov/17962021/>`_. Its website,
  chain.igs.umaryland.edu, is no longer available; the Internet Archive holds
  a `copy from 2022
  <https://web.archive.org/web/20220418025310/http://chain.igs.umaryland.edu/>`_.


Should I use CMA in my own work?
--------------------------------

Unless you're working with MAPGAPS or CHAIN, no.


Development
-----------

Two development workflows are supported from the same ``pyproject.toml``.
Both install the package in editable mode together with the tools in the
``dev`` dependency group, pytest, ruff and mypy.

uv_ takes every package from PyPI, at the versions recorded in ``uv.lock``::

    uv sync
    uv run pytest
    uv run ruff check
    uv run ruff format
    uv run mypy

pixi_ takes Python, Biopython, pytest, ruff and mypy from conda-forge, at the
versions recorded in ``pixi.lock``. The ``default`` environment uses the latest
Python that conda-forge provides, and ``py312`` the oldest supported one::

    pixi run test
    pixi run -e py312 test
    pixi run lint
    pixi run fmt
    pixi run typecheck

The block and sequence records that ``biocma.cma`` reads and writes are plain
dicts, described by the ``TypedDict`` classes ``cma.Block`` and
``cma.Sequence``. mypy checks that every function annotated to return one sets
exactly the declared keys, and checks a value's type where the type of the
expression giving it is known; values from unannotated helpers are not
checked.

Code is formatted with ``ruff format``, and the commit that first applied it is
listed in ``.git-blame-ignore-revs``. To skip that commit in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

.. _uv: https://docs.astral.sh/uv/
.. _pixi: https://pixi.sh/
