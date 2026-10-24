BioCMA
======

Consensus Multiple Alignment format, using Biopython alignments.

I/O support and relevant functionality for the Consensus Alignment Format (CMA).
This format represents protein sequence alignments. It is used by a few tools by
Dr. Andrew F. Neuwald, notably CHAIN_ and MAPGAPS_.

.. _CHAIN: http://chain.igs.umaryland.edu/
.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/

Biopython objects and conventions are used where possible.


Installation
------------

This is an ordinary Python package. You can install it from source with the
setup.py script::

    python setup.py build
    python setup.py install

Stable releases are uploaded to PyPI as well, so you can install BioFrills with
Python package managers::

    pip install biofrills

Or::

    easy_install biofrills

Some scripts that I find useful are in the scripts/ directory. By default these
are not installed, but you can include them by uncommenting the line in setup.py
that starts with ``scripts=glob``...

Alternatively, you can just copy those scripts into another directory in your
``$PATH``.


What can the CMA format do for me?
----------------------------------

- Like the A2M and Stockholm formats, alignments are shown with insertions as
  lowercase characters and deletions are dashes
- Like the A3M format, alignments are pairwise versus a profile (or "consensus"
  sequence), which also dictates which sites are indels. By compressing the
  insert columns, a large of alignment of many divergent (but related) sequences
  can be shown without filling it will mostly gap characters, as Stockholm can.
- Like Stockholm, but unlike A2M and A3M, more than one alignment can be
  contained in a single file.
- Typically, an ungapped consensus sequence will be included as the first
  sequence.
- A FASTA-like header contains additional, optional fields for the number of
  leading and trailing sites and NBCI taxonomy codes


Who uses CMA?
-------------

The CMA format appears to have been invented by Dr. Andrew Neuwald at the
University of Maryland, and is used in these programs:

- MAPGAPS: http://mapgaps.igs.umaryland.edu/
- CHAIN: http://chain.igs.umaryland.edu/


Should I use CMA in my own work?
--------------------------------

Unless you're working with MAPGAPS or CHAIN, no.

