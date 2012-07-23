#!/usr/bin/env python

"""Convert a CMA file to another alignment format.

Usage:
    cma-convert.py to_format [cma file] > outfile
"""
# TODO - auto-generate a list of supported AlignIO formats
#   - display these w/ '-h'

import sys
from Bio import AlignIO
from biocma import biocma

if len(sys.argv) == 2:
    infile = sys.stdin
    do_close = False
elif len(sys.argv) == 3:
    infile = open(sys.argv[2])
    do_close = True
else:
    sys.exit(__doc__)

outfmt = sys.argv[1]
records = biocma.parse(infile)
count = AlignIO.write(records, sys.stdout, outfmt)
print >>sys.stderr, "cma-convert.py: Converted %d alignment(s)" % count
if do_close:
    infile.close()

