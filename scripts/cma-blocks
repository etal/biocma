#!/usr/bin/env python

"""Remove unaligned 'insert' columns from a CMA file, print as FASTA alignment.

Usage:
    cma-blocks.py input.cma > output.fasta

"""

def parse_cmafile(cmafname):
    with open(cmafname) as infile:
        lines = iter(infile)
        for line in lines:
            if line.startswith('>'):
                header = line.strip()
                seq = next(lines).strip()
                yield (header, seq[3:-4])


if __name__ == '__main__':
    import sys
    if not len(sys.argv) == 2 or sys.argv[1] in ('-h', '--help'):
        sys.exit(__doc__)

    infname = sys.argv[1]
    for header, seq in parse_cmafile(infname):
        print header
        print ''.join((c for c in seq if not c.islower()))

