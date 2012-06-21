#!/usr/bin/env python

"""I/O and utilities for the Consensus Multiple Alignment (CMA) file format.

This format represents protein sequence alignments. It is used by a few tools by
Dr. Andrew F. Neuwald, notably CHAIN and MAPGAPS.

See:
    - http://chain.igs.umaryland.edu/
    - http://mapgaps.igs.umaryland.edu/
"""

from glob import glob
try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(  name='BioCMA',
        version='dev',
        description=__doc__,
        author='Eric Talevich',
        author_email='eric.talevich@gmail.com',
        url='http://github.com/etal/biocma',
        packages=['biocma'],
        scripts=glob('scripts/*.py'),
        )

