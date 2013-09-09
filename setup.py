#!/usr/bin/env python

"""I/O and utilities for the Consensus Multiple Alignment (CMA) file format."""

from glob import glob
from os.path import join
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='biocma',
    version='0.2.1',
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@gmail.com',
    url='http://github.com/etal/biocma',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['biocma'],
    #scripts=glob(join('scripts', '*')),
)

