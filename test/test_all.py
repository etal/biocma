"""Unit tests for BioCMA."""

import unittest
from cStringIO import StringIO

from biocma import cma
from biocma import biocma

class IOTests(unittest.TestCase):
    """Tests for parsing and writing the CMA format."""
    def test_read(self):
        pass

    def test_parse(self):
        pass

class UtilTests(unittest.TestCase):
    """Tests for utility functions."""


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

