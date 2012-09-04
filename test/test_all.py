"""Unit tests for BioCMA."""

import unittest
from cStringIO import StringIO

from biocma import biocma, cma

class IOTests(unittest.TestCase):
    """Tests for parsing and writing the CMA format."""
    def test_read(self):
        pass

    def test_parse(self):
        pass

class UtilTests(unittest.TestCase):
    """Tests for utility functions."""

    def test_iron(self):
        seq = '-a-a-LLLGQPIFPGDSGVDQLVEIIKVLgtptre---qiremnpnyteFKFPQIK' \
                '---ahpwtkvfrprtPPEAIALCSRLLEYTPTARLT-----PLEACAHSFF-'
        iseq = cma.iron(seq)
        self.assertEqual(len(seq.replace('-', '')),
                         len(iseq.replace('-', '')))



# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

