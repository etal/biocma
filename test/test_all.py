"""Unit tests for BioCMA."""

import unittest
from cStringIO import StringIO

from Bio import SeqIO

from biocma import biocma, cma, utils

SINGLE_CMA = "fikk.cma"

class IOTests(unittest.TestCase):
    """Tests for parsing and writing the CMA format."""

    def test_read(self):
        block = cma.read(SINGLE_CMA)
        self.assertEqual(len(block['sequences']), 24)
        self.assertEqual(block['query_length'], block['sequences'][0]['length'])

    def test_parse(self):
        blocks = list(cma.parse(SINGLE_CMA))
        self.assertEqual(len(blocks), 1)


class BioTests(unittest.TestCase):
    """Tests for the Biopython wrapper."""

    def test_read(self):
        aln = biocma.read(SINGLE_CMA)


class UtilTests(unittest.TestCase):
    """Tests for utility functions."""

    def test_iron(self):
        seq = ('-a-a-LLLGQPIFPGDSGVDQLVEIIKVLgtptre---qiremnpnyteFKFPQIK'
               '---ahpwtkvfrprtPPEAIALCSRLLEYTPTARLT-----PLEACAHSFF-')
        iseq = cma.iron(seq)
        self.assertEqual(len(seq.replace('-', '')),
                         len(iseq.replace('-', '')))

    def test_get_inserts(self):
        block = cma.read(SINGLE_CMA)
        inserts = utils.get_inserts(block)
        self.assertEqual(len(inserts), len(block['sequences']))
        fullseqs = SeqIO.to_dict(SeqIO.parse("fikk-full.fasta", 'fasta'))
        for sequence, targets in (
            (block['sequences'][1], ['n', 't', 'fyklyllkkydsntlfnv']),
            (block['sequences'][-1], ['altkl', 'nkl',
                                      'siptvgfskdgdrlqemykasvcsyteecqg',
                                      'ndndgeylldge', 'eh', 'p',
                                      'epecancneedknmsennhkkdskhkgdsnhksdsnhksdsnhksdsnhksgsnhksdcnhksgsnhksdsnhqsdcnhmsdhnhksdnnhksdsshksdsshksdsshksgsnhksdnnhksdsshksgsnhksdhnhksdsnhksdsnhknesnhknesnhknesnhknesnhknesnhkndsnhksdsnhmsdhnhksdnnhksdhnhmsdhnhksdnnhksdnnhmsdhnhksdnnhksdnnhksdnnhksdhnhmsdhnhksdnnhksdhnhksdsnhmsdhnhmsdhnhksdhnhksdhnhksdnnhksdsnhksdsnhksdhnhksdsnhmsdhnhmsdhnhksdhnhksdnnhksdsnhksdsnhksdhnhksdsnhmsdhnhmsdhnhmsdhnhksdhnhksdnnhksdsnhksdsnhksdsnhksdhnhksdhkhmsdnnhksdnnhksdhnhksdnnhksdhnhksdsnhksdsnhksdsnhksdsnhksdnnhksdhnhnsdsnhmsdhnhksdhnhksdhnhksdnnhksdnnhksdhnhksdhkknnnnnkdnknddnddsdasdavhediellesysdlnkfnemlteqln',
                                      'vt', 'edtrv', 'pmythnl', 'g',
                                      'sfqscqpcv', 'iirehiklkidnpfehlstitdqee',
                                      'yfd', 'ra', 'fqlak'])):
            full = fullseqs[sequence['id']]
            ins_ranges = [str(full.seq)[start-1:end]
                       for start, end in inserts[sequence['id']]]
            print sequence['id'], ins_ranges
            self.assertEqual(len(ins_ranges), len(targets))
            for ins, tgt in zip(ins_ranges, targets):
                self.assertEqual(ins.lower(), tgt)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

