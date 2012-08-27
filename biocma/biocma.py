"""Use CMA blocks as a subclass of Biopython's MSA.

This module converts the raw dictionaries produced by the CMA parser into
Biopython-compatible alignments.

See:
    Bio.AlignIO
    Bio.Align.MultipleSeqAlignment
"""
# ENH: column_annotations, like SeqRecord letter_annotations but for whole
# column positions -- use for conservation & pattern info

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import extended_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import cma
from . import sugar


def parse(infile):
    for block in cma.parse(infile):
        yield ChainMultiAlignment.from_block(block)

read = sugar.make_reader(parse)


def write(msas, outfile, do_iron=True):
    # ENH: keep the other block info through round-trip
    blocks = [cma.collapse_to_consensus(msa, do_iron=do_iron) for msa in msas]
    cma.write(blocks, outfile)


def convert(in_file, in_fmt, out_file, out_fmt):
    records = list(parse(in_file)
            if in_fmt == 'cma'
            else AlignIO.parse(in_file, in_fmt))
    if out_fmt == 'cma':
        write(records, out_file)
    else:
        AlignIO.write(records, out_file, out_fmt)


class ChainMultiAlignment(MultipleSeqAlignment):
    """Based on Biopython's class for a multiple sequence alignment."""

    def __init__(self, records,
            # CMA-specific attributes here
            # one=1,
            level=0,
            name=None,
            params=None,
            query_length=None,
            query_chars=None,
            ):
        # NB: alphabet is always protein; X is OK
        MultipleSeqAlignment.__init__(self, records, extended_protein)
        self.level = level
        self.name = name
        self.params = params
        self.query_length = query_length
        self.query_chars = query_chars

    @classmethod
    def from_block(cls, block):
        """Instantiate this class given a raw block (see parse_raw)."""
        rseqs = cma.realign_seqs(block)
        records = (SeqRecord(Seq(rseq, extended_protein),
                    id=bseq['id'],
                    description=bseq['description'],
                    dbxrefs=bseq['dbxrefs'].values(),   # list of strings
                    annotations=dict(
                        index=bseq['index'],
                        length=bseq['length'],
                        dbxrefs=bseq['dbxrefs'],
                        phylum=bseq['phylum'],
                        taxchar=bseq['taxchar'],
                        head_seq=bseq['head_seq'],
                        tail_seq=bseq['tail_seq'],
                        head_len=bseq['head_len'],
                        tail_len=bseq['tail_len'],      # dict
                        ),
                    # ENH: annotate with conservation levels
                    # letter_annotations=bseq['x'],
                    )
                for bseq, rseq in zip(block['sequences'], rseqs))
        return cls(records,
                # CMA attributes
                # block['one'],
                block['level'],
                block['name'],
                block['params'],
                block['query_length'],
                block['query_chars'],
                )

