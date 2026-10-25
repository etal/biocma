"""Use CMA blocks as a subclass of Biopython's MSA.

This module converts the raw dictionaries produced by the CMA parser into
Biopython-compatible alignments.

See:
    Bio.AlignIO
    Bio.Align.MultipleSeqAlignment
"""
# ENH: column_annotations, like SeqRecord letter_annotations but for whole
# column positions -- use for conservation & pattern info

import contextlib
import logging

from Bio import AlignIO, StreamModeError
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import cma, sugar

logger = logging.getLogger(__name__)


def parse(infile):
    """Parse each block of a CMA file as a ChainMultiAlignment.

    A block without sequences is skipped with a warning, since an alignment
    of no sequences cannot be written in most formats.
    """
    with contextlib.closing(cma.parse(infile)) as blocks:
        for block in blocks:
            if not block["sequences"]:
                logger.warning("Skipping block %s: no sequences", block["name"])
                continue
            yield ChainMultiAlignment.from_block(block)


read = sugar.make_reader(parse)


def write(msas, outfile, do_iron=True):
    blocks = [cma.collapse_to_consensus(msa, do_iron=do_iron) for msa in msas]
    cma.write(blocks, outfile)


def convert(in_file, in_fmt, out_file, out_fmt):
    """Convert alignments between CMA and a Bio.AlignIO format.

    A format that Biopython writes as bytes, such as seqxml, is written to the
    byte stream under a text handle, or to a path opened in binary mode. Such
    a file holds one alignment: Biopython would write a complete document for
    each, which together are not a valid file, so in these two cases more
    than one alignment raises ValueError. After that error a path is left
    empty, since it has already been opened for writing. A binary handle is
    given to Biopython without this check, and several alignments are then
    written as one invalid file.
    """
    records = list(
        parse(in_file) if in_fmt == "cma" else AlignIO.parse(in_file, in_fmt)
    )
    if out_fmt == "cma":
        write(records, out_file)
        return
    try:
        AlignIO.write(records, out_file, out_fmt)
    except StreamModeError:
        # Biopython checks the handle's mode before it writes anything, so the
        # write can be repeated on a binary handle
        with sugar.maybe_open(out_file, "wb") as stream:
            if stream is out_file:
                # No byte stream to retry on, or a binary handle given a text
                # format: the error is Biopython's to report
                raise
            if len(records) > 1:
                raise ValueError(
                    f"A {out_fmt} file holds one alignment; "
                    f"the input has {len(records)}"
                ) from None
            AlignIO.write(records, stream, out_fmt)


class ChainMultiAlignment(MultipleSeqAlignment):
    """Based on Biopython's class for a multiple sequence alignment."""

    def __init__(
        self,
        records,
        # CMA-specific attributes here
        level=0,
        name=None,
        params=None,
        query_length=None,
        query_chars=None,
    ):
        MultipleSeqAlignment.__init__(self, records)
        self.level = level
        self.name = name
        self.params = params
        self.query_length = query_length
        self.query_chars = query_chars

    @classmethod
    def from_block(cls, block):
        """Instantiate this class given a raw block from cma.parse."""
        rseqs = cma.realign_seqs(block)
        records = (
            SeqRecord(
                Seq(rseq),
                id=bseq["id"],
                description=bseq["description"],
                dbxrefs=[f"{key}:{value}" for key, value in bseq["dbxrefs"].items()],
                annotations=dict(
                    # NB: always protein; X is OK
                    molecule_type="protein",
                    index=bseq["index"],
                    length=bseq["length"],
                    dbxrefs=bseq["dbxrefs"],
                    phylum=bseq["phylum"],
                    taxchar=bseq["taxchar"],
                    head_seq=bseq["head_seq"],
                    tail_seq=bseq["tail_seq"],
                    head_len=bseq["head_len"],
                    tail_len=bseq["tail_len"],
                ),
                # ENH: annotate with conservation levels
            )
            for bseq, rseq in zip(block["sequences"], rseqs, strict=True)
        )
        return cls(
            records,
            # CMA attributes
            block["level"],
            block["name"],
            block["params"],
            block["query_length"],
            block["query_chars"],
        )
