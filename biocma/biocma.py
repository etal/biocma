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
import copy
import io
import logging
from xml.etree.ElementTree import ParseError
from xml.sax import SAXParseException

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

    A format that Biopython reads or writes as bytes, such as seqxml, is read
    from or written to the byte stream under a text handle, or a path opened
    in binary mode. A text handle given for such input must not have been
    read from, since its text layer may hold data read ahead from the byte
    stream. A file in such a format holds one alignment: Biopython would
    write a complete document for each, which together are not a valid file,
    so more than one alignment raises ValueError before the output is opened.
    Malformed input in an XML format raises ValueError, as it does in other
    formats.
    """
    if in_fmt == "cma":
        records = list(parse(in_file))
    else:
        records = _read_alignments(in_file, in_fmt)
    if out_fmt == "cma":
        write(records, out_file)
    elif _writes_bytes(out_fmt):
        if len(records) > 1:
            raise ValueError(
                f"A {out_fmt} file holds one alignment; the input has {len(records)}"
            )
        with sugar.maybe_open(out_file, "wb") as stream:
            AlignIO.write(records, stream, out_fmt)
    else:
        AlignIO.write(records, out_file, out_fmt)


def _writes_bytes(fmt):
    """Return whether Bio.AlignIO writes fmt as bytes rather than text.

    Biopython checks the handle's mode when it constructs a writer, before
    the writer looks at any alignment, so an empty alignment written to a
    throwaway text handle raises StreamModeError for exactly these formats.
    Text formats either accept the empty alignment or reject it with another
    ValueError, as does an unknown format, which the real write reports.
    The probe uses a throwaway handle, rather than a first attempt on the
    output as for input, because Bio.AlignIO opens a path in text mode, and
    so truncates it, before it checks the mode.
    """
    try:
        AlignIO.write([MultipleSeqAlignment([])], io.StringIO(), fmt)
    except StreamModeError:
        return True
    except ValueError:
        pass
    return False


def _read_alignments(in_file, fmt):
    """Read the alignments in a Bio.AlignIO format into a list."""
    try:
        try:
            return list(AlignIO.parse(in_file, fmt))
        except StreamModeError:
            # Biopython checks the handle's mode before it reads anything, so
            # the read can be repeated on a binary handle; a binary handle
            # given a text format raises the same error again
            with sugar.maybe_open(in_file, "rb") as stream:
                return list(AlignIO.parse(stream, fmt))
    except (SAXParseException, ParseError) as exc:
        # The SeqXML and UniProt XML parsers raise these, which are not
        # ValueError
        raise ValueError(f"Malformed {fmt} input: {exc}") from exc


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
        super().__init__(records)
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
            level=block["level"],
            name=block["name"],
            params=block["params"],
            query_length=block["query_length"],
            query_chars=block["query_chars"],
        )

    def _dots_as_dashes(self):
        """Copy this alignment as a MultipleSeqAlignment with '-' for '.'.

        Biopython reads only '-' as a gap, so the insert gaps that from_block
        writes as '.' would otherwise count as residues. The records are deep
        copies, since Biopython's conversion empties the letter annotations of
        a shallow copy in place, and so those of this alignment's records.
        """
        records = []
        for record in self:
            record = copy.deepcopy(record)
            record.seq = Seq(str(record.seq).replace(".", "-"))
            records.append(record)
        return MultipleSeqAlignment(
            records,
            annotations=self.annotations,
            column_annotations=self.column_annotations,
        )

    @property
    def alignment(self):
        """Return a Bio.Align.Alignment, reading '.' as a gap character.

        A column with a gap in every record is dropped by Biopython's
        conversion, so the result can have fewer columns than
        get_alignment_length().
        """
        return self._dots_as_dashes().alignment

    @property
    def substitutions(self):
        """Return the substitution counts, reading '.' as a gap character."""
        return self._dots_as_dashes().substitutions
