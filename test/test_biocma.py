"""Tests for the Biopython interface (biocma.biocma)."""

import io

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocma import biocma, cma


def round_trip(aln):
    buf = io.StringIO()
    biocma.write([aln], buf)
    buf.seek(0)
    return biocma.read(buf)


def ungapped(aln):
    return [str(rec.seq).replace("-", "").replace(".", "").upper()
            for rec in aln]


def test_read(fikk_cma):
    aln = biocma.read(fikk_cma)
    assert isinstance(aln, biocma.ChainMultiAlignment)
    assert len(aln) == 24
    assert aln[0].id == "FIKK"
    assert (aln.level, aln.name, aln.query_length) == (0, "FIKK", 255)
    assert aln[1].dbxrefs == ["MAL7P1.127"]


def test_export_to_genbank(fikk_cma):
    # The GenBank writer requires molecule_type in each record's annotations
    AlignIO.write(biocma.read(fikk_cma), io.StringIO(), "genbank")


def test_round_trip_preserves_ids_and_residues(fikk_cma):
    aln = biocma.read(fikk_cma)
    again = round_trip(aln)
    assert [rec.id for rec in again] == [rec.id for rec in aln]
    assert ungapped(again) == ungapped(aln)


def test_round_trip_is_stable_after_first_pass(fikk_cma):
    # The fixture's consensus row has inserts and deletions, which the first
    # write collapses, so only later passes reproduce the aligned sequences.
    first = round_trip(biocma.read(fikk_cma))
    second = round_trip(first)
    assert [str(rec.seq) for rec in second] == [str(rec.seq) for rec in first]


def test_write_does_not_modify_alignment(fikk_cma):
    aln = biocma.read(fikk_cma)
    before = [(rec.id, str(rec.seq)) for rec in aln]
    outputs = []
    for _ in range(2):
        buf = io.StringIO()
        biocma.write([aln], buf)
        outputs.append(buf.getvalue())
    assert [(rec.id, str(rec.seq)) for rec in aln] == before
    assert outputs[0] == outputs[1]


def test_write_irons_indels_by_default():
    # Column 3 is an insert ('.' in the consensus); r1's insert 'x' there sits
    # beside its deletion in column 2, so ironing turns it into a match.
    aln = MultipleSeqAlignment([SeqRecord(Seq("AB.CDE"), id="cons"),
                                SeqRecord(Seq("A-xC-E"), id="r1")])
    buf = io.StringIO()
    biocma.write([aln], buf)
    block = cma.read(io.StringIO(buf.getvalue()))
    assert [seq["seq"] for seq in block["sequences"]] == ["ABCDE", "AXC-E"]
