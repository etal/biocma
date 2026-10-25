"""Tests for the Biopython interface (biocma.biocma)."""

import io
import logging

import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocma import biocma, cma


def write_to_str(msas):
    buf = io.StringIO()
    biocma.write(msas, buf)
    return buf.getvalue()


def ungapped(aln):
    return [str(rec.seq).replace("-", "").replace(".", "").upper() for rec in aln]


def test_read(fikk_cma):
    aln = biocma.read(fikk_cma)
    assert isinstance(aln, biocma.ChainMultiAlignment)
    assert len(aln) == 24
    assert aln[0].id == "FIKK"
    assert (aln.level, aln.name, aln.query_length) == (0, "FIKK", 255)
    assert aln[1].dbxrefs == ["MAL7P1.127"]


def test_parse_skips_block_without_sequences(aaaplus_sma, caplog):
    # Block MiscRfc has only wrapped sequence lines, which cma.parse skips
    with caplog.at_level(logging.WARNING, logger="biocma.biocma"):
        alignments = list(biocma.parse(aaaplus_sma))
    assert len(alignments) == len(list(cma.parse(aaaplus_sma))) - 1
    assert "MiscRfc" not in [aln.name for aln in alignments]
    assert "Skipping block MiscRfc: no sequences" in caplog.text


def test_read_rejects_file_without_sequences(fikk_cma):
    buf = io.StringIO()
    cma.write(dict(cma.read(fikk_cma), sequences=[]), buf)
    buf.seek(0)
    with pytest.raises(ValueError, match="no elements"):
        biocma.read(buf)


def test_export_to_genbank(fikk_cma):
    # The GenBank writer requires molecule_type in each record's annotations
    AlignIO.write(biocma.read(fikk_cma), io.StringIO(), "genbank")


def test_round_trip_preserves_ids_and_residues(fikk_cma):
    aln = biocma.read(fikk_cma)
    again = biocma.read(io.StringIO(write_to_str([aln])))
    assert [rec.id for rec in again] == [rec.id for rec in aln]
    assert ungapped(again) == ungapped(aln)


@pytest.mark.parametrize("fixture", ["fikk_cma", "aaaplus_sma"])
def test_write_is_stable_after_first_pass(fixture, request):
    # The first write collapses inserts and deletions in consensus rows, so
    # only later passes reproduce the written text.
    first = write_to_str(biocma.parse(request.getfixturevalue(fixture)))
    second = write_to_str(biocma.parse(io.StringIO(first)))
    assert second == first


def test_write_keeps_descriptions(fikk_cma):
    text = write_to_str([biocma.read(fikk_cma)])
    assert "\n>FIKK consensus\n" in text
    assert "\n>psu|MAL7P1.127 | organism=Plasmodium_falciparum_3D7 | " in text


def test_write_does_not_modify_alignment(fikk_cma):
    aln = biocma.read(fikk_cma)
    before = [(rec.id, str(rec.seq)) for rec in aln]
    outputs = [write_to_str([aln]) for _ in range(2)]
    assert [(rec.id, str(rec.seq)) for rec in aln] == before
    assert outputs[0] == outputs[1]


def test_write_irons_indels_by_default():
    # Column 3 is an insert ('.' in the consensus); r1's insert 'x' there sits
    # beside its deletion in column 2, so ironing turns it into a match.
    aln = MultipleSeqAlignment(
        [SeqRecord(Seq("AB.CDE"), id="cons"), SeqRecord(Seq("A-xC-E"), id="r1")]
    )
    block = cma.read(io.StringIO(write_to_str([aln])))
    assert [seq["seq"] for seq in block["sequences"]] == ["ABCDE", "AXC-E"]
