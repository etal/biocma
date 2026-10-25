"""Tests for the Biopython interface (biocma.biocma)."""

import io
import logging

import pytest
from Bio import AlignIO, StreamModeError
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
    assert aln[1].dbxrefs == ["psu:MAL7P1.127"]


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


def test_export_to_seqxml(fikk_cma):
    # The SeqXML writer requires each cross-reference in 'source:id' form
    AlignIO.write(biocma.read(fikk_cma), io.BytesIO(), "seqxml")


def test_convert_writes_binary_format_to_path(fikk_cma, tmp_path):
    # Bio.AlignIO opens a path in text mode, which seqxml rejects
    out = tmp_path / "fikk.xml"
    biocma.convert(fikk_cma, "cma", out, "seqxml")
    assert len(AlignIO.read(io.BytesIO(out.read_bytes()), "seqxml")) == 24


def test_convert_reads_binary_format(fikk_cma, tmp_path):
    # Bio.AlignIO opens a path in text mode, which seqxml rejects; Bio.SeqIO
    # reads the CMA fields back as one-element lists of strings
    xml = tmp_path / "fikk.xml"
    biocma.convert(fikk_cma, "cma", xml, "seqxml")
    out = io.StringIO()
    biocma.convert(xml, "seqxml", out, "cma")
    out.seek(0)
    expect = cma.read(io.StringIO(write_to_str([biocma.read(fikk_cma)])))
    assert cma.read(out)["sequences"] == expect["sequences"]


@pytest.mark.parametrize("fmt", ["seqxml", "uniprot-xml"])
@pytest.mark.parametrize("as_path", [False, True], ids=["handle", "path"])
def test_convert_reports_malformed_xml_input(fmt, as_path, tmp_path):
    # Biopython's XML parsers raise exceptions that are not ValueError. The
    # seqXML root lets the SeqXML reader reach the truncation, where its SAX
    # parser fails; a different root is reported by the reader as ValueError
    path = tmp_path / "truncated.xml"
    path.write_bytes(b'<?xml version="1.0"?>\n<seqXML seqXMLversion="0.4"><entry')
    in_file = path if as_path else io.BytesIO(path.read_bytes())
    with pytest.raises(ValueError, match=f"Malformed {fmt} input: "):
        biocma.convert(in_file, fmt, io.StringIO(), "cma")


def test_convert_refuses_several_alignments_in_binary_format(aaaplus_sma, tmp_path):
    out = tmp_path / "out.xml"
    out.write_text("earlier content")
    with pytest.raises(ValueError, match="A seqxml file holds one alignment"):
        biocma.convert(aaaplus_sma, "cma", out, "seqxml")
    # Refused before the path was opened for writing
    assert out.read_text() == "earlier content"


def test_convert_refuses_several_alignments_to_binary_handle(aaaplus_sma):
    # Biopython accepts a binary handle without error, and would write one
    # seqxml document for each block
    handle = io.BytesIO()
    with pytest.raises(ValueError, match="A seqxml file holds one alignment"):
        biocma.convert(aaaplus_sma, "cma", handle, "seqxml")
    assert handle.getvalue() == b""


@pytest.mark.parametrize(
    "handle_type, fmt", [(io.StringIO, "seqxml"), (io.BytesIO, "fasta")]
)
@pytest.mark.parametrize("direction", ["input", "output"])
def test_convert_reports_handle_of_wrong_mode(direction, handle_type, fmt, fikk_cma):
    # Biopython's StreamModeError, also where the handle has no byte stream
    # to repeat the read or write on
    with pytest.raises(StreamModeError):
        if direction == "output":
            biocma.convert(fikk_cma, "cma", handle_type(), fmt)
        else:
            biocma.convert(handle_type(), fmt, io.StringIO(), "cma")


def test_alignment_reads_dots_as_gaps(fikk_cma):
    aln = biocma.read(fikk_cma)
    rows = [str(rec.seq) for rec in aln]
    alignment = aln.alignment
    assert [str(rec.seq) for rec in alignment.sequences] == [
        row.replace("-", "").replace(".", "") for row in rows
    ]
    assert "." not in aln.substitutions.alphabet
    # The records of the ChainMultiAlignment keep their '.' gaps
    assert [str(rec.seq) for rec in aln] == rows


def test_alignment_leaves_letter_annotations_of_records():
    record = SeqRecord(Seq("AB.D"), id="a", letter_annotations={"q": [1, 2, 3, 4]})
    aln = biocma.ChainMultiAlignment([record, SeqRecord(Seq("ABcD"), id="b")])
    assert aln.alignment.sequences[0].letter_annotations == {"q": [1, 2, 4]}
    assert record.letter_annotations == {"q": [1, 2, 3, 4]}


def test_read_back_residues_past_last_consensus_column():
    # The write keeps the residues as a trailing insert of the second record
    records = [
        SeqRecord(Seq("ABCDE..."), id="cons"),
        SeqRecord(Seq("ABCDEfgh"), id="longer"),
        SeqRecord(Seq("ABCDE---"), id="short"),
    ]
    text = write_to_str([biocma.ChainMultiAlignment(records, name="t")])
    aln = biocma.read(io.StringIO(text))
    assert [str(rec.seq) for rec in aln] == ["ABCDE...", "ABCDEfgh", "ABCDE..."]


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
    assert (
        "\n>psu|MAL7P1.127 {|888(935)|}| organism=Plasmodium_falciparum_3D7 | " in text
    )


@pytest.mark.parametrize("fixture", ["fikk_cma", "aaaplus_sma"])
def test_write_keeps_cma_fields(fixture, request):
    fields = (
        "id",
        "length",
        "head_len",
        "tail_len",
        "head_seq",
        "tail_seq",
        "phylum",
        "taxchar",
    )
    path = request.getfixturevalue(fixture)
    original = [block for block in cma.parse(path) if block["sequences"]]
    written = list(cma.parse(io.StringIO(write_to_str(biocma.parse(path)))))
    assert [block["params"] for block in written] == [
        block["params"] for block in original
    ]
    assert [
        [tuple(seq[f] for f in fields) for seq in block["sequences"]]
        for block in written
    ] == [
        [tuple(seq[f] for f in fields) for seq in block["sequences"]]
        for block in original
    ]


@pytest.mark.parametrize(
    "query_chars, expect",
    [("!**", "!**"), ("!***", "***")],
    ids=["kept", "wrong-length"],
)
def test_write_keeps_block_params(query_chars, expect):
    aln = biocma.ChainMultiAlignment(
        [SeqRecord(Seq("ACD"), id="cons"), SeqRecord(Seq("AC-"), id="r1")],
        params="go=19,gx=2,pn=5.0,lf=0,rf=0",
        query_chars=query_chars,
    )
    block = cma.read(io.StringIO(write_to_str([aln])))
    assert (block["params"], block["query_chars"]) == (
        "go=19,gx=2,pn=5.0,lf=0,rf=0",
        expect,
    )


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
