"""Tests for parsing and writing the CMA format (biocma.cma)."""

import copy
import gc
import io
import logging
from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocma import cma


def write_to_str(blocks):
    buf = io.StringIO()
    cma.write(blocks, buf)
    return buf.getvalue()


def without_index(blocks):
    """Copy blocks without sequence 'index', which written files renumber."""
    blocks = copy.deepcopy(blocks)
    for block in blocks:
        for seq in block["sequences"]:
            del seq["index"]
    return blocks


def test_read(fikk_cma):
    block = cma.read(fikk_cma)
    assert len(block["sequences"]) == 24
    assert block["query_length"] == block["sequences"][0]["length"] == 255


def test_read_accepts_path_objects(fikk_cma):
    assert cma.read(Path(fikk_cma)) == cma.read(fikk_cma)


# An unclosed file is reported as a ResourceWarning when it is collected, and
# the test configuration turns that warning into an error.
def test_read_closes_file_of_several_blocks(aaaplus_sma):
    with pytest.raises(ValueError, match="multiple elements"):
        cma.read(aaaplus_sma)
    gc.collect()


def test_parse_closes_file_when_stopped_early(aaaplus_sma):
    blocks = cma.parse(aaaplus_sma)
    next(blocks)
    blocks.close()
    gc.collect()


def truncation_offsets():
    """Offsets that cut test/fikk.cma inside each kind of line."""
    # Read here rather than through a fixture, since parametrize needs them
    text = (Path(__file__).parent / "fikk.cma").read_text()
    preheader = text.index("\n$2=") + 1
    return {
        "block header": 20,
        "query line": text.index("\n(") + 10,
        "after preheader": text.index("\n", preheader) + 1,
        "sequence line": text.index("{()", preheader) + 20,
        "before terminator": text.index("_0]."),
    }


TRUNCATION_OFFSETS = truncation_offsets()


@pytest.mark.parametrize(
    "offset", TRUNCATION_OFFSETS.values(), ids=TRUNCATION_OFFSETS.keys()
)
def test_parse_rejects_truncated_input(offset, fikk_text, tmp_path):
    truncated = tmp_path / "truncated.cma"
    truncated.write_text(fikk_text[:offset])
    with pytest.raises(ValueError, match=r"^Malformed CMA block '\[0_\(1\)=FIKK"):
        cma.read(truncated)
    gc.collect()


@pytest.mark.parametrize("kept", [1, 2, 3])
def test_parse_rejects_block_without_terminator(kept, fikk_text):
    # The block cut after the first 1, 2 or 3 lines of its second sequence,
    # then a whole block, whose header comes where the first block expects a
    # sequence's header line, its sequence line, or the next sequence
    offset = fikk_text.index("\n$2=") + 1
    for _ in range(kept):
        offset = fikk_text.index("\n", offset) + 1
    with pytest.raises(ValueError, match="does not end before the next block"):
        list(cma.parse(io.StringIO(fikk_text[:offset] + fikk_text)))


def test_parse_passes_encoding_errors_through(fikk_cma):
    with open(fikk_cma, "rb") as handle:
        data = handle.read()
    # A byte that is not valid UTF-8, inside a sequence ID past the first
    # 8192-byte chunk that TextIOWrapper decodes, so the error arises while a
    # block is being parsed
    offset = data.index(b"\n>", 8192) + 3
    data = data[:offset] + b"\xfc" + data[offset + 1 :]
    with pytest.raises(UnicodeDecodeError):
        cma.read(io.TextIOWrapper(io.BytesIO(data), encoding="utf-8"))


def test_parse_skips_sequence_line_without_braces(fikk_text, caplog):
    lines = fikk_text.splitlines(keepends=True)
    bodies = [i for i, line in enumerate(lines) if line.startswith("{()")]
    lines[bodies[1]] = "NOBRACES\n"
    with caplog.at_level(logging.WARNING):
        block = cma.read(io.StringIO("".join(lines)))
    assert len(block["sequences"]) == 23
    assert "Unparseable sequence" in caplog.text


def test_parse_skips_wrapped_sequences(aaaplus_sma, caplog):
    with caplog.at_level(logging.WARNING):
        blocks = list(cma.parse(aaaplus_sma))
    assert [len(b["sequences"]) for b in blocks] == [
        4,
        8,
        12,
        3,
        0,
        14,
        9,
        1,
        31,
        13,
        19,
        4,
        31,
    ]
    messages = [rec.getMessage() for rec in caplog.records]
    skipped = [m.split()[2] for m in messages if m.startswith("Unparseable sequence:")]
    assert skipped == ["RfcE_human", "RfcBCD", "RfcA", "RfcE"]
    assert "Expected 4 sequences in block RfcE, found 3" in messages
    assert "Expected 3 sequences in block MiscRfc, found 0" in messages
    assert {rec.name for rec in caplog.records} == {"biocma.cma"}


def test_parse_header_and_flank_fields(fikk_cma, aaaplus_sma):
    # Expected values are copied from the fixture text
    fikk = cma.read(fikk_cma)
    assert fikk["params"] == "go=19,gx=2,pn=5.0,lf=0,rf=0"
    seq = fikk["sequences"][1]
    assert seq["id"] == "psu|MAL7P1.127"
    assert seq["dbxrefs"] == {"psu": "MAL7P1.127"}
    assert (seq["head_len"], seq["tail_len"]) == (888, 935)
    assert seq["description"].startswith(
        "| organism=Plasmodium_falciparum_3D7 | product=conserved"
    )

    records = {
        seq["id"]: seq for block in cma.parse(aaaplus_sma) for seq in block["sequences"]
    }
    rfca = records["RfcA_plafa"]
    assert (rfca["index"], rfca["length"]) == (15, 232)
    assert (rfca["head_len"], rfca["tail_len"]) == (349, 323)
    assert (rfca["phylum"], rfca["taxchar"]) == ("Apicomplexa", "E")
    assert rfca["description"] == (
        "gi|11999114|gb|AAG43050.1|  replication "
        "factor C subunit 1 [Plasmodium falciparum]"
    )
    rfcb = records["RfcB_human"]
    assert (rfcb["phylum"], rfcb["taxchar"]) == ("Chordata", "m")
    # Flanking sequences in parentheses; their lengths are inferred
    assert rfcb["head_seq"] == "HY"
    assert rfcb["tail_seq"].startswith("INSENVFKVCDEP")
    assert rfcb["tail_seq"].endswith("ARLCQKTMA")
    assert (rfcb["head_len"], rfcb["tail_len"]) == (2, len(rfcb["tail_seq"]))
    assert records["gi|255554168|ref|XP_002518124.1|"]["dbxrefs"] == {
        "gi": "255554168",
        "ref": "XP_002518124.1",
    }


def test_write_reproduces_source_text(fikk_cma, fikk_text):
    written = write_to_str(cma.read(fikk_cma))
    # The writer puts one space after each sequence ID; the fixture has two
    assert written.split() == fikk_text.split()


@pytest.mark.parametrize("fixture", ["fikk_cma", "aaaplus_sma"])
def test_round_trip_preserves_every_field(fixture, request):
    blocks = list(cma.parse(request.getfixturevalue(fixture)))
    expect = without_index(blocks)
    reparsed = list(cma.parse(io.StringIO(write_to_str(blocks))))
    assert without_index(reparsed) == expect
    for block in reparsed:
        assert [seq["index"] for seq in block["sequences"]] == list(
            range(1, len(block["sequences"]) + 1)
        )


@pytest.mark.parametrize("level", [2, 10, -1])
def test_round_trip_synthetic_block(level):
    # Header values that neither fixture varies: level, the "one" field and a
    # query_chars line that is not all '*'
    block = {
        "level": level,
        "one": 3,
        "name": "synth",
        "params": "go=19,gx=2,pn=5.0,lf=0,rf=0",
        "query_length": 4,
        "query_chars": "*!**",
        "sequences": [
            {
                "index": 1,
                "id": "cons",
                "description": "consensus row",
                "dbxrefs": {},
                "phylum": "",
                "taxchar": "",
                "head_len": 0,
                "tail_len": 0,
                "head_seq": "",
                "tail_seq": "",
                "length": 4,
                "seq": "ACDE",
            },
            {
                "index": 2,
                "id": "sp|P12345",
                "description": "second row",
                "dbxrefs": {"sp": "P12345"},
                "phylum": "Chordata",
                "taxchar": "M",
                "head_len": 2,
                "tail_len": 3,
                "head_seq": "MK",
                "tail_seq": "GHW",
                "length": 4,
                "seq": "AcC-E",
            },
        ],
    }
    expect = copy.deepcopy(block)
    assert cma.read(io.StringIO(write_to_str(block))) == expect


@pytest.mark.parametrize("fixture", ["fikk_cma", "aaaplus_sma"])
def test_write_is_idempotent(fixture, request):
    first = write_to_str(list(cma.parse(request.getfixturevalue(fixture))))
    second = write_to_str(list(cma.parse(io.StringIO(first))))
    assert second == first


def test_write_does_not_modify_blocks(aaaplus_sma):
    # Sequence indices in this fixture do not start at 1 in every block, and
    # the writer renumbers them in its output
    blocks = list(cma.parse(aaaplus_sma))
    before = copy.deepcopy(blocks)
    first = write_to_str(blocks)
    assert blocks == before
    assert write_to_str(blocks) == first


@pytest.mark.parametrize(
    "rec_id, description, expect",
    [
        # Read from CMA: the description does not begin with the ID
        ("psu|X", "| organism=Plasmodium", "| organism=Plasmodium"),
        ("FIKK", "consensus", "consensus"),
        ("abc", "abcdef", "abcdef"),
        # Read by Bio.SeqIO from FASTA: the description begins with the ID
        ("abc", "abc rest of it", "rest of it"),
        ("abc", "abc\trest", "rest"),
        # Read by Bio.AlignIO from Clustal: the description is the ID
        ("abc", "abc", " "),
        # No description at all
        ("abc", "", " "),
    ],
)
def test_seqrecord2sequence_description(rec_id, description, expect):
    record = SeqRecord(Seq("ACD"), id=rec_id, description=description)
    assert cma.seqrecord2sequence(record, 3, 1)["description"] == expect


def test_id_only_description_writes_two_spaces():
    # Records read from Clustal keep the header layout of earlier releases,
    # which the MAPGAPS template fixtures of etal/fammer pin byte for byte
    block = cma.consensus2block(SeqRecord(Seq("ACD"), id="abc", description="abc"))
    assert "\n>abc  \n" in write_to_str(block)


def test_collapse_moves_leading_residues_to_head_flank():
    # Residues before the first consensus column extend the head flank: its
    # length always, its sequence only where the record has one. The length
    # field counts only residues in the flank and aligned sequences, so the
    # moved residues add to it only for the record with a head flank sequence.
    aln = MultipleSeqAlignment(
        [
            SeqRecord(Seq("--ACD"), id="cons"),
            SeqRecord(Seq("MKACD"), id="plain"),
            SeqRecord(Seq("MKA-D"), id="lengths", annotations={"head_len": 5}),
            SeqRecord(
                Seq("-KACD"), id="flank", annotations={"head_len": 2, "head_seq": "PQ"}
            ),
        ]
    )
    block = cma.collapse_to_consensus(aln)
    assert [
        (seq["head_len"], seq["head_seq"], seq["seq"], seq["length"])
        for seq in block["sequences"][1:]
    ] == [(2, "", "ACD", 3), (7, "", "A-D", 2), (3, "PQK", "ACD", 6)]


def test_seqrecord2sequence_treats_none_annotations_as_absent():
    keys = ("head_len", "tail_len", "head_seq", "tail_seq", "phylum", "taxchar")
    record = SeqRecord(Seq("ACD"), id="r1", annotations=dict.fromkeys(keys))
    seq = cma.seqrecord2sequence(record, 3, 1)
    assert [seq[key] for key in keys] == [0, 0, "", "", "", ""]


def test_realign_seqs_of_empty_block():
    assert cma.realign_seqs({"sequences": []}) == []


def test_realign_seqs_pads_trailing_inserts():
    # Inserts before the first, between and after the last aligned residues
    block = {
        "sequences": [
            {"id": "cons", "seq": "ABC"},
            {"id": "tail", "seq": "ABCde"},
            {"id": "head", "seq": "aABC"},
            {"id": "inner", "seq": "ABcC"},
        ]
    }
    assert cma.realign_seqs(block) == [".AB.C..", ".AB.Cde", "aAB.C..", ".ABcC.."]


def test_realign_seqs_rejects_unequal_aligned_lengths():
    block = {"sequences": [{"id": "a", "seq": "ABC"}, {"id": "b", "seq": "ABc"}]}
    with pytest.raises(ValueError, match="Sequence b has 2 aligned characters"):
        cma.realign_seqs(block)


def test_iron():
    assert cma.iron("AAAAbc--de-f--gAAA") == "AAAABCDEFgAAA"
    seq = (
        "-a-a-LLLGQPIFPGDSGVDQLVEIIKVLgtptre---qiremnpnyteFKFPQIK"
        "---ahpwtkvfrprtPPEAIALCSRLLEYTPTARLT-----PLEACAHSFF-"
    )
    ironed = cma.iron(seq)
    assert ironed != seq
    assert cma.iron(ironed) == ironed
