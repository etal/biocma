"""Tests for parsing and writing the CMA format (biocma.cma)."""

import copy
import io
import logging

import pytest

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


def test_parse_skips_wrapped_sequences(aaaplus_sma, caplog):
    with caplog.at_level(logging.WARNING):
        blocks = list(cma.parse(aaaplus_sma))
    assert [len(b["sequences"]) for b in blocks] == [
        4, 8, 12, 3, 0, 14, 9, 1, 31, 13, 19, 4, 31]
    messages = [rec.getMessage() for rec in caplog.records]
    skipped = [m.split()[2] for m in messages
               if m.startswith("Unparseable sequence:")]
    assert skipped == ["RfcE_human", "RfcBCD", "RfcA", "RfcE"]
    assert "Expected 4 sequences in block RfcE, found 3" in messages
    assert "Expected 3 sequences in block MiscRfc, found 0" in messages


def test_parse_header_and_flank_fields(fikk_cma, aaaplus_sma):
    # Expected values are copied from the fixture text
    fikk = cma.read(fikk_cma)
    assert fikk["params"] == "go=19,gx=2,pn=5.0,lf=0,rf=0"
    seq = fikk["sequences"][1]
    assert seq["id"] == "psu|MAL7P1.127"
    assert seq["dbxrefs"] == {"psu": "MAL7P1.127"}
    assert (seq["head_len"], seq["tail_len"]) == (888, 935)
    assert seq["description"].startswith(
        "| organism=Plasmodium_falciparum_3D7 | product=conserved")

    records = {seq["id"]: seq
               for block in cma.parse(aaaplus_sma)
               for seq in block["sequences"]}
    rfca = records["RfcA_plafa"]
    assert (rfca["index"], rfca["length"]) == (15, 232)
    assert (rfca["head_len"], rfca["tail_len"]) == (349, 323)
    assert (rfca["phylum"], rfca["taxchar"]) == ("Apicomplexa", "E")
    assert rfca["description"] == ("gi|11999114|gb|AAG43050.1|  replication "
                                   "factor C subunit 1 [Plasmodium falciparum]")
    rfcb = records["RfcB_human"]
    assert (rfcb["phylum"], rfcb["taxchar"]) == ("Chordata", "m")
    # Flanking sequences in parentheses; their lengths are inferred
    assert rfcb["head_seq"] == "HY"
    assert rfcb["tail_seq"].startswith("INSENVFKVCDEP")
    assert rfcb["tail_seq"].endswith("ARLCQKTMA")
    assert (rfcb["head_len"], rfcb["tail_len"]) == (2, len(rfcb["tail_seq"]))
    assert records["gi|255554168|ref|XP_002518124.1|"]["dbxrefs"] == {
        "gi": "255554168", "ref": "XP_002518124.1"}


def test_write_reproduces_source_text(fikk_cma):
    with open(fikk_cma) as handle:
        source = handle.read()
    written = write_to_str(cma.read(fikk_cma))
    # The writer puts one space after each sequence ID; the fixture has two
    assert written.split() == source.split()


@pytest.mark.parametrize("fixture", ["fikk_cma", "aaaplus_sma"])
def test_round_trip_preserves_every_field(fixture, request):
    blocks = list(cma.parse(request.getfixturevalue(fixture)))
    expect = without_index(blocks)
    reparsed = list(cma.parse(io.StringIO(write_to_str(blocks))))
    assert without_index(reparsed) == expect
    for block in reparsed:
        assert ([seq["index"] for seq in block["sequences"]]
                == list(range(1, len(block["sequences"]) + 1)))


def test_round_trip_synthetic_block():
    # Header values that neither fixture varies: level, the "one" field and a
    # query_chars line that is not all '*'
    block = {
        "level": 2, "one": 3, "name": "synth",
        "params": "go=19,gx=2,pn=5.0,lf=0,rf=0",
        "query_length": 4, "query_chars": "*!**",
        "sequences": [
            {"index": 1, "id": "cons", "description": "consensus row",
             "dbxrefs": {}, "phylum": "", "taxchar": "",
             "head_len": 0, "tail_len": 0, "head_seq": "", "tail_seq": "",
             "length": 4, "seq": "ACDE"},
            {"index": 2, "id": "sp|P12345", "description": "second row",
             "dbxrefs": {"sp": "P12345"}, "phylum": "Chordata", "taxchar": "M",
             "head_len": 2, "tail_len": 3, "head_seq": "MK", "tail_seq": "GHW",
             "length": 4, "seq": "AcC-E"},
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


def test_iron():
    assert cma.iron("AAAAbc--de-f--gAAA") == "AAAABCDEFgAAA"
    seq = ("-a-a-LLLGQPIFPGDSGVDQLVEIIKVLgtptre---qiremnpnyteFKFPQIK"
           "---ahpwtkvfrprtPPEAIALCSRLLEYTPTARLT-----PLEACAHSFF-")
    ironed = cma.iron(seq)
    assert ironed != seq
    assert cma.iron(ironed) == ironed
