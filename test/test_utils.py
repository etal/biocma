"""Tests for utility functions (biocma.utils)."""

import copy

import pytest
from Bio import SeqIO

from biocma import cma, utils


def make_block(*rows):
    return {
        "query_length": sum(not c.islower() for c in rows[0]),
        "sequences": [
            {"id": f"s{i}", "seq": row, "head_len": 0} for i, row in enumerate(rows)
        ],
    }


def with_consensus(block):
    block = copy.deepcopy(block)
    block["sequences"].insert(0, utils.get_consensus(block))
    return block


def test_get_inserts(fikk_cma, fikk_fasta):
    block = cma.read(fikk_cma)
    inserts = utils.get_inserts(block)
    assert len(inserts) == len(block["sequences"])
    with open(fikk_fasta) as handle:
        fullseqs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    for sequence, targets in (
        (block["sequences"][1], ["n", "t", "fyklyllkkydsntlfnv"]),
        (
            block["sequences"][-1],
            [
                "altkl",
                "nkl",
                "siptvgfskdgdrlqemykasvcsyteecqg",
                "ndndgeylldge",
                "eh",
                "p",
                "epecancneedknmsennhkkdskhkgdsnhksdsnhksdsnhksdsnhksgsnhksdcnhksgsnhksdsnhqsdcnhmsdhnhksdnnhksdsshksdsshksdsshksgsnhksdnnhksdsshksgsnhksdhnhksdsnhksdsnhknesnhknesnhknesnhknesnhknesnhkndsnhksdsnhmsdhnhksdnnhksdhnhmsdhnhksdnnhksdnnhmsdhnhksdnnhksdnnhksdnnhksdhnhmsdhnhksdnnhksdhnhksdsnhmsdhnhmsdhnhksdhnhksdhnhksdnnhksdsnhksdsnhksdhnhksdsnhmsdhnhmsdhnhksdhnhksdnnhksdsnhksdsnhksdhnhksdsnhmsdhnhmsdhnhmsdhnhksdhnhksdnnhksdsnhksdsnhksdsnhksdhnhksdhkhmsdnnhksdnnhksdhnhksdnnhksdhnhksdsnhksdsnhksdsnhksdsnhksdnnhksdhnhnsdsnhmsdhnhksdhnhksdhnhksdnnhksdnnhksdhnhksdhkknnnnnkdnknddnddsdasdavhediellesysdlnkfnemlteqln",
                "vt",
                "edtrv",
                "pmythnl",
                "g",
                "sfqscqpcv",
                "iirehiklkidnpfehlstitdqee",
                "yfd",
                "ra",
                "fqlak",
            ],
        ),
    ):
        full = fullseqs[sequence["id"]]
        ins_ranges = [
            str(full.seq)[start - 1 : end] for start, end in inserts[sequence["id"]]
        ]
        assert [ins.lower() for ins in ins_ranges] == targets


def test_get_consensus_breaks_ties_by_first_residue_seen(fikk_cma):
    assert utils.get_consensus(make_block("AD", "CD"))["seq"] == "AD"
    assert utils.get_consensus(make_block("CD", "AD"))["seq"] == "CD"
    # Ties in the fixture: column 9 is S:4, N:4, C:4; column 23 is T:6, K:6
    cons = utils.get_consensus(cma.read(fikk_cma))["seq"]
    assert (cons[8], cons[22]) == ("S", "T")


def test_get_consensus_skips_inserts_and_marks_gaps():
    # Counting the leading inserts would make column 1 'q'
    cons = utils.get_consensus(make_block("A-", "qC-", "qC-"))
    assert (cons["seq"], cons["length"]) == ("CX", 2)


@pytest.mark.parametrize("rows", [("AC", "ACD"), ("AC", "A")], ids=["long", "short"])
def test_get_consensus_rejects_rows_of_another_length(rows):
    with pytest.raises(ValueError, match="s1 has [13] aligned characters"):
        utils.get_consensus(make_block(*rows))


def test_get_consensus_rejects_empty_block():
    with pytest.raises(ValueError, match="no sequences"):
        utils.get_consensus({"query_length": 2, "sequences": []})


def test_get_conservation():
    block = with_consensus(make_block("AC", "A-", "GgC", "AD"))
    assert block["sequences"][0]["seq"] == "AC"
    assert utils.get_conservation(block) == {1: 0.75, 2: 0.5}


def test_get_equivalent_positions():
    block = make_block("ABE", "AxyB-", "-BE")
    block["sequences"][1]["head_len"] = 10
    assert utils.get_equivalent_positions(block) == {
        1: {"s1": 11},
        2: {"s1": 14, "s2": 1},
        3: {"s2": 2},
    }


def test_equivalent_positions_agree_with_number_letters(fikk_cma):
    block = with_consensus(cma.read(fikk_cma))
    equivalent = utils.get_equivalent_positions(block)
    assert len(equivalent) == block["query_length"]
    for rec in block["sequences"][1:]:
        letters = utils.number_letters(rec)
        aligned = [c for c in rec["seq"] if not c.islower()]
        for posn, resnums in equivalent.items():
            if aligned[posn - 1] == "-":
                assert rec["id"] not in resnums
            else:
                assert letters[resnums[rec["id"]]] == aligned[posn - 1]
