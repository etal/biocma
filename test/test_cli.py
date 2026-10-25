"""Tests for the command-line programs in biocma.cli."""

import gc
import io
import os
import subprocess
import sys

import pytest
from Bio import AlignIO

from biocma import biocma, cli, cma

COMMANDS = [
    ("cma-blocks", cli.blocks_main, []),
    ("cma-convert", cli.convert_main, ["-o", "fasta"]),
    ("cma-dedupe", cli.dedupe_main, []),
    ("cma-tweak", cli.tweak_main, []),
]


def read_alignment(text):
    """Parse FASTA output; Biopython rejects rows of unequal length."""
    return AlignIO.read(io.StringIO(text), "fasta")


@pytest.mark.parametrize("prog, main", [command[:2] for command in COMMANDS])
def test_help(prog, main, capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.startswith(f"usage: {prog} ")


@pytest.mark.parametrize("prog, main, options", COMMANDS)
def test_missing_file_is_reported(prog, main, options, tmp_path, capsys):
    missing = str(tmp_path / "missing.cma")
    assert main(options + [missing]) == 1
    out, err = capsys.readouterr()
    assert out == ""
    assert err.startswith(f"{prog}: error: ")
    assert missing in err


@pytest.mark.parametrize("prog, main, options", COMMANDS)
def test_malformed_input_is_reported(prog, main, options, fikk_text, tmp_path, capsys):
    truncated = tmp_path / "truncated.cma"
    truncated.write_text(fikk_text[: len(fikk_text) // 2])
    assert main(options + [str(truncated)]) == 1
    out, err = capsys.readouterr()
    assert out == ""
    assert err.startswith(f"{prog}: error: Malformed CMA block ")
    gc.collect()


class ClosedPipe(io.StringIO):
    """Standard output whose reader has exited, at a write or at the flush."""

    def __init__(self, fails_at):
        super().__init__()
        self.fails_at = fails_at

    def write(self, text):
        if self.fails_at == "write":
            raise BrokenPipeError(32, "Broken pipe")
        return super().write(text)

    def flush(self):
        if self.fails_at == "flush":
            raise BrokenPipeError(32, "Broken pipe")


@pytest.mark.parametrize("fails_at", ["write", "flush"])
@pytest.mark.parametrize("prog, main, options", COMMANDS)
def test_closed_pipe_exits_quietly(
    prog, main, options, fails_at, fikk_cma, monkeypatch, capsys
):
    monkeypatch.setattr(sys, "stdout", ClosedPipe(fails_at))
    assert main(options + [fikk_cma]) == 1
    err = capsys.readouterr().err
    # cma-dedupe reports its count before writing
    assert err in ("", "Dropped 0 duplicate sequences\n")


@pytest.mark.parametrize(
    "command", [["blocks"], ["convert", "-o", "seqxml"]], ids=["blocks", "seqxml"]
)
def test_closed_pipe_in_a_real_process(command, fikk_cma, tmp_path, monkeypatch):
    # Standard output is a pipe whose read end is already closed, so the first
    # write fails; the process exits with status 1 and prints nothing. Output
    # still buffered must not be flushed again at interpreter exit, which
    # would print 'Exception ignored ... BrokenPipeError' and exit with status
    # 120. PYTHONUNBUFFERED would hide that case, so the child runs without it.
    # Two records keep the output smaller than the buffer, also as seqxml,
    # which is written to the byte stream under sys.stdout.
    block = cma.read(fikk_cma)
    small = tmp_path / "small.cma"
    cma.write(dict(block, sequences=block["sequences"][:2]), small)
    monkeypatch.delenv("PYTHONUNBUFFERED", raising=False)
    read_end, write_end = os.pipe()
    os.close(read_end)
    name, *options = command
    script = (
        f"import sys; from biocma import cli; sys.exit(cli.{name}_main(sys.argv[1:]))"
    )
    with os.fdopen(write_end, "wb") as stdout:
        result = subprocess.run(
            [sys.executable, "-c", script, *options, str(small)],
            stdout=stdout,
            stderr=subprocess.PIPE,
            text=True,
        )
    assert (result.returncode, result.stderr) == (1, "")


def test_tweak_reports_invalid_block(fikk_cma, tmp_path, capsys):
    block = cma.read(fikk_cma)
    # One aligned column more than the query length
    block["sequences"][1]["seq"] += "A"
    ragged = tmp_path / "ragged.cma"
    cma.write(block, ragged)
    assert cli.tweak_main([str(ragged), "-c"]) == 1
    err = capsys.readouterr().err
    assert err.startswith("cma-tweak: error: Sequence ")
    assert err.endswith("query length is 255\n")
    gc.collect()


def test_blocks_writes_aligned_columns(fikk_cma, capsys):
    assert cli.blocks_main([fikk_cma]) == 0
    alignment = read_alignment(capsys.readouterr().out)
    assert len(alignment) == 24
    assert alignment.get_alignment_length() == 255
    assert alignment[0].description == "FIKK consensus"
    assert alignment[1].description.startswith(
        "psu|MAL7P1.127 | organism=Plasmodium_falciparum_3D7"
    )


@pytest.mark.parametrize("args", [[], ["-"]], ids=["default", "dash"])
def test_blocks_reads_standard_input(args, fikk_cma, monkeypatch, capsys):
    cli.blocks_main([fikk_cma])
    expected = capsys.readouterr().out
    with open(fikk_cma) as handle:
        monkeypatch.setattr("sys.stdin", handle)
        assert cli.blocks_main(args) == 0
    assert capsys.readouterr().out == expected


def test_blocks_excludes_flanks_across_blocks(aaaplus_sma, capsys):
    assert cli.blocks_main([aaaplus_sma]) == 0
    alignment = read_alignment(capsys.readouterr().out)
    # 153 records in 13 blocks, less 4 whose bodies are wrapped over lines
    assert len(alignment) == 149
    assert alignment.get_alignment_length() == 209
    # The body of RfcB_human begins with an N-terminal flank, "{(HY)ELPW...";
    # the flank is not part of the alignment.
    rfcb = next(rec for rec in alignment if rec.id == "RfcB_human")
    assert str(rfcb.seq).startswith("ELPWVEKYRP")


def test_convert_to_and_from_fasta(fikk_cma, tmp_path, capsys):
    ids = [seq["id"] for seq in cma.read(fikk_cma)["sequences"]]
    assert cli.convert_main(["-o", "fasta", fikk_cma]) == 0
    fasta = capsys.readouterr().out
    assert [rec.id for rec in read_alignment(fasta)] == ids
    aligned = tmp_path / "fikk.fasta"
    aligned.write_text(fasta)

    assert cli.convert_main(["-i", "fasta", str(aligned)]) == 0
    converted = cma.read(io.StringIO(capsys.readouterr().out))
    assert [seq["id"] for seq in converted["sequences"]] == ids


def test_convert_to_binary_format(fikk_cma, capsysbinary):
    # Biopython writes seqxml only to a binary handle
    assert cli.convert_main(["-o", "seqxml", fikk_cma]) == 0
    out = capsysbinary.readouterr().out
    records = AlignIO.read(io.BytesIO(out), "seqxml")
    assert [rec.id for rec in records] == [
        seq["id"] for seq in cma.read(fikk_cma)["sequences"]
    ]


def test_convert_refuses_several_alignments_in_binary_format(aaaplus_sma, capsys):
    # Biopython would write one seqxml document for each of the blocks
    assert cli.convert_main(["-o", "seqxml", aaaplus_sma]) == 1
    out, err = capsys.readouterr()
    assert out == ""
    assert err.startswith("cma-convert: error: A seqxml file holds one alignment")


def test_convert_from_binary_format_on_stdin(fikk_cma, tmp_path, monkeypatch, capsys):
    # Biopython reads seqxml only from a binary handle, here sys.stdin.buffer
    xml = tmp_path / "fikk.xml"
    biocma.convert(fikk_cma, "cma", xml, "seqxml")
    monkeypatch.setattr(sys, "stdin", io.TextIOWrapper(io.BytesIO(xml.read_bytes())))
    assert cli.convert_main(["-i", "seqxml", "-"]) == 0
    converted = cma.read(io.StringIO(capsys.readouterr().out))
    assert [seq["id"] for seq in converted["sequences"]] == [
        seq["id"] for seq in cma.read(fikk_cma)["sequences"]
    ]


def test_convert_reports_malformed_binary_input(fikk_cma, capsys):
    # The SeqXML parser raises SAXParseException, which is not a ValueError
    assert cli.convert_main(["-i", "seqxml", fikk_cma]) == 1
    out, err = capsys.readouterr()
    assert out == ""
    assert err.startswith("cma-convert: error: Malformed seqxml input: ")


@pytest.mark.parametrize("option", ["-o", "-i"])
def test_convert_reports_unknown_format(option, fikk_cma, capsys):
    assert cli.convert_main([option, "nosuchformat", fikk_cma]) == 1
    out, err = capsys.readouterr()
    assert out == ""
    assert err.startswith("cma-convert: error: ")
    assert "nosuchformat" in err


@pytest.mark.parametrize("options", [[], ["-o", "fasta", "-i", "fasta"]])
def test_convert_requires_one_direction(options, fikk_cma):
    with pytest.raises(SystemExit) as exc:
        cli.convert_main(options + [fikk_cma])
    assert exc.value.code == 2


def test_dedupe_keeps_best_and_last_of_ties(fikk_cma, tmp_path, capsys):
    block = cma.read(fikk_cma)
    first, second, third, *rest = block["sequences"]
    worse = dict(second, seq=second["seq"].lower(), description="worse")
    tie = dict(third, description="tie")
    block["sequences"] = [first, worse, second, third, tie] + rest
    infile = str(tmp_path / "dupes.cma")
    cma.write(block, infile)

    assert cli.dedupe_main([infile]) == 0
    out, err = capsys.readouterr()
    assert err == "Dropped 2 duplicate sequences\n"
    kept = cma.read(io.StringIO(out))["sequences"]
    assert len(kept) == 24
    assert [(seq["id"], seq["description"]) for seq in kept[:3]] == [
        (first["id"], first["description"]),
        (second["id"], second["description"]),
        (third["id"], "tie"),
    ]


def test_dedupe_compares_ids_within_each_block(fikk_cma, tmp_path, capsys):
    block = cma.read(fikk_cma)
    # Both blocks hold the same 24 IDs; only the first repeats one of them.
    repeated = dict(
        block, name="repeated", sequences=block["sequences"] + block["sequences"][1:2]
    )
    infile = str(tmp_path / "blocks.cma")
    cma.write([repeated, block], infile)

    assert cli.dedupe_main([infile]) == 0
    out, err = capsys.readouterr()
    assert err == "Dropped 1 duplicate sequences\n"
    blocks = list(cma.parse(io.StringIO(out)))
    assert [len(b["sequences"]) for b in blocks] == [24, 24]


def test_dedupe_reports_write_error(fikk_cma, monkeypatch, capsys):
    def fail(blocks, outfile):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(cma, "write", fail)
    assert cli.dedupe_main([fikk_cma]) == 1
    assert capsys.readouterr().err.endswith(
        "cma-dedupe: error: [Errno 28] No space left on device\n"
    )


def test_tweak_sets_level_and_name(fikk_cma, tmp_path):
    raised = str(tmp_path / "raised.cma")
    assert cli.tweak_main([fikk_cma, "--level", "10", "-o", raised]) == 0
    assert cma.read(raised)["level"] == 10
    # Level 0 must also be settable, to reset a non-zero level.
    reset = str(tmp_path / "reset.cma")
    assert (
        cli.tweak_main([raised, "--level", "0", "--name", "Testfam", "-o", reset]) == 0
    )
    block = cma.read(reset)
    assert (block["level"], block["name"]) == (0, "Testfam")


@pytest.mark.parametrize("name", ["Foo(1)", "two\nlines", ""])
def test_tweak_rejects_unreadable_name(name, fikk_cma, capsys):
    with pytest.raises(SystemExit) as exc:
        cli.tweak_main([fikk_cma, "--name", name])
    assert exc.value.code == 2
    assert "--name" in capsys.readouterr().err


def test_tweak_adds_consensus(fikk_cma, capsys):
    assert cli.tweak_main([fikk_cma, "-c"]) == 0
    sequences = cma.read(io.StringIO(capsys.readouterr().out))["sequences"]
    assert len(sequences) == 25
    assert sequences[0]["id"] == "consensus"


def test_tweak_irons_indel_transitions(fikk_cma, tmp_path, capsys):
    block = cma.read(fikk_cma)
    # The example in cma.iron's docstring: 13 residues in 12 aligned columns.
    # test/fikk.cma has nothing to iron.
    bent = dict(
        block["sequences"][0],
        id="bent",
        description="",
        length=13,
        seq="AAAAbc--de-f--gAAA",
    )
    infile = str(tmp_path / "bent.cma")
    cma.write(
        dict(block, query_length=12, query_chars="*" * 12, sequences=[bent]), infile
    )

    assert cli.tweak_main([infile, "--iron"]) == 0
    ironed = cma.read(io.StringIO(capsys.readouterr().out))["sequences"]
    assert [seq["seq"] for seq in ironed] == ["AAAABCDEFgAAA"]
