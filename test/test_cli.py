"""Tests for the command-line programs in biocma.cli."""

import pytest
from Bio import AlignIO

from biocma import cli, cma


def read_fasta(text):
    """Return (header, sequence) pairs from FASTA text with one-line bodies."""
    lines = text.splitlines()
    assert all(line.startswith('>') for line in lines[::2])
    return [(head[1:], seq) for head, seq in zip(lines[::2], lines[1::2])]


@pytest.mark.parametrize('main', [cli.blocks_main, cli.convert_main,
                                  cli.dedupe_main, cli.tweak_main])
def test_help(main, capsys):
    with pytest.raises(SystemExit) as exc:
        main(['--help'])
    assert exc.value.code == 0
    assert capsys.readouterr().out.startswith('usage: cma-')


def test_missing_file_is_reported(tmp_path, capsys):
    missing = str(tmp_path / 'missing.cma')
    assert cli.blocks_main([missing]) == 1
    out, err = capsys.readouterr()
    assert out == ''
    assert err.startswith('cma-blocks: error: ')
    assert missing in err


def test_blocks_writes_aligned_columns(fikk_cma, capsys):
    assert cli.blocks_main([fikk_cma]) == 0
    records = read_fasta(capsys.readouterr().out)
    assert len(records) == 24
    assert records[0][0] == 'FIKK consensus'
    assert records[1][0].startswith(
        'psu|MAL7P1.127 | organism=Plasmodium_falciparum_3D7')
    assert {len(seq) for _head, seq in records} == {255}


def test_blocks_reads_standard_input(fikk_cma, monkeypatch, capsys):
    cli.blocks_main([fikk_cma])
    expected = capsys.readouterr().out
    with open(fikk_cma) as handle:
        monkeypatch.setattr('sys.stdin', handle)
        assert cli.blocks_main([]) == 0
    assert capsys.readouterr().out == expected


def test_blocks_excludes_flanks_across_blocks(aaaplus_sma, capsys):
    assert cli.blocks_main([aaaplus_sma]) == 0
    records = dict(read_fasta(capsys.readouterr().out))
    # 153 records in 13 blocks, less 4 whose bodies are wrapped over lines
    assert len(records) == 149
    assert {len(seq) for seq in records.values()} == {209}
    # The body of RfcB_human begins with an N-terminal flank, "{(HY)ELPW...";
    # the flank is not part of the alignment.
    rfcb = next(seq for head, seq in records.items()
                if head.startswith('RfcB_human'))
    assert rfcb.startswith('ELPWVEKYRP')


def test_convert_to_and_from_fasta(fikk_cma, tmp_path, capsys):
    ids = [seq['id'] for seq in cma.read(fikk_cma)['sequences']]
    assert cli.convert_main(['-o', 'fasta', fikk_cma]) == 0
    aligned = tmp_path / 'fikk.fasta'
    aligned.write_text(capsys.readouterr().out)
    assert [rec.id for rec in AlignIO.read(str(aligned), 'fasta')] == ids

    assert cli.convert_main(['-i', 'fasta', str(aligned)]) == 0
    converted = tmp_path / 'fikk.cma'
    converted.write_text(capsys.readouterr().out)
    assert [seq['id'] for seq in cma.read(str(converted))['sequences']] == ids


@pytest.mark.parametrize('option', ['-o', '-i'])
def test_convert_reports_unknown_format(option, fikk_cma, capsys):
    assert cli.convert_main([option, 'nosuchformat', fikk_cma]) == 1
    out, err = capsys.readouterr()
    assert out == ''
    assert err == "cma-convert: error: Unknown format 'nosuchformat'\n"


@pytest.mark.parametrize('options', [[], ['-o', 'fasta', '-i', 'fasta']])
def test_convert_requires_one_direction(options, fikk_cma):
    with pytest.raises(SystemExit) as exc:
        cli.convert_main(options + [fikk_cma])
    assert exc.value.code == 2


def test_dedupe_keeps_best_and_last_of_ties(fikk_cma, tmp_path, capsys):
    block = cma.read(fikk_cma)
    first, second, third, *rest = block['sequences']
    worse = dict(second, seq=second['seq'].lower(), description='worse')
    tie = dict(third, description='tie')
    block['sequences'] = [first, worse, second, third, tie] + rest
    infile = str(tmp_path / 'dupes.cma')
    cma.write(block, infile)

    assert cli.dedupe_main([infile]) == 0
    out, err = capsys.readouterr()
    assert err == "Dropped 2 duplicate sequences\n"
    outfile = tmp_path / 'deduped.cma'
    outfile.write_text(out)
    kept = cma.read(str(outfile))['sequences']
    assert [(seq['id'], seq['description']) for seq in kept[:3]] == [
        (first['id'], first['description']),
        (second['id'], second['description']),
        (third['id'], 'tie')]
    assert len(kept) == 24


def test_dedupe_handles_multiple_blocks(aaaplus_sma, capsys):
    assert cli.dedupe_main([aaaplus_sma]) == 0
    out, err = capsys.readouterr()
    assert err == "Dropped 0 duplicate sequences\n"
    assert out.count('_0].\n') == 13


def test_tweak_sets_level_and_name(fikk_cma, tmp_path):
    raised = str(tmp_path / 'raised.cma')
    assert cli.tweak_main([fikk_cma, '--level', '2', '-o', raised]) == 0
    assert cma.read(raised)['level'] == 2
    # Level 0 must also be settable, to reset a non-zero level.
    reset = str(tmp_path / 'reset.cma')
    assert cli.tweak_main(
        [raised, '--level', '0', '--name', 'Testfam', '-o', reset]) == 0
    block = cma.read(reset)
    assert (block['level'], block['name']) == (0, 'Testfam')


def test_tweak_adds_consensus_and_irons(fikk_cma, tmp_path):
    outfile = str(tmp_path / 'tweaked.cma')
    assert cli.tweak_main([fikk_cma, '-c', '--iron', '-o', outfile]) == 0
    sequences = cma.read(outfile)['sequences']
    assert len(sequences) == 25
    assert sequences[0]['id'] == 'consensus'
    assert all(cma.iron(seq['seq']) == seq['seq'] for seq in sequences)
