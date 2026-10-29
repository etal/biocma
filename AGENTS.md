# BioCMA

Two filenames reach this one file: OMP and Codex read it as `AGENTS.md`, while Claude Code reads it as `CLAUDE.md`, a symlink recorded in git.
Where symlinks do not survive checkout (Git for Windows outside Developer Mode, or a DrvFS mount such as `/mnt/c`), `CLAUDE.md` arrives as a stub naming this file; read `AGENTS.md` directly rather than pasting content into the stub.

Reading and writing the Consensus Multiple Alignment (CMA) format of Andrew Neuwald's MAPGAPS and CHAIN programs, with Biopython alignments, and four `cma-*` command-line programs.
The 2013 Python 2 code base has been ported to Python 3.12+ and current Biopython; `README.rst` is the user and developer documentation.

## Consumers

[Fammer](https://github.com/etal/fammer) pins biocma to a git commit, and [BioFrills](https://github.com/etal/biofrills) imports `biocma.cma` when it is installed, without declaring it as a dependency.
The names that either one calls are the compatibility contract.

- Fammer writes its MAPGAPS profiles with `biocma.biocma.ChainMultiAlignment(aln, level=level)` followed by `biocma.biocma.write([cmaln], path, do_iron=True)`, where the first record of `aln` is the consensus.
  Treat any change to that call sequence as a change of output, and check it against Fammer's reference files with a checkout of Fammer beside this repository.
  For each case in `../fammer/tests/data/py27/library/mapgaps.json`, read the Clustal alignment that `alignment` names, relative to `../fammer/tests/data/py27/`.
  Prepend `consensus` as a record whose ID is the stem of `target` and whose description is that ID followed by ` consensus`, and write the result with `level`.
  Compare the output byte for byte with the file that `target` names under `../fammer/tests/data/py27/library/mapgaps/`.
- BioFrills' `seq-diff` and `seq-purge` call `cma.read` and `cma.write`, and read the `id`, `description` and `seq` keys of each sequence record.
  After a change to any of these, run `PYTHONPATH=$PWD/../biocma uv run pytest tests/test_cli.py` in a checkout of BioFrills beside this repository; its CMA cases are skipped when biocma cannot be imported.

## Layout

- `biocma/cma.py`: the CMA parser and writer, and the block operations `realign_seqs`, `collapse_to_consensus` and `iron`.
  Blocks and sequence records are plain dicts described by the `TypedDict` classes `Block` and `Sequence`.
- `biocma/biocma.py`: `ChainMultiAlignment`, a subclass of `MultipleSeqAlignment`, and `parse`, `read`, `write` and `convert` between CMA and the `Bio.AlignIO` formats.
- `biocma/utils.py`: consensus, conservation, equivalent-position and insert functions over blocks.
- `biocma/cli.py`: one `*_main` function per console script registered in `[project.scripts]`; the `_console` decorator gives every command the same handling of input errors and broken pipes.
- `biocma/sugar.py`: `maybe_open`, which accepts a path or a handle, and `make_reader`, which builds `cma.read` and `biocma.read` from the matching `parse`.
- `test/`: the pytest suite and its input files `fikk.cma`, `fikk-full.fasta` and `AAAPlus.sma` (several blocks), which `test/conftest.py` provides as fixtures.

`cma.py`, `utils.py` and `sugar.py` import no Biopython, and `cli.py` imports `biocma.biocma` only inside `convert_main`, so the other commands start without loading Biopython.

## Commands

```sh
pixi run test              # newest Python in pixi.lock
pixi run -e py312 test     # Python 3.12, the oldest supported version
pixi run lint              # ruff check
pixi run fmt               # ruff format
pixi run typecheck         # mypy
uv run pytest              # the same suite with the PyPI dependencies of uv.lock
```

The tests call the `*_main` functions directly, so after a change to `[project.scripts]` run `pixi install` or `uv sync` and check that the installed scripts resolve, as CI does with `--help`.
CI also runs the suite at the declared dependency floors with `UV_RESOLUTION=lowest-direct` and `UV_PYTHON=3.12`; set these variables on the command, never with `export` in a persistent shell, where they turn the next `uv run --locked` into a lock failure.

## The CMA format

- No written specification exists; the authority is what MAPGAPS reads.
  A block header is `[LEVEL_(ONE)=NAME(NSEQ){go=…,gx=…,pn=…,lf=…,rf=…}:` and a block ends with `_LEVEL].`; the level is a signed integer of any length, and the name runs to the first `(`.
- In a sequence record's `$N=LENGTH(QLEN):` line, `LENGTH` counts the residues of the head flank text, the aligned row without its `-` characters, and the tail flank text.
  The head and tail lengths in the `{|HEAD(TAIL)|}` header are the authority, and the head length is the offset from which residues are numbered.
  MAPGAPS files often give a flank length with no flank text, or flank text of another length; the parser's logged warning "Conflicting head flank lengths" reflects an assumption of biocma, not a rule of the format.
- In an aligned row, uppercase letters are match columns, lowercase letters are inserts and `-` is a deletion; `realign_seqs` fills insert columns with `.`.
  Biopython reads only `-` as a gap, so `ChainMultiAlignment.alignment` and `.substitutions` replace `.` with `-` first.

## Testing

- pytest turns messages of the `warnings` module into errors; the parser's validation messages are `logging` records, which tests read with `caplog`.
- A round trip cannot detect a parser that drops a field, because the writer reproduces the loss and the second parse agrees with the first.
  `test_write_reproduces_source_text` therefore compares the whitespace-separated tokens of the written output with those of the fixture text, and a field that no fixture sets to a value other than its default, such as the block level, needs a small hand-built input.
- A shell that exports `PYTHONUNBUFFERED=1`, as some agent harnesses do, hides the broken-pipe failure at interpreter exit that `cli._console` guards against.
  Measure that behavior with `env -u PYTHONUNBUFFERED`, and keep `test_closed_pipe_in_a_real_process` removing the variable from its child process.
- The pre-port source is commit `eac0c57`.
  The 2020 pull request that first ported the code, commit `8f97c3b`, was merged as `33a600f` with the tree of its first parent `8fbf1ad`, so the merge discarded that pull request's changes, and the current port builds on `8fbf1ad`.
  To compare behavior with `eac0c57` on Apple Silicon, extract it with `git archive eac0c57 biocma/` into a scratch directory and run it with `pixi exec -p osx-64 -s python=2.7 -s biopython python …`, with that directory first on `sys.path`.
  `utils.get_consensus` breaks ties by the first residue seen, while Python 2 broke them in string-hash order, so consensus strings can differ from the pre-port output at tied columns.

## Pitfalls

- `Bio.AlignIO` opens an output path in text mode, and so truncates it before it raises `StreamModeError` for a binary format such as seqxml.
  `biocma.biocma.convert` decides the mode of an output format with `_writes_bytes` before it opens the output; keep that order.
- mypy checks `biocma/` against `python_version = "3.12"`, set explicitly rather than taken from the interpreter that runs it.
- Library modules log through `logging.getLogger(__name__)`.
- Runtime dependency floors default to the versions that the current Ubuntu LTS (26.04) ships; raise one only for a concrete problem, such as a security advisory or a conflict.
  The Python floor, 3.12, is a separate decision.
  Biopython's floor is 1.87, one version above the LTS: Biopython 1.86 reads UniProt XML from a text handle with a deprecation warning instead of raising `StreamModeError`, so `convert` never repeats the read in binary mode and reads such input as text.
  After a change to a dependency, update both `uv.lock` and `pixi.lock`.

## Conventions

- Format with `ruff format`, and list each formatting-only commit in `.git-blame-ignore-revs`.
- Commit directly to `master`, without a feature branch or pull request.
- Commit messages: an imperative subject line, then a short paragraph on what changed and why; American English.
- In Markdown files, put each sentence on its own line.
