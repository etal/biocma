from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent


@pytest.fixture
def fikk_cma():
    # str, not Path: sugar.maybe_open treats any non-str argument as a handle
    return str(DATA_DIR / "fikk.cma")


@pytest.fixture
def fikk_fasta():
    return str(DATA_DIR / "fikk-full.fasta")


@pytest.fixture
def aaaplus_sma():
    return str(DATA_DIR / "AAAPlus.sma")
