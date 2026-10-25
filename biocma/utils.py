from collections import Counter

from .cma import Block, Sequence


def find_seq_rec(block, name, case_sensitive=True):
    """Given part of a sequence ID, find the first matching record."""
    if case_sensitive:

        def test(name, rec):
            return name in rec["id"]
    else:

        def test(name, rec):
            return name.upper() in rec["id"].upper()

    for rec in block["sequences"]:
        if test(name, rec):
            return rec
    raise ValueError(f"No sequence ID matches {name!r}")


def find_seq_id(block, name, case_sensitive=True):
    """Given part of a sequence ID, find the first actual ID that contains it.

    Example::

        >>> find_seq_id(block, '2QG5')
        'gi|158430190|pdb|2QG5|A'

    Raise a ValueError if no matching key is found.
    """
    rec = find_seq_rec(block, name, case_sensitive)
    return rec["id"]


def get_consensus(block: Block) -> Sequence:
    """Calculate a simple consensus sequence for the block.

    Raise a ValueError if the block has no sequences, or if the number of
    aligned (non-insert) characters in any sequence differs from the block's
    query length.
    """
    if not block["sequences"]:
        raise ValueError("Block has no sequences")
    # Take aligned (non-insert) chars from all rows; transpose
    rows = []
    for row in block["sequences"]:
        aligned = [c for c in row["seq"] if not c.islower()]
        if len(aligned) != block["query_length"]:
            raise ValueError(
                f"Sequence {row['id']} has {len(aligned)} aligned characters, "
                f"query length is {block['query_length']}"
            )
        rows.append(aligned)
    columns = zip(*rows, strict=True)
    cons_chars = [Counter(col).most_common()[0][0] for col in columns]
    cons_chars = [c if c != "-" else "X" for c in cons_chars]
    return {
        "index": 1,
        "id": "consensus",
        "description": "",
        "dbxrefs": {},
        "phylum": "",
        "taxchar": "",
        "head_len": 0,
        "tail_len": 0,
        "head_seq": "",
        "tail_seq": "",
        "length": block["query_length"],
        "seq": "".join(cons_chars),
    }


def get_conservation(block):
    """Calculate conservation levels at each consensus position.

    Return a dict of {position: float conservation}
    """
    consensus = block["sequences"][0]["seq"]
    assert all(c.isupper() for c in consensus), "So-called consensus contains indels!"
    # remove all non-consensus positions -- now alignment is easy
    cleaned = [[c for c in s["seq"] if not c.islower()] for s in block["sequences"][1:]]
    height = float(len(cleaned))
    # validation
    for row in cleaned:
        if len(row) != len(consensus):
            raise ValueError(
                f"Aligned sequence length ({len(row)}) doesn't match "
                f"consensus ({len(consensus)})"
            )
    # transpose & go
    columns = list(zip(*cleaned, strict=True))
    return dict(
        (idx + 1, columns[idx].count(cons_char) / height)
        for idx, cons_char in enumerate(consensus)
    )


def get_equivalent_positions(block):
    """Create a mapping of equivalent residue positions to consensus.

    Build a dict-of-dicts::

        {consensus-posn: {id: equiv-posn, id: equiv-posn, ...}, ...}

    The first sequence in the alignment is assumed to be the (gapless) consensus
    sequence.
    """
    consensus = block["sequences"][0]["seq"]
    rest = block["sequences"][1:]

    # Validation
    if "-" in consensus or "." in consensus:
        raise ValueError("First sequence (consensus?) contains gaps")
    # Check for duplicate sequence IDs
    seen = set()
    dupes = set()
    for rec in rest:
        if rec["id"] in seen:
            dupes.add(rec["id"])
        else:
            seen.add(rec["id"])
    if dupes:
        raise ValueError("Duplicate sequences:\n" + "\n".join(dupes))

    curr_shift = {}
    curr_resn = {}
    # NB: consensus doesn't have head/tail, but other sequences may
    for rec in rest:
        # Count inserts seen so far -- shift string indexes by this far ahead to
        # get the "equivalent" location in the sequence string
        #   - as in, how far ahead in the current seq do we need to jump to get
        #     to a position equivalent to what's in the consensus?
        #   - can this ever be less than 0 (==consensus)?  No, because that's
        #     where gaps come from. Good.
        curr_shift[rec["id"]] = 0
        # Residue number in the actual sequence at the current (shifted)
        # location
        #   curr_posn[id] = current equivalent res.num in `id` to cons[i]
        curr_resn[rec["id"]] = rec["head_len"]

    equivalencies: dict[int, dict[str, int]] = dict(
        (i + 1, {}) for i in range(len(consensus))
    )
    # Map each character position i in the consensus sequence
    # to equivalent residues in each of the other sequences
    #   i = index in the consensus string (== consensus res.num - 1)
    for i, char in enumerate(consensus):
        assert char.isupper()
        for rec in rest:
            rid = rec["id"]
            strposn = i + curr_shift[rid]
            if rec["seq"][strposn].isupper():
                # Match
                curr_resn[rid] += 1
            elif rec["seq"][strposn].islower():
                # Insert
                while rec["seq"][strposn].islower():
                    # Count the whole insert size
                    curr_shift[rid] += 1
                    curr_resn[rid] += 1
                    strposn += 1
                curr_resn[rid] += 1  # Count the next match, too
            else:
                # Deletion / gap
                assert rec["seq"][strposn] in ".-"
                continue
            equivalencies[i + 1][rid] = curr_resn[rid]

    return equivalencies


def get_inserts(block):
    """Identify the inserts in sequence in a block.

    Inserts are relative to the consensus (theoretically), and identified by
    lowercase letters in the sequence. The returned integer pairs represent the
    insert start and end positions in the full-length sequence, using one-based
    numbering.

    The first sequence of the CMA block is included, though it may just be the
    consensus sequence, which shouldn't have any inserts.

    Output:
        {id1: [(start, end), (start, end), ...], id2: ..., ...}

    """

    def find_inserts(seq, head_len):
        """Locate the lowercase regions in a character sequence.

        Yield the insert ranges as tuples using 1-based numbering, shifted by
        head_len.
        """
        in_insert = False
        curr_start = None
        deletions = 0
        for idx, is_lower in enumerate(map(str.islower, seq)):
            if is_lower:
                if not in_insert:
                    # Start of a new insert region
                    curr_start = head_len + idx + 1 - deletions
                    in_insert = True
            else:
                if in_insert:
                    # End of the current insert region
                    yield (curr_start, head_len + idx - deletions)
                    in_insert = False
                if seq[idx] == "-":
                    deletions += 1

    return dict(
        (record["id"], list(find_inserts(record["seq"], record["head_len"])))
        for record in block["sequences"]
    )


def number_letters(record):
    """Return a dict of {posn: restype} for each letter in the sequence."""
    # Ungapped sequence -- accounts for dels vs. consensus
    seq = record["seq"].replace("-", "").replace(".", "")
    # Difference between string positions and residue numbers
    shift = record["head_len"] + 1
    return dict((idx + shift, letter) for idx, letter in enumerate(seq))
