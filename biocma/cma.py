"""Lower-level functionality for parsing CMA (.cma) files.

This handles .chn (CHAIN) files as a collection of CMA blocks, i.e. alignments.
"""

import itertools
import logging
import re
from collections.abc import Generator
from typing import TypedDict

from . import sugar

logger = logging.getLogger(__name__)


class Sequence(TypedDict):
    """One sequence record of a CMA block, as cma.parse gives it."""

    index: int
    id: str
    description: str
    dbxrefs: dict[str, str]
    phylum: str
    taxchar: str
    head_len: int
    tail_len: int
    head_seq: str
    tail_seq: str
    length: int
    seq: str


class Block(TypedDict):
    """One CMA alignment block, as cma.parse gives it."""

    level: int
    one: int
    name: str
    params: str
    query_length: int
    query_chars: str
    sequences: list[Sequence]


def parse(infile) -> Generator[Block, None, None]:
    """Parse each block of a CMA file, given as a path or an open handle.

    Raise ValueError, naming the block, if the input ends inside a block or a
    block's structure is malformed. Unparseable sequence lines are skipped
    with a warning.
    """
    with sugar.maybe_open(infile) as instream:
        yield from _parse_blocks(instream)


read = sugar.make_reader(parse)

# --------------------------------------------------------------------
# Parse the cma into raw chunks

# Flow:
# if START:
#   parse line 1
#   parse line 2
#   pass control to _parse_sequences
#
#       _parse_sequences:
#       read a non-blank line;
#           if END: return (w/ all sequences accumulated)
#           else:
#               read 2 more non-blank lines
#               parse the 3 read lines & bundle up data
#               yield a seq
#


def _parse_blocks(instream) -> Generator[Block, None, None]:
    """Parse an alignment block from the given file handle.

    Block looks like:
    [0_(1)=fa2cma(8){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:
    (209)***********************************************...

    ... sequences, numbered 1-8 ...

    _0].

    """
    ilines = sugar.unblank(instream)
    for line in ilines:
        if line.startswith("["):
            # Start of block
            try:
                level, one, name, seqcount, params = _parse_block_header(line)
                qlen, qchars = _parse_block_postheader(
                    _next_line(ilines, "the query line")
                )
                # Pass control to the sequence parser
                sequences = list(_parse_sequences(ilines, qlen))
            except UnicodeDecodeError:
                # An encoding error, not a malformed block
                raise
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Malformed CMA block {line!r}: {exc}") from exc
            # Validation
            if not len(sequences) == seqcount:
                logger.warning(
                    "Expected %d sequences in block %s, found %d",
                    seqcount,
                    name,
                    len(sequences),
                )
            yield {
                "level": level,
                "one": one,
                "name": name,
                "params": params,
                "query_length": qlen,
                "query_chars": qchars,
                "sequences": sequences,
            }


def _parse_sequences(ilines, expect_qlen) -> Generator[Sequence, None, None]:
    """Parse the sequences in the current block.

    Sequence looks like:

    $3=227(209):
    >gi|15606894|ref|NP_214275.1|  {|2(244)|<Aquificae(B)>}DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|2984127|gb|AAC07663.1| DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|75
    {()YVPFARKYRPKFFREVIGQEAPVRILKNAIKNDRVAHaYLFAGPRGVGKTTIARILAKALNcknpskgepcgecencreiDRGVFPDLIEMDAASNRGIDDVRA-LKEAVNYKPIKG-KYKVYIIDEAHMLTKEAFNALLKTLEEPPPRTVFVLCTTEYDKILPTILSRCQRIIFSKVRKEKVIEYLKKICEKEGIECEEGALEVLAHASEGCMRDAASLLDQASVYGE()}*

    """  # noqa: E501
    while True:
        first = _next_line(ilines, "a sequence or the end of the block")
        if first.startswith("_") and first.endswith("]."):
            # End of sequences & end of block
            break

        # ENH: handle wrapped lines?
        try:
            index, this_len, query_len = _parse_seq_preheader(first)
        except ValueError:
            logger.warning("Unparseable line (SKIPPING):\n%s", first)
            continue
        (rec_id, dbxrefs, headlen, taillen, phylum, taxchar, description) = (
            _parse_seq_header(_next_line(ilines, f"the header line after {first}"))
        )
        body = _next_line(ilines, f"the sequence line of {rec_id}")
        try:
            headseq, molseq, tailseq = _parse_seq_body(body)
        except (IndexError, ValueError):
            logger.warning("Unparseable sequence: %s -- SKIPPING", rec_id)
            continue

        # Validation
        if expect_qlen != query_len:
            logger.warning(
                "Query length in %s given as %d; expected %d",
                rec_id,
                query_len,
                expect_qlen,
            )
        if not headseq and not headlen:
            headlen = 0
        if not tailseq and not taillen:
            taillen = 0
        if headseq:
            if headlen is None:
                headlen = len(headseq)
            elif headlen != len(headseq):
                logger.warning(
                    "Conflicting head flank lengths in %s: %d, %d",
                    rec_id,
                    headlen,
                    len(headseq),
                )
        if tailseq:
            if taillen is None:
                taillen = len(tailseq)
            elif taillen != len(tailseq):
                logger.warning(
                    "Conflicting tail flank lengths in %s: %d, %d",
                    rec_id,
                    taillen,
                    len(tailseq),
                )

        yield {
            "index": index,
            "id": rec_id,
            "description": description,
            "dbxrefs": dbxrefs,
            "phylum": phylum,
            "taxchar": taxchar,
            "head_len": headlen,
            "tail_len": taillen,
            "head_seq": headseq,
            "tail_seq": tailseq,
            "length": this_len,
            "seq": molseq,
        }


def _next_line(ilines, expected):
    """Return the next line of the current block.

    Raise ValueError naming what was expected if the input ends, or if the
    line is the header of another block.
    """
    try:
        line = next(ilines)
    except StopIteration:
        raise ValueError(f"the input ends where {expected} was expected") from None
    if line.startswith("["):
        raise ValueError(f"the block does not end before the next block {line!r}")
    return line


# Microparsing


def _parse_block_header(line) -> tuple[int, int, str, int, str]:
    """
    [0_(1)=fa2cma(8){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:

    The level is a signed integer of any length. The name ends at the first
    '(', as in MAPGAPS, which reads it with the scanf conversion %[^(].
    """
    level, _rest = line[1:].split("_", 1)
    one, _rest = _rest[1:].split(")=", 1)
    name, _rest = _rest.split("(", 1)
    seqcount, _rest = _rest.split(")", 1)
    params = _rest.strip("{}:")
    return int(level), int(one), name, int(seqcount), params


def _parse_block_postheader(line) -> tuple[int, str]:
    """
    (209)**************!*****************!!*************...
    """
    parts = line[1:].split(")", 1)
    qlen = int(parts[0])
    if not len(parts[1]) == qlen:
        logger.warning(
            "postheader expected %d-long query, found %d", qlen, len(parts[1])
        )
    return qlen, parts[1]


def _parse_seq_preheader(line) -> tuple[int, int, int]:
    """
    $3=227(209):
    """
    match = re.match(r"\$ (\d+) = (\d+) \( (\d+) \):", line, re.VERBOSE)
    if not match:
        raise ValueError("Unparseable header: " + line)
    index, this_len, query_len = match.groups()
    return int(index), int(this_len), int(query_len)


def _parse_seq_header(line):
    """Unique ID, head/tail lengths and taxonomy info from a sequence header.

    The description is the part of the FASTA/CMA sequence header starting after
    the first space (i.e. excluding ID), to the end of the line.

    This function looks inside the first '{...}' pair to extract info.

    Ex:

    >consensus seq
    >gi|15606894|ref|NP_214275.1|  {|2(244)|<Aquificae(B)>}DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|2984127|gb|AAC07663.1| DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|75
    >gi|3212262|pdb|1A2K|C  {<Chordata(M)>}Chain C, Gdpran-Ntf2 Complex >gi|3212263|pdb|1A2K|D Chain D, Gdpran-Ntf2 Complex >gi|3212264|pdb|1A2K|E Chain E, Gdpran-Ntf2 Complex >gi|5542273|pdb|1IBR|A C

    """  # noqa: E501
    _parts = line[1:].split(None, 1)
    rec_id = _parts[0]
    descr = _parts[1] if _parts[1:] else ""
    # Database cross references
    dbxrefs = {}
    if "|" in rec_id:
        id_gen = iter(rec_id.rstrip("|").split("|"))
        for key in id_gen:
            try:
                dbxrefs[key] = next(id_gen)
            except StopIteration:
                break

    # Head/tail lengths and taxonomy codes
    headlen = taillen = None
    phylum = taxchar = ""
    if descr.startswith("{"):
        _deets, description = descr[1:].split("}", 1)
        match = re.search(
            r"""
        (?:
        \| (?P<headlen> \d+)
            \( (?P<taillen> \d+)
            \)
        \|
        )?
        (?:
        <   (?P<phylum> .+?)
            \( (?P<taxchar> \w)
            \)
        >
        )?
        """,
            _deets,
            re.VERBOSE,
        )
        if match:
            headlen, taillen, phylum, taxchar = match.groups()
            if headlen is not None:
                headlen = int(headlen)
            if taillen is not None:
                taillen = int(taillen)
            if phylum is None:
                phylum = ""
            if taxchar is None:
                taxchar = ""
        else:
            logger.warning("Couldn't match head/tail: %s", _deets)
    else:
        description = descr

    # TODO - return a dictionary here, update it in _parse_sequences
    return rec_id, dbxrefs, headlen, taillen, phylum, taxchar, description


def _parse_seq_body(line):
    """
    Ex:

    {()YVPFARKYRPKFFREVIGQEAPVRILKALNcknpskgepcgereiDRGVFPDVRA-LKLLDQASVYGE()}*
    MENINNI{()----------FKLILVGDGKFFSSSGEIIFNIWDTKFGGLRDGYYRLTYKNEDDM()}*

    Or:

    {(HY)ELPWVEKYR...

    The sequence fragments in parentheses represent N- or C-terminal flanking
    regions that are not part of the alignment block (I think). Most tools don't
    include them, but some do, apparently.
    """
    line = line.rstrip("*")
    if "{()" in line:
        head, _rest = line.split("{()", 1)
    else:
        # Match parens
        _rest = line.split("{(", 1)[1]
        head, _rest = _rest.split(")", 1)
    if "()}" in _rest:
        molseq, tail = _rest.split("()}", 1)
    else:
        # Match parens
        molseq, _rest = _rest.split("(", 1)
        tail = _rest.split(")}", 1)[0]
    return (head, molseq, tail)


# --------------------------------------------------------------------
# Write


def write(blocks, outfile):
    if isinstance(blocks, dict) and "sequences" in blocks:
        blocks = [blocks]
    with sugar.maybe_open(outfile, "w") as outstream:
        outstream.writelines(itertools.chain(*map(_format_block, blocks)))


def _format_block(block):
    # Block header
    # [0_(1)=structs.seq(12){go=19,gx=2,pn=5.0,lf=0,rf=0}:
    # (69)*****!***************************!!!******!******************!*******
    yield """\
[{level}_({one})={name}({seqcount}){{{params}}}:
({query_length}){query_chars}

""".format(seqcount=len(block["sequences"]), **block)
    # Sequences
    yield from _format_sequences(block["sequences"], block["query_length"])
    # Block close
    yield f"_{block['level']}].\n"


def _format_sequences(sequences, query_length):
    """Format each sequence record as CMA text.

    A formatted sequence looks like:

    $1=254(255):
    >2QG5|A  {|27(6)|}
    {()YTLENTIGRGSWGEVKIAVQKGTRIRRAAKKIPKYFV---EDVDRFKQEIEIMKSLDHPNIIRLYETFEDNTDIYLVMELCTGGELFERVVHKRVFRESDAARIMKDVLSAVAYCHKLNVAHRDLKPENFLFltdSPDSPLKLIDFGLAARFkpGKMMRTKVGTPYYVSPQVLEGL-YGPECDEWSAGVMMYVLLCGYPPFSAPTDXEVMLKIREGTFtfpeKDWLNVSPQAESLIRRLLTKSPKQRIT-----SLQALEHEW-()}*

    """
    for idx, seq in enumerate(sequences):
        head_tail = (
            "|{head_len}({tail_len})|".format(**seq)
            if seq["head_len"] or seq["tail_len"]
            else ""
        )
        taxonomy = (
            "<{phylum}({taxchar})>".format(**seq)
            if seq["phylum"] and seq["taxchar"]
            else ""
        )
        special_header = (
            "{" + head_tail + taxonomy + "}" if head_tail or taxonomy else ""
        )
        yield """\
${index}={length}({0}):
>{id} {1}{description}
{head_seq}{{(){seq}()}}{tail_seq}*

""".format(query_length, special_header, **dict(seq, index=idx + 1))


# --------------------------------------------------------------------
# Utilities


def realign_seqs(block, gap_char=".", align_indels=False):
    """Add gaps to a block so all residues in a column are equivalent.

    Given a block, containing a list of "sequences" (dicts) each containing a
    "seq" (actual string sequence, where upper=match, lower=insert, dash=gap),
    insert gaps (- or .) into the sequences s.t.

    1. columns line up properly, and
    2. all resulting sequences have the same length


    The reason this needs to be done is that the query/consensus sequence is not
    assigned gaps to account for inserts in the other sequences. We need to add
    the gaps back to obtain a normal alignment.

    `return`: a list of realigned sequence strings, empty for a block without
    sequences.
    """
    # ENH: align inserts using an external tool (if align_indels)
    all_chars = [list(sq["seq"]) for sq in block["sequences"]]
    if not all_chars:
        return []
    # NB: If speed is an issue here, consider Numpy or Cython
    #   main problem: list.insert is O(n) -- would OrderedDict help?
    nrows = len(all_chars)
    i = 0
    while i < len(all_chars[0]):
        rows_need_gaps = [r for r in all_chars if not r[i].islower()]
        if len(rows_need_gaps) != nrows:
            for row in rows_need_gaps:
                row.insert(i, gap_char)
        i += 1
    return ["".join(row) for row in all_chars]


def consensus2block(
    record, level=0, name=None, molseq=None, params=None, query_chars=None
) -> Block:
    """Convert a Biopython SeqRecord to a cma block.

    Ungapping is handled here.

    `molseq`: if given, a string used in place of the record's own sequence.
    The record is not modified.

    `params`, `query_chars`: block values to keep, as a ChainMultiAlignment
    holds them. Without `params` the block gets the parameters of MAPGAPS's
    fa2cma; `query_chars` is kept only if its length equals the number of
    residues in the consensus, and is otherwise all '*'.
    """
    if molseq is None:
        molseq = str(record.seq)
    cons_ungap = molseq.replace("-", "").replace(".", "").upper()
    qlen = len(cons_ungap)
    if not query_chars or len(query_chars) != qlen:
        query_chars = "*" * qlen
    return dict(
        level=level,
        one=1,
        name=name or record.id,
        params=params or "go=10000,gx=2000,pn=1000.0,lf=0,rf=0",
        query_length=qlen,
        query_chars=query_chars,
        sequences=[seqrecord2sequence(record, qlen, 1, molseq=cons_ungap)],
    )


def seqrecord2sequence(record, qlen, index, molseq=None, head_extra="") -> Sequence:
    """Convert a Biopython SeqRecord to a cma sequence record.

    Indels (gaps, casing) must have already been handled in the record or, if
    given, in `molseq`, a string used in place of the record's own sequence.

    A record read by Bio.SeqIO from FASTA has a description that begins with
    its ID, while a record read from CMA has one without it. As in Biopython's
    FASTA writer, the first word of the description is removed only when it
    equals the record's ID, so a CMA description that happens to begin with
    its own ID loses that word. An empty remainder, as from a Clustal record
    whose description is its ID, is written as a single space, which keeps
    the header bytes of earlier releases.

    The flank lengths, flank sequences and taxonomy that
    ChainMultiAlignment.from_block stores in the record's annotations are
    kept; a record without them, or with None for them, has empty flanks and
    no taxonomy. The `length` field counts the residues in the flank
    sequences as well as in `molseq`.

    `head_extra`: residues of the record that precede its aligned region but
    are not in `molseq`. They are counted in the head flank length, and are
    appended to the head flank sequence if the record has one.
    """
    if molseq is None:
        molseq = str(record.seq)
    words = record.description.split(None, 1)
    if words and words[0] == record.id:
        description = words[1] if len(words) > 1 else ""
    else:
        description = record.description
    annotations = record.annotations
    head_seq = annotations.get("head_seq") or ""
    if head_seq:
        head_seq += head_extra
    tail_seq = annotations.get("tail_seq") or ""
    return dict(
        index=index,
        id=record.id,
        description=description or " ",
        dbxrefs={},
        phylum=annotations.get("phylum") or "",
        taxchar=annotations.get("taxchar") or "",
        head_len=(annotations.get("head_len") or 0) + len(head_extra),
        tail_len=annotations.get("tail_len") or 0,
        head_seq=head_seq,
        tail_seq=tail_seq,
        length=len(head_seq) + len(molseq) - molseq.count("-") + len(tail_seq),
        seq=molseq,
    )


def replace_asterisks(seq, label=None):
    if "*" in seq:
        logger.warning(
            "Sequence %scontains '*' character; replacing with 'X'",
            str(label) + " " if label else "",
        )
    return str(seq).replace("*", "X")


def collapse_to_consensus(seqrecords, strict=False, do_iron=True) -> Block:
    """Opposite of realign_seqs.

    Input sequences should all be the same length.

    The first record must be the consensus.

    The block takes its level, name, params and query_chars from attributes
    of `seqrecords` where it has them, as a ChainMultiAlignment does. The
    residues of a record that precede the first consensus column are removed
    from its sequence and counted in its head flank.
    """
    level = getattr(seqrecords, "level", 0)
    name = getattr(seqrecords, "name", None)
    params = getattr(seqrecords, "params", None)
    query_chars = getattr(seqrecords, "query_chars", None)
    # rows is a new list; the caller's input is not modified
    consensus, *rows = seqrecords
    cons_seq = str(consensus.seq)
    cons_length = len(cons_seq)
    for index, s in enumerate(rows, 2):
        if len(s) != cons_length:
            raise ValueError(
                f"Sequence #{index} has length {len(s)}, consensus is {cons_length}"
            )

    if "." in cons_seq:
        # Strict -- error if there's a '-'
        if "-" in cons_seq:
            if strict:
                raise ValueError("Consensus contains '-' gap characters")
            logger.warning(
                "Consensus sequence contains both '.' and '-' gap "
                "characters -- is it really the consensus?"
            )
            aligned_cols = [(c not in ".-") for c in cons_seq]
        else:
            aligned_cols = [c != "." for c in cons_seq]
    else:
        # A little more ambiguous...
        aligned_cols = [c != "-" for c in cons_seq]

    # Start a block with the consensus sequence
    block = consensus2block(
        consensus,
        level=level,
        name=name,
        molseq=replace_asterisks(cons_seq, "consensus"),
        params=params,
        query_chars=query_chars,
    )
    qlen = block["query_length"]

    # Collapse & add remaining sequences to the block
    for index, rec in enumerate(rows, 2):
        # Collapse rec.seq down to aligned size
        new_mol_seq = []
        head_extra = []
        is_beginning = True
        for aligned_col, char in zip(
            aligned_cols, replace_asterisks(rec.seq, index), strict=True
        ):
            if aligned_col:
                is_beginning = False
                if char in "-.":
                    # deletion
                    new_mol_seq.append("-")
                else:
                    # aligned character
                    new_mol_seq.append(char.upper())
            elif char not in "-.":
                # An insert, or a residue before the first consensus column,
                # which belongs to the head flank
                if is_beginning:
                    head_extra.append(char.upper())
                else:
                    new_mol_seq.append(char.lower())
        molseq = "".join(new_mol_seq)
        if do_iron:
            molseq = iron(molseq)
        block["sequences"].append(
            seqrecord2sequence(
                rec, qlen, index, molseq=molseq, head_extra="".join(head_extra)
            )
        )

    return block


def iron(sequence):
    """'Iron out' indel regions in the aligned sequence.

    Any inserts next to deletions are converted to matches (uppercase).

    Given a CMA string like:
        AAAAbc--de-f--gAAA
    Result:
        AAAABCDEFgAAA
    """
    r_indel = re.compile(r"(-[a-y]|[a-y]-)")
    orig_sequence = sequence
    while r_indel.search(sequence):
        in_insert = False
        in_gap = False
        seen_gaps = 0
        inserts = []
        outchars = []

        for char in sequence:
            if in_insert:
                if char.islower():
                    # Extend the insert
                    inserts.append(char)
                elif char.isupper():
                    # Indel is over; 'iron' out & emit inserts, then gaps
                    in_insert = False
                    outchars.extend(inserts)
                    inserts = []
                    outchars.append("-" * seen_gaps)
                    seen_gaps = 0
                    outchars.append(char)
                else:
                    # Convert a preceding indel char to a 'match' (uppercase)
                    # If the indel and gap are both multiple chars, this will
                    # capitalize the insert left-to-right, then leave any gap
                    # remainer as-is.
                    assert char == "-"
                    if not inserts:
                        in_insert = False
                        in_gap = True
                        seen_gaps += 1
                    else:
                        outchars.append(inserts.pop(0).upper())
                        # NB: Only leave the insert region if we've finished
                        # converting all the insert chars
                        if not inserts:
                            in_insert = False
                            in_gap = True

            elif in_gap:
                if char.islower():
                    in_insert = True
                    in_gap = False
                    # If some inserts previously seen, emit them now
                    # If no inserts have been seen yet, we'll iron this indel
                    if inserts:
                        outchars.extend(inserts)
                        outchars.append("-" * seen_gaps)
                        seen_gaps = 0
                    inserts = [char]
                elif char.isupper():
                    in_gap = False
                    # End of the gap -- emit
                    if inserts:
                        outchars.extend(inserts)
                        inserts = []
                    outchars.append("-" * seen_gaps)
                    seen_gaps = 0
                    outchars.append(char)
                else:
                    # Extend the gap
                    assert char == "-"
                    seen_gaps += 1

            else:
                assert not inserts and not seen_gaps, (
                    f"Inserts: {inserts}, gaps: {seen_gaps}, seq: {sequence}, "
                    f"in_ins={in_insert}, in_gap={in_gap}"
                )
                # Coming from Match state
                if char.isupper():
                    # Extend the match
                    outchars.append(char)
                elif char.islower():
                    inserts.append(char)
                    in_insert = True
                else:
                    assert char == "-"
                    seen_gaps += 1
                    in_gap = True

        # Emit any trailing indel
        if inserts:
            outchars.extend(inserts)
        if seen_gaps:
            outchars.append("-" * seen_gaps)
        sequence = "".join(outchars)
        assert (
            sequence.replace("-", "").upper() == orig_sequence.replace("-", "").upper()
        ), "\nOrig: " + orig_sequence + "\nIron: " + sequence
    return sequence
