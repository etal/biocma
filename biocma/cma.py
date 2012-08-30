#!/usr/bin/env python

"""Lower-level functionality for parsing CMA (.cma) files.

This handles .chn (CHAIN) files as a collection of CMA blocks, i.e. alignments.
"""

import itertools
import logging
import re

from . import sugar


def parse(infile):
    with sugar.maybe_open(infile) as instream:
        for block in _parse_blocks(instream):
            yield block

read = sugar.make_reader(parse)

# Removed
# parse_raw = parse

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

def _parse_blocks(instream):
    """Parse an alignment block from the given file handle.

    Block looks like:
    [0_(1)=fa2cma(8){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:
    (209)***********************************************...

    ... sequences, numbered 1-8 ...

    _0].

    """
    ilines = sugar.unblank(instream)
    for line in ilines:
        if line.startswith('['):
            # Start of block
            level, one, name, seqcount, params = _parse_block_header(line)
            qlen, qchars = _parse_block_postheader(next(ilines))
            # Pass control to the sequence parser
            sequences = list(_parse_sequences(ilines, qlen))
            # Validation
            if not len(sequences) == seqcount:
                logging.warn("Expected %d sequences in block, found %d",
                             seqcount, len(sequences))
            yield {'level': level,
                   'one': one,
                   'name': name,
                   # 'seqcount': seqcount,
                   'params': params,
                   'query_length': qlen,
                   'query_chars': qchars,
                   'sequences': sequences,
                   }


def _parse_sequences(ilines, expect_qlen):
    """Parse the sequences in the current block.

    Sequence looks like:

    $3=227(209):
    >gi|15606894|ref|NP_214275.1|  {|2(244)|<Aquificae(B)>}DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|2984127|gb|AAC07663.1| DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|75
    {()YVPFARKYRPKFFREVIGQEAPVRILKNAIKNDRVAHaYLFAGPRGVGKTTIARILAKALNcknpskgepcgecencreiDRGVFPDLIEMDAASNRGIDDVRA-LKEAVNYKPIKG-KYKVYIIDEAHMLTKEAFNALLKTLEEPPPRTVFVLCTTEYDKILPTILSRCQRIIFSKVRKEKVIEYLKKICEKEGIECEEGALEVLAHASEGCMRDAASLLDQASVYGE()}*

    """
    while True:
        first = next(ilines)
        if first.startswith('_') and first.endswith('].'):
            # End of sequences & end of block
            break

        index, this_len, query_len = _parse_seq_preheader(first)
        (rec_id, dbxrefs, headlen, taillen, phylum, taxchar, description
                ) = _parse_seq_header(next(ilines))
        headseq, molseq, tailseq = _parse_seq_body(next(ilines))

        # Validation
        if expect_qlen != query_len:
            logging.warn("Query length in %s given as %d; expected %d",
                         rec_id, query_len, expect_qlen)
        if not headseq and not headlen:
            headlen = 0
        if not tailseq and not taillen:
            taillen = 0
        if headseq:
            if headlen is None:
                headlen = len(headseq)
            elif headlen != len(headseq):
                logging.warn("Conflicting head flank lengths in %s: %d, %d",
                             rec_id, headlen, len(headseq))
        if tailseq:
            if taillen is None:
                taillen = len(tailseq)
            elif taillen != len(tailseq):
                logging.warn("Conflicting tail flank lengths in %s: %d, %d",
                             rec_id, taillen, len(tailseq))

        yield {'index': index,
               'id': rec_id,
               'description': description,
               'dbxrefs': dbxrefs,
               'phylum': phylum,
               'taxchar': taxchar,
               'head_len': headlen,
               'tail_len': taillen,
               'head_seq': headseq,
               'tail_seq': tailseq,
               'length': this_len,
               'seq': molseq,
               }


# Microparsing

def _parse_block_header(line):
    """
    [0_(1)=fa2cma(8){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:
    """
    level = line[1]
    one, _rest = line[4:].split(')=', 1)
    name, _rest = _rest.split('(', 1)
    seqcount, _rest = _rest.split(')', 1)
    params = _rest.strip('{}:')
    # try:
    #     params = dict((key, float(val))
    #             for key, val in (pair.split('=')
    #                 for pair in _rest[1:-2].split(',')))
    # except ValueError:
    #     # Couldn't convert params to key-val pairs, for whatever reason
    #     logging.warn("Failed to parse CMA params: %s", _rest[1:-2])
    #     params = {}
    return int(level), int(one), name, int(seqcount), params


def _parse_block_postheader(line):
    """
    (209)**************!*****************!!*************...
    """
    parts = line[1:].split(')', 1)
    qlen = int(parts[0])
    if not len(parts[1]) == qlen:
        logging.warn("postheader expected %d-long query, found %d",
                     qlen, len(parts[1]))
    return qlen, parts[1]


def _parse_seq_preheader(line):
    """
    $3=227(209):
    """
    match = re.match(r"\$ (\d+) = (\d+) \( (\d+) \):", line, re.VERBOSE)
    index, this_len, query_len = match.groups()
    return map(int, (index, this_len, query_len))


def _parse_seq_header(line):
    """Unique ID, head/tail lengths and taxonomy info from a sequence header.

    The description is the part of the FASTA/CMA sequence header starting after
    the first space (i.e. excluding ID), to the end of the line.

    This function looks inside the first '{...}' pair to extract info.

    Ex:

    >consensus seq
    >gi|15606894|ref|NP_214275.1|  {|2(244)|<Aquificae(B)>}DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|2984127|gb|AAC07663.1| DNA polymerase III gamma subunit [Aquifex aeolicus VF5] >gi|75
    >gi|3212262|pdb|1A2K|C  {<Chordata(M)>}Chain C, Gdpran-Ntf2 Complex >gi|3212263|pdb|1A2K|D Chain D, Gdpran-Ntf2 Complex >gi|3212264|pdb|1A2K|E Chain E, Gdpran-Ntf2 Complex >gi|5542273|pdb|1IBR|A C

    """
    # ENH: use the two functions in esbglib.parseutils
    # or, move one or both of those functions into here
    _parts = line[1:].split(None, 1)
    rec_id = _parts[0]
    descr = _parts[1] if _parts[1:] else ''
    # Database cross references
    dbxrefs = {}
    if '|' in rec_id:
        id_gen = iter(rec_id.rstrip('|').split('|'))
        for key in id_gen:
            try:
                dbxrefs[key] = next(id_gen)
            except StopIteration:
                break

    # Head/tail lengths and taxonomy codes
    headlen = taillen = None
    phylum = taxchar = ''
    if descr.startswith('{'):
        _deets, description = descr[1:].split('}', 1)
        match = re.search(r"""
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
        """, _deets, re.VERBOSE)
        if match:
            headlen, taillen, phylum, taxchar = match.groups()
            if headlen is not None:
                headlen = int(headlen)
            if taillen is not None:
                taillen = int(taillen)
            if phylum is None:
                phylum = ''
            if taxchar is None:
                taxchar = ''
        else:
            logging.warn("Couldn't match head/tail: %s", _deets)
    else:
        description = descr

    # TODO - return a dictionary here, update it in _parse_sequences
    return rec_id, dbxrefs, headlen, taillen, phylum, taxchar, description


def _parse_seq_body(line):
    """
    Ex:

    {()YVPFARKYRPKFFREVIGQEAPVRILKALNcknpskgepcgereiDRGVFPDVRA-LKLLDQASVYGE()}*
    MENINNI{()----------FKLILVGDGKFFSSSGEIIFNIWDTKFGGLRDGYYRLTYKNEDDM()}*
    """
    head, _rest = line.rstrip('*').split('{()', 1)
    molseq, tail = _rest.split('()}', 1)
    return (head, molseq, tail)


# --------------------------------------------------------------------
# Write

def write(blocks, outfile):
    if isinstance(blocks, dict) and 'sequences' in blocks:
        blocks = [blocks]
    with sugar.maybe_open(outfile, 'w+') as outstream:
        outstream.writelines(
                itertools.chain(*map(_format_block, blocks)))


def _format_block(block):
    # Block header
    # [0_(1)=structs.seq(12){go=19,gx=2,pn=5.0,lf=0,rf=0}:
    # (69)*****!***************************!!!******!******************!*******
    yield """\
[{level}_({one})={name}({seqcount}){{{params}}}:
({query_length}){query_chars}

""".format(
        seqcount=len(block['sequences']),
        # fixed_params=','.join('{0}={1}'.format(key, val)
        #                       for key, val in block['params'].iteritems()),
        **block)
    # Sequences
    for s in _format_sequences(block['sequences'], block['query_length']):
        yield s
    # Block close
    yield "_%s].\n" % block['level']


def _format_sequences(sequences, query_length):
    # $1=254(255):
    # >2QG5|A  {|27(6)|}
    # {()YTLENTIGRGSWGEVKIAVQKGTRIRRAAKKIPKYFV---EDVDRFKQEIEIMKSLDHPNIIRLYETFEDNTDIYLVMELCTGGELFERVVHKRVFRESDAARIMKDVLSAVAYCHKLNVAHRDLKPENFLFltdSPDSPLKLIDFGLAARFkpGKMMRTKVGTPYYVSPQVLEGL-YGPECDEWSAGVMMYVLLCGYPPFSAPTDXEVMLKIREGTFtfpeKDWLNVSPQAESLIRRLLTKSPKQRIT-----SLQALEHEW-()}*
    for idx, seq in enumerate(sequences):
        head_tail = ("|{head_len}({tail_len})|".format(**seq)
                    if seq['head_len'] or seq['tail_len'] else "")
        taxonomy=("<{phylum}({taxchar})>".format(**seq)
                    if seq['phylum'] and seq['taxchar'] else "")
        special_header = ("{{{0}{1}}}".format(head_tail, taxonomy)
                    if head_tail or taxonomy else "")
        seq['index'] = idx + 1
        yield """\
${index}={length}({0}):
>{id} {1}{description}
{head_seq}{{(){seq}()}}{tail_seq}*

""".format(query_length, special_header, **seq)


# --------------------------------------------------------------------
# Utilities

def realign_seqs(block, gap_char='.', align_indels=False):
    """Add gaps to a block so all residues in a column are equivalent.

    Given a block, containing a list of "sequences" (dicts) each containing a
    "seq" (actual string sequence, where upper=match, lower=insert, dash=gap),
    insert gaps (- or .) into the sequences s.t.

    1. columns line up properly, and
    2. all resulting sequences have the same length


    The reason this needs to be done is that the query/consensus sequence is not
    assigned gaps to account for inserts in the other sequences. We need to add
    the gaps back to obtain a normal alignment.

    `return`: a list of realigned sequence strings.
    """
    # ENH: align inserts using an external tool (if align_indels)
    all_chars = [list(sq['seq']) for sq in block['sequences']]
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
    return [''.join(row) for row in all_chars]


def consensus2block(record, level=0, name=None):
    """Convert a Biopython SeqRecord to a esbglib.cma block.

    Ungapping is handled here.
    """
    cons_ungap = str(record.seq).replace('-', '').replace('.', '').upper()
    record.seq = cons_ungap
    return dict(
            level=level, #record.annotations.get('level', 0),
            one=1,
            name=name or record.id,
            params='go=10000,gx=2000,pn=1000.0,lf=0,rf=0',
            query_length=len(cons_ungap),
            query_chars='*'*len(cons_ungap),
            sequences=[seqrecord2sequence(record, len(cons_ungap), 1)]
            )


def seqrecord2sequence(record, qlen, index):
    """Convert a Biopython SeqRecord to a esbglib.cma block.

    Indels (gaps, casing) must have already been handled in the record.
    """
    # aligned_len = sum(map(str.isupper, str(record.seq)))
    # assert qlen == aligned_len, (
    #         "Aligned sequence length %s, query %s\n%s"
    #         % (aligned_len, qlen, str(record)))

    description = (record.description.split(' ', 1)[1]
            if ' ' in record.description
            else ' ')
    return dict(index=index,
                id=record.id,
                description=description,
                dbxrefs={},
                phylum='',
                taxchar='',
                head_len=None,
                tail_len=None,
                head_seq='',
                tail_seq='',
                length=len(record.seq) - record.seq.count('-'),
                seq=str(record.seq),
                )

def replace_asterisks(seq, label=None):
    if '*' in seq:
        logging.warn("Sequence %scontains '*' character; replacing with 'X'",
                str(label) + ' ' if label else '')
    return str(seq).replace('*', 'X')


def collapse_to_consensus(seqrecords, strict=False, do_iron=True):
    """Opposite of realign_seqs.

    Input sequences should all be the same length.

    The first record must be the consensus.
    """
    level = 0
    name = seqrecords[0].id
    # If this is a CMA alignment, extract additional info:
    if hasattr(seqrecords, '_records'):
        if hasattr(seqrecords, 'level'):
            level = seqrecords.level
        if hasattr(seqrecords, 'name'):
            name = seqrecords.name
        seqrecords = seqrecords._records

    consensus = seqrecords.pop(0)
    cons_length = len(consensus)
    for i, s in enumerate(seqrecords):
        if len(s) != cons_length:
            raise ValueError(
            "Sequence #%d has length %d, consensus is %d"
            % (i+2, len(s), cons_length))

    if '.' in str(consensus.seq):
        # Strict -- error if there's a '-'
        if '-' in str(consensus.seq):
            if strict:
                raise ValueError("Consensus contains '-' gap characters")
            logging.warn("Consensus sequence contains both '.' and '-' gap "
                         "characters -- is it really the consensus?")
            aligned_cols = [(c not in '.-') for c in str(consensus.seq)]
        else:
            aligned_cols = [c != '.' for c in str(consensus.seq)]
    else:
        # A little more ambiguous...
        aligned_cols = [c != '-' for c in str(consensus.seq)]

    consensus.seq = replace_asterisks(consensus.seq, 'consensus')
    # Start a block with the consensus sequence
    block = consensus2block(consensus, level=level, name=name)
    qlen = block['query_length']

    # Collapse & add remaining sequences to the block
    for index, rec in zip(xrange(2, len(seqrecords)+2), seqrecords):
        # Collapse rec.seq down to aligned size
        new_mol_seq = []
        is_beginning = True
        for aligned_col, char in zip(aligned_cols,
                                     replace_asterisks(rec.seq, index)):
            if aligned_col:
                is_beginning = False
                if char in '-.':
                    # deletion
                    new_mol_seq.append('-')
                else:
                    # aligned character
                    new_mol_seq.append(char.upper())
            else:
                # it's an insert or nothing
                # (also, skip any left-side inserts)
                if char not in '-.' and not is_beginning:
                    new_mol_seq.append(char.lower())
        rec.seq = ''.join(new_mol_seq)
        if do_iron:
            rec.seq = iron(rec.seq)
        block['sequences'].append(seqrecord2sequence(rec, qlen, index))

    return block


def iron(sequence):
    """'Iron out' indel regions in the aligned sequence.

    Any inserts next to deletions are converted to matches (uppercase).

    Given a CMA string like:
        AAAAbc--de-f--gAAA
    Result:
        AAAABCDEFgAAA
    """
    r_indel = re.compile(r'(-[a-y]|[a-y]-)')
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
                    outchars.append('-' * seen_gaps)
                    seen_gaps = 0
                    outchars.append(char)
                else:
                    # Convert a preceding indel char to a 'match' (uppercase)
                    # If the indel and gap are both multiple chars, this will
                    # capitalize the insert left-to-right, then leave any gap
                    # remainer as-is.
                    assert char == '-'
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
                        outchars.append('-' * seen_gaps)
                        seen_gaps = 0
                    inserts = [char]
                elif char.isupper():
                    in_gap = False
                    # End of the gap -- emit
                    if inserts:
                        outchars.extend(inserts)
                        inserts = []
                    outchars.append('-' * seen_gaps)
                    seen_gaps = 0
                    outchars.append(char)
                else:
                    # Extend the gap
                    assert char == '-'
                    seen_gaps += 1

            else:
                assert not inserts and not seen_gaps, (
                    "Inserts: %s, gaps: %s, seq: %s, in_ins=%s, in_gap=%s"
                    % (inserts, seen_gaps, sequence, in_insert, in_gap))
                # Coming from Match state
                if char.isupper():
                    # Extend the match
                    outchars.append(char)
                elif char.islower():
                    inserts.append(char)
                    in_insert = True
                else:
                    assert char == '-'
                    seen_gaps += 1
                    in_gap = True

        # Emit any trailing indel
        if inserts:
            outchars.extend(inserts)
        if seen_gaps:
            outchars.append('-' * seen_gaps)
        sequence = ''.join(outchars)
        # logging.info(sequence)
        assert (sequence.replace('-', '').upper()
                ==
                orig_sequence.replace('-', '').upper()), \
                '\nOrig: ' + orig_sequence + \
                '\nIron: ' + sequence
    return sequence


# --------------------------------------------------------------------

if __name__ == '__main__':
    # Test
    import sys
    from .utils import get_equivalent_positions

    cmafiles = sys.argv[1:]
    if not cmafiles:
        sys.exit("usage: cma.py <cmafile1> [<cmafile2> ...]")

    for fname in cmafiles:
        block = next(parse(fname))
        print block['query_length'], "aa query"
        print len(block['sequences']), "sequences"
        print "of lengths:",
        for seq in block['sequences']:
            print "%d/%d(%d)" % (
                len(seq['seq']), seq['length'], seq['seq'].count('-')),
        print
        print "  Equivalencies:"
        for idx in xrange(1, 20):
            print idx, '=', get_equivalent_positions(block)[idx]
        print
        print "  Writing the file back out:"
        print
        if len(block['sequences']) < 60:
            write(block, sys.stdout)

