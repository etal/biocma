BioCMA
======

Consensus Multiple Alignment format, using Biopython alignments.

I/O support and relevant functionality for the Consensus Alignment Format (CMA).

Biopython objects and conventions are used where possible.

Features of the CMA format:

- Sequence headers indicate #head and #tail residues
- Alignments are pairwise versus the first "consensus" sequence (consensus
  generally has no gaps)
- Insertions are lowercase; deletions are dashes

