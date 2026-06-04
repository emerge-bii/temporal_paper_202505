# GraftM analysis scripts for the EMERGE temporal paper

GraftM analysis pipeline for distinguishing between genes with different functions but the same KO (e.g. pmoA vs amoA).

## Workflow

Inputs:
Seed sequences (ideally with experimentally confirmed function)
Emerge DRAM sequences (matching KO)

1. mmseqs search seed sequences against Uniref90 (max 300 sequences per seed)
2. GraftM create package (including Emerge DRAM sequences)
3. Manual annotation in arb of functional clades based on experimentally confirmed sequences
4. GraftM create package (using annotated tree)
5. Classify Emerge sequences using new GraftM package
6. Encode function as KO id + 0/1 (see below)
7. Replace KO id with functional KO id for downstream analyses

## New encodings

Final gene counts across EMERGE MAGs in brackets.

narG: K003700 (76)
nxrA: K003701 (2)
narH: K003710 (65)
nxrB: K003711 (1)
norB: K045610 (222)
norZ: K045611 (0)
pmoA: K109440 (12)
amoA: K109441 (1)
pmoB: K109450 (12)
amoB: K109451 (1)
pmoC: K109460 (22)
amoC: K109461 (2)
dsrA: K111800 (47)
rdsrA: K111801 (18)
dsrB: K111810 (45)
rdsrB: K111811 (19)
