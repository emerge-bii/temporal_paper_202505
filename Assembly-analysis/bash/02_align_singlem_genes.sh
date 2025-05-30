#!/bin/bash

# This script creates two phylogenetic tree from two of the singlm marker genes to compare them
# to the gtdb tree created from the MAGs

# Step 0: Identify variables (assumes working directory is temporal_paper/Assembly-analysis/bash)
DATA_DIR=../outputs
OUTPUT_DIR=../outputs

# Step 1: Align sequences using mafft
mafft --auto --thread 24 $DATA_DIR/S3.5.repset.fasta > $OUTPUT_DIR/S3.5.repset_aln.fasta

mafft --auto --thread 24 $DATA_DIR/S3.11.repset.fasta > $OUTPUT_DIR/S3.11.repset_aln.fasta

# Step 2: Create tree using fasttree
# We will use -gtr -gamma options

FastTree -nt -gtr -gamma $OUTPUT_DIR/S3.5.repset_aln.fasta > $OUTPUT_DIR/S3.5.repset_aln.tre

FastTree -nt -gtr -gamma $OUTPUT_DIR/S3.11.repset_aln.fasta > $OUTPUT_DIR/S3.11.repset_aln.tre  



