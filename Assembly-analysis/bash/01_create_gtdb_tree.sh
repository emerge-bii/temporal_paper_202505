#!/bin/bash
# 0: Set up variables and initialize directories

PROJ_DIR=/mnt/home/ernakovich/heh1030/EMERGE/00_ref_data/Phylogenies/GTDB_TK_Tree
MAGS_DIR=/mnt/home/ernakovich/heh1030/EMERGE/00_ref_data/Emerge_MAGs_v1 # location of the MAG db
ARCH_GENOMES_LIST=$PROJ_DIR/archaea.txt
BAC_GENOMES_LIST=$PROJ_DIR/bacteria.txt                                     
ARCH_SYM_LNKS=$PROJ_DIR/Archaea_genomes
BAC_SYM_LNKS=$PROJ_DIR/Bacteria_genomes
ARCH_TREE_DIR=$PROJ_DIR/Archaea_Tree
BAC_TREE_DIR=$PROJ_DIR/Bacteria_Tree


# Directories to hold trees
mkdir -p $ARCH_TREE_DIR
mkdir -p $BAC_TREE_DIR


# 1: Create a EMERGE MAGs Sym link directory for archaea and bacterial MAGs
echo "Creating directories and gathering genomes into sym links"
mkdir -p $ARCH_SYM_LNKS
mkdir -p $BAC_SYM_LNKS

# Archaeal sym links
cat $ARCH_GENOMES_LIST | while read genome; 
   do
     find $MAGS_DIR -name $genome -print -exec ln -s {} $ARCH_SYM_LNKS ";"
   done

# Bacterial Sym links
cat $BAC_GENOMES_LIST | while read genome;
   do
     find $MAGS_DIR -name $genome -print -exec ln -s {} $BAC_SYM_LNKS ";"
   done


# 2: run gtdb denovo workflow

# On bacteria
echo "Running gtdb-tk de novo work flow on Bacteria"
# note: I choose to root at Patescibacteria b/c this is the CPR radiation
# and Zhu et al. 2019 suggests that the archaeal and bacterial branch out from here.
gtdbtk de_novo_wf \
   --genome_dir $BAC_SYM_LNKS \
   --bacteria \
   --outgroup_taxon p__Patescibacteria \
   --out_dir $BAC_TREE_DIR \
   --cpus 24

# On Archaea
echo "Running gtdb-tk de novo work flow	on Archaea"
gtdbtk de_novo_wf \
   --genome_dir $ARCH_SYM_LNKS \
   --archaea \
   --outgroup_taxon p__Altarchaeota \
   --out_dir $ARCH_TREE_DIR \
   --cpus 24

# NOTE: The gtdb newick format cannot be read into R as is. You'll need to use
# FigTree (available in Bioconda) to view and export a new newick tree version that can be
# read in in R

# FigTree can be installed like this:
# conda create -n figtree
# conda activate figtree
# conda install -c bioconda -c conda-forge figtree

