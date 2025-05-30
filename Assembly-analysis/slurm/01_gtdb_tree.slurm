#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## NTasks is not thread count, be sure to leave it set at 1
#SBATCH --ntasks=1

## If your program will be using less than 24 threads, or you
## require more than 24 threads, set cpus-per-task to the 
## desired threadcount.  Leave this commented out for the
## default 24 threads.
##SBATCH --cpus-per-task=2

## You will need to specify a minimum amount of memory in the
## following situaitons:
##   1. If you require more than 128GB of RAM, specify either:
##      a. "--mem=512000" for at least 512GB of RAM (6 possible nodes)
##      b. "--mem=1000000" for at least 1TB of RAM (2 possible nodes)
##   2. If you are running a job with less than 24 threads, you will
##      normally be given your thread count times 5.3GB in RAM.  So
##      a single thread would be given about 5GB of RAM.  If you
##      require more, please specify it as a "--mem=XXXX" option,
##      but avoid using all available RAM so others may share the node.
#SBATCH --mem=512000

## Normally jobs will restart automatically if the cluster experiences
## an unforeseen issue.  This may not be desired if you want to retain
## the work that's been performed by your script so far.   
## --no-requeue

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="gtdb_denovo"
#SBATCH --output=01_gtdb_denovo.output

## Load the appropriate modules first.  Linuxbrew/colsa contains most
## programs, though some are contained within the anaconda/colsa
## module.  Refer to http://premise.sr.unh.edu for more info.
module purge
module load anaconda/colsa
conda find gtdbtk-2.1.0 -a

## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.


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
   --outgroup_taxon p__Altiarchaeota \
   --out_dir $ARCH_TREE_DIR \
   --cpus 24

# NOTE: The gtdb newick format cannot be read into R as is. You'll need to use
# FigTree (available in Bioconda) to view and export a new newick tree version that can be
# read in in R

# FigTree can be installed like this:
# conda create -n figtree
# conda activate figtree
# conda install -c bioconda -c conda-forge figtree
