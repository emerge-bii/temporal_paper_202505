#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## NTasks is not thread count, be sure to leave it set at 1
#SBATCH --ntasks=96

## If your program will be using less than 24 threads, or you
## require more than 24 threads, set cpus-per-task to the 
## desired threadcount.  Leave this commented out for the
## default 24 threads.
#SBATCH --cpus-per-task=1

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
##SBATCH --mem=512000

## Normally jobs will restart automatically if the cluster experiences
## an unforeseen issue.  This may not be desired if you want to retain
## the work that's been performed by your script so far.   
#SBATCH --no-requeue

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="calculate_RF_assembly.R"
#SBATCH --output=11_calculate_RF_assembly_%j.output

#SBATCH --exclude=node117

## Load the appropriate modules first.  Linuxbrew/colsa contains most
## programs, though some are contained within the anaconda/colsa
## module.  Refer to http://premise.sr.unh.edu for more info.
module purge
module load anaconda/colsa
eval "$(conda shell.bash hook)"
#source /mnt/home/ernakovich/shared/software/miniconda3/bin/activate
conda activate /mnt/home/ernakovich/shared/software/miniconda3/envs/temporal_paper_v0.14


## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.



# 0: Set up variables and initialize directories

PROJ_DIR=/mnt/home/ernakovich/heh1030/EMERGE/decadal/temporal_paper/Assembly-analysis/

# move into project directory
cd $PROJ_DIR
# Make Rout dir if it doesnt' exist
mkdir -p R/Rout

echo "Output will be in"
echo $PROJ_DIR/R/Rout

# Run the R script
#Rscript R/00_all_MAGs_assembly/05_calculate_rf_assembly.R
mpiexec -n  1 Rscript R/00_all_MAGs_assembly/05_calculate_rf_assembly.R --show-progress
