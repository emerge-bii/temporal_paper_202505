#!/bin/bash

# Install the conda environment:
conda env create -f temporal_paper.yml --force

source $(conda info --base)/etc/profile.d/conda.sh
# Since the name changes with each version, make sure to pull the version
# name out of the yml file
ENV_NAME=$(awk '/name: temporal_paper_v/{print $NF}' temporal_paper.yml)
conda activate $ENV_NAME

# Install doMPI package in R (we use the cloud-0 mirror so that the download
# will be fast in any part of the world)
Rscript -e 'install.packages("doMPI", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("iCAMP", repos="https://cloud.r-project.org")' No longer necessary
Rscript -e 'devtools::install_github("m-jahn/fluctuator")'
Rscript -e 'install.packages("analogue", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("fastshap", repos="https://cloud.r-project.org")'

# Deactivate the conda environment
conda deactivate
