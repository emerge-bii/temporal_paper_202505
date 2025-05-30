# Phylofactorization


library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)
library(ggrepel)
library(ape)
library(tidytree) # For tree plotting
library(ggtree)
library(ggnewscale)
library(ggstance)
library(viridis)
library(phylofactor)


# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("Assembly-analysis", "outputs")
figures.fp <- here("Assembly-analysis", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### ====================================================================== ####

#### Define some helpful functions:
#### ====================================================================== ####
#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Read in input
input <- input_counts

#### ====================================================================== ####
# 


#### ====================================================================== ####
data("FTmicrobiome")

OTUTable <- FTmicrobiome$OTUTable   #Our OTU table. Rows are OTUs and columns are samples
OTUTable[1:5,1:5]


body.site <- FTmicrobiome$X         #Our independent variable

#### ====================================================================== ####
#### ====================================================================== ####
# Too make this faster we will subset this to only include those in the lower quartile of CV
CV_dat_new <- CV_dat %>%
  filter(Lower_quartile == TRUE)
# Change tree genome names so that numbers start with X
tree <- keep.tip(input_ra$tree, unique(CV_dat_new$genome))
tree$tip.label <- gsub("^([0-9])(.{1,})", "X\\1\\2",tree$tip.label)



Data <- CV_dat_new %>%
  select(genome, Habitat, CV_recip, Depth) %>%
  pivot_wider(id_cols = genome, names_from = c("Depth", "Habitat"), 
              values_from = "CV_recip", values_fill = 0) %>%
  mutate(genome = gsub("^([0-9])(.{1,})", "X\\1\\2",genome)) %>%
  column_to_rownames(var = "genome")
Data <- Data[match(rownames(Data), tree$tip.label),] # reorder to match tree

# Prepare Independant variable
depth.habitat <- names(Data)
depth.habitat <- factor(depth.habitat, levels = depth.habitat)

depth.habitat_sep <- data.frame(depth.habitat = depth.habitat) %>%
  separate(depth.habitat, sep = "_", into = c("Depth", "Habitat"), remove = F) %>%
  mutate(Depth = factor(Depth, levels = c("0-9", "10-19", "20-29", "30-39")),
         Habitat = factor(Habitat, levels =c("Palsa", "Bog", "Fen")))



test.phyfac <- PhyloFactor(Data,
                           tree = tree,
                           depth.habitat_sep,
                           frmla = Data ~ Depth*Habitat, choice = "var", 
                           ncores = 7)



CV_recip <- CV_dat$CV_recip
names(CV_recip) <- CV_dat$genome

CV_recip <- CV_recip[tree$tip.label]
CV_dat[]



