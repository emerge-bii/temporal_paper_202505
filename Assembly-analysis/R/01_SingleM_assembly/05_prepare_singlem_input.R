#' ## Interpret Assembly Analysis results
#' This step takes the output from 01_calc_betaNTI.R and 02_calc_rcbc.R and runs initial figures and analyses on them.

#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = FALSE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')
#'
#+ include=TRUE
# Loading necessary packages and data
# Load DADA2 and required packages
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(cowplot)
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("Assembly-analysis", "outputs")
figures.fp <- here("Assembly-analysis", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### =================================================================================================== ####
# User defined functions
#### =================================================================================================== ####
# Reads in a newick format tree and roots it with the midpoint root
read_singlem_tree <- function(tree_file_path, root = TRUE, 
                              visualize = FALSE) {
  require(ape)
  require(phytools)
  
  tree.raw <- phytools::read.newick(tree_file_path)
  
  # check tree merge
  if(visualize) {
    writeLines("Plotting initial tree")
    plot(tree.raw)
  }
  
  # Root the tree
  if(root){
    tree.rt <- phytools::midpoint.root(tree.raw)
  } else {
    tree.rt <- tree.raw
  }
  
  
  # check tree
  if(visualize) {
    writeLines("Plotting rooted tree")
    # plot resulting tree
    plot(tree.rt)
  }
  return(tree.rt)
}

read_singlem_taxonomy <- function(tax_file_path, DistinctTally = FALSE) {
  tax_table <- read_csv(tax_file_path) %>%
    select(OTU, taxonomy) %>%
    separate(taxonomy, 
             into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = "; ") %>%
    select(OTU, Root, Domain:Species)
  
  # Replace conflicting taxonomy levels with NAs for each OTU
  # Essentially, taxonomy is the highest level common across all taxa
  writeLines("Collapsing taxonomy by least common method")
  tax_table <- tax_table %>%
    group_by(OTU, Root) %>%
    mutate(DomainDistinct = n_distinct(Domain),
           Domain = ifelse(DomainDistinct != 1, NA, Domain)) %>% 
    mutate(across(Domain:Species, ~ifelse(is.na(Domain), NA, .x))) %>%
    ungroup() %>%
    group_by(OTU, Domain) %>%
    mutate(PhylumDistinct = n_distinct(Phylum),
           Phylum = ifelse(PhylumDistinct != 1, NA, Phylum)) %>% # if there are mismatches in phylum within the same OTU, make it NA
    mutate(across(Phylum:Species, ~ifelse(is.na(Phylum), NA, .x))) %>% # Make any otu that is NA in Phylum, NA across it's taxonomy string below phylum
    ungroup() %>%
    group_by(OTU, Phylum) %>%
    mutate(ClassDistinct = n_distinct(Class),
           Class = ifelse(ClassDistinct != 1, NA, Class)) %>%
    mutate(across(Class:Species, ~ifelse(is.na(Class), NA, .x))) %>%
    ungroup() %>%
    group_by(OTU, Class) %>%
    mutate(OrderDistinct = n_distinct(Order),
           Order = ifelse(OrderDistinct != 1, NA, Order)) %>%
    mutate(across(Order:Species, ~ifelse(is.na(Order), NA, .x))) %>%
    ungroup() %>%
    group_by(OTU, Order) %>%
    mutate(FamilyDistinct = n_distinct(Family),
           Family = ifelse(FamilyDistinct != 1, NA, Family)) %>%
    mutate(across(Family:Species, ~ifelse(is.na(Family), NA, .x))) %>%
    ungroup() %>%
    group_by(OTU, Family) %>%
    mutate(GenusDistinct = n_distinct(Genus),
           Genus = ifelse(GenusDistinct != 1, NA, Genus)) %>%
    mutate(across(Genus:Species, ~ifelse(is.na(Genus), NA, .x))) %>%
    ungroup() %>%
    group_by(OTU, Genus) %>%
    mutate(SpeciesDistinct = n_distinct(Species),
           Species = ifelse(SpeciesDistinct != 1, NA, Species)) %>%
    ungroup() %>%
    group_by(OTU) %>%
    distinct()
  
  if(DistinctTally == F) {
    tax_table <- tax_table %>%
      select(!contains("Distinct"))
  }

  return(tax_table)
}

# Check for matching names between datasets
check_sample_metadata_singlem <- function(otu_table, sample_metadata, 
                                  sample_id_column = "temporal_sample_id",
                                  taxonomy = NULL,
                                  tree = NULL) {
  # Checking OTU table and sample metadata
  # Find the missing samples for otu_table and sample metadata
  missing_genomic_data <- na.omit(sample_metadata$temporal_sample_id[!(sample_metadata$temporal_sample_id %in% names(otu_table[-1]))]) # Samples missing from OTU table
  writeLines(paste("Samples missing from the OTU table that are present in the",
                   "metadata:", paste(missing_genomic_data, collapse = ", ")))
  
  missing_env_data <- names(otu_table[-1])[!(names(otu_table[-1]) %in% sample_metadata$temporal_sample_id )] # Samples missing from env data
  writeLines(paste("Samples missing from the metadata that are present in the",
                   "OTU table:", paste(missing_env_data, collapse = ", ")))
  writeLines(paste0("Now filtering out missing sample(s)..."))
  
  # Filter out the samples that don't occur in both
  otu_table <- otu_table %>% select(OTU, !all_of(missing_env_data))
  
  sample_metadata <- sample_metadata %>% 
    filter((!!as.name(sample_id_column) %in% names(otu_table[-1])))
  
  # After filtering samples, remove any columns which are entirely NAs (or 0s for MAGs)
  sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]
  
  # Identify env_data columns that are all NA
  all_na_cols <- names(sample_metadata[,colSums(is.na(sample_metadata))==nrow(sample_metadata)])
  # filter out metadata variables that are all NAs
  sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]
  
  
  # Reorder otu table columns to match order in sample metadata
  otu_table <- otu_table[,c("OTU", sample_metadata$temporal_sample_id)]
  
  writeLines(paste(nrow(sample_metadata), "samples in sample_metadata"))
  writeLines(paste(ncol(otu_table[-1]), "samples in otu_table"))
  
  
  # Checking taxonomy, tree, and otu table (if taxonomy and tree are present)
  # and reordering to match tree order
  if(!is.null(tree) & !is.null(taxonomy)) {
    # reorder otu table to match order of tips in tree
    match.phylo.otu <- picante::match.phylo.data(tree, 
                                                 otu_table %>%
                                                   column_to_rownames(var = "OTU"))
    tree <- match.phylo.otu$phy
    otu_table <- match.phylo.otu$data %>%
      rownames_to_column(var = "OTU") %>%
      mutate(rowname = OTU) %>%
      column_to_rownames()
    
    # Check that otu table OTUs are all in the tree and taxonomy file
    OTUs_missing_from_taxonomy <- otu_table$OTU[!otu_table$OTU %in% taxonomy$OTU]
    writeLines(paste("OTUs missing from the taxonomy that are present in the",
                     "OTU table:", paste(OTUs_missing_from_taxonomy, collapse = ", ")))
    
    OTUs_missing_from_tree <- otu_table$OTU[!otu_table$OTU %in% tree$tip.label]
    writeLines(paste("OTUs missing from the tree that are present in the",
                     "OTU table:", paste(OTUs_missing_from_tree, collapse = ", ")))
    
    writeLines(paste0(nrow(otu_table), " taxa in otu table"))
    writeLines(paste0(length(tree$tip.label), " tips in tree"))
    writeLines(paste0(nrow(taxonomy), " taxa in taxonomy"))
    
    # reorder taxonomy to match otu_table and tree tip order
    taxonomy <- taxonomy %>% column_to_rownames("OTU")
    taxonomy <- taxonomy[rownames(match.phylo.otu$data),] %>%
      rownames_to_column(var = "OTU") %>%
      mutate(rowname = OTU) %>%
      column_to_rownames()
    
  }
  
  
  # Create list of outputs to return:
  return_list <- list(sample_metadata = sample_metadata, otu_table = otu_table)
  
  if(!is.null(taxonomy)) {
    return_list$taxonomy <- taxonomy
  }
  
  if(!is.null(tree)) {
    return_list$tree <- tree
  }
  
  return(return_list)
}

#### =================================================================================================== ####
#### Read in data 
#### =================================================================================================== ####
# Read in OTU table and env data
singlem_otu_S3.11_otu <- read_csv(here(outputs.fp, "singlem_otu_S3.11.rel_abund.csv"))
singlem_otu_S3.11_taxa <- read_singlem_taxonomy(here(outputs.fp, "singlem_otu_S3.11.taxonomy.csv"))

# Read in singlem Tree
singlem_otu_S3.11_tree <- read_singlem_tree(here(outputs.fp, "S3.11.repset_aln.tre"), root = FALSE)


S3.11_input_ra <- check_sample_metadata_singlem(otu_table = singlem_otu_S3.11_otu,
                                  sample_metadata = input_ra$sample_metadata,
                                  taxonomy = singlem_otu_S3.11_taxa,
                                  tree = singlem_otu_S3.11_tree)

# Read in OTU table and env data
singlem_otu_S3.5_otu <- read_csv(here(outputs.fp, "singlem_otu_S3.5.rel_abund.csv"))
#singlem_otu_S3.5_otu_count <- read_csv(here(outputs.fp, ""))
singlem_otu_S3.5_taxa <- read_singlem_taxonomy(here(outputs.fp, "singlem_otu_S3.5.taxonomy.csv"))

# Read in singlem Tree
singlem_otu_S3.5_tree <- read_singlem_tree(here(outputs.fp, "S3.5.repset_aln.tre"), root = FALSE)


S3.5_input_ra <- check_sample_metadata_singlem(otu_table = singlem_otu_S3.5_otu,
                                                sample_metadata = input_ra$sample_metadata,
                                                taxonomy = singlem_otu_S3.5_taxa,
                                                tree = singlem_otu_S3.5_tree)

#### =================================================================================================== ####

#### Compare filtered taxa data to unfiltered
#### =================================================================================================== ####
singlem_otu_S3.11_taxa.dist <- read_singlem_taxonomy(here(outputs.fp, "singlem_otu_S3.11.taxonomy.csv"),
                                                     DistinctTally = T)

DistinctTally.S3.11 <- singlem_otu_S3.11_taxa.dist %>%
  ungroup() %>%
  select(ends_with("Distinct")) %>%
  mutate(OTUDistinct = 1) %>%
  summarize(across(ends_with("Distinct"), ~sum(.))) %>%
  mutate(across(ends_with("Distinct"), ~.x-OTUDistinct, .names = "NumNonDist_{.col}")) %>%
  mutate(across(!contains("NumNon"), ~100*(.x-OTUDistinct)/OTUDistinct, 
                .names = "PercNonDist_{.col}")) %>%
  pivot_longer(cols = contains("_")) %>%
  separate(name, into = c("NonDist", "Dist")) %>%
  select(!ends_with("Distinct")) %>%
  pivot_wider(names_from = "Dist", values_from = "value")  

#DistinctTally.S3.11_seq <- singlem_otu_S3.11_taxa.dist %>%
#  left_join(singlem_otu_S3.11.wf$taxa_table %>% 
#              select(OTU, sequence) %>%
#              distinct(), by = "OTU")
  

singlem_otu_S3.5_taxa.dist <- read_singlem_taxonomy(here(outputs.fp, "singlem_otu_S3.5.taxonomy.csv"),
                                                     DistinctTally = T)

DistinctTally.S3.5 <- singlem_otu_S3.5_taxa.dist %>%
  ungroup() %>%
  select(ends_with("Distinct")) %>%
  mutate(OTUDistinct = 1) %>%
  summarize(across(ends_with("Distinct"), ~sum(.))) %>%
  mutate(across(ends_with("Distinct"), ~.x-OTUDistinct, .names = "NumNonDist_{.col}")) %>%
  mutate(across(!contains("NumNon"), ~100*(.x-OTUDistinct)/OTUDistinct, 
                .names = "PercNonDist_{.col}")) %>%
  pivot_longer(cols = contains("_")) %>%
  separate(name, into = c("NonDist", "Dist")) %>%
  select(!ends_with("Distinct")) %>%
  pivot_wider(names_from = "Dist", values_from = "value")  

#DistinctTally.S3.5_seq <- singlem_otu_S3.5_taxa.dist %>%
#  left_join(singlem_otu_S3.5.wf$taxa_table %>% 
#              select(OTU, sequence) %>%
#              distinct(), by = "OTU")
 
#### =================================================================================================== ####
#### Singleton analyses
#### =================================================================================================== ####
S3.11_specaccum <- specaccum(t(S3.11_input_ra$otu_table[,-1]))


#### =================================================================================================== ####
# Save outputs
#### =================================================================================================== ####
saveRDS(S3.11_input_ra, paste0(outputs.fp, "/S3.11_input_ra.RDS"))

saveRDS(S3.5_input_ra, paste0(outputs.fp, "/S3.5_input_ra.RDS"))

#write.csv(DistinctTally.S3.11_seq,
#          paste0(outputs.fp, "/singlem_otu_S3.11.taxa_distinct.csv"),
#          quote=FALSE, row.names = FALSE)

#write.csv(DistinctTally.S3.5_seq,
#          paste0(outputs.fp, "/singlem_otu_S3.5.taxa_distinct.csv"),
#          quote=FALSE, row.names = FALSE)
#### =================================================================================================== ####
