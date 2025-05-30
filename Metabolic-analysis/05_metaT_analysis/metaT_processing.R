library(tidyverse)
library(here)
# Set memory_saving to FALSE
source(here("setup.R"))

# Update directory here first
metaT_processed_directory <- here("data", "Emerge_metaTs_processed_v7")

# See setup.R for functions
metaT_data <- raw_sample_metadata %>%
  read_metaT_data()

metaT_data %>%
  apply_pathway_labels_by_ko(pathway_definitions, cazy_definitions) %>%
  write_tsv(here(metaT_processed_directory, "metaT_pathways.tsv"))

metaT_data %>%
  complete(genome, SampleID__, fill = list(tpm = 0)) %>%
  group_by(genome, SampleID__) %>%
  summarise(tpm = sum(tpm)) %>%
  write_tsv(here(metaT_processed_directory, "metaT_genomes.tsv"))
