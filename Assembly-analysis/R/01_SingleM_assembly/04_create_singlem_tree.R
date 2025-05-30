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

#### Read in data 
#### =================================================================================================== ####
# Read in OTU table and env data
input <- input_ra
input$otu_table <- input$otu_table[-1]

# Turn Habitat into a factor
input$sample_metadata <- input$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, 
                            levels = c("Palsa", "Collapsed Palsa", "Bog", "Fen")))

# Filter singlm otu talbes
singlem_otu_S3.11 <- singlem_otu_tables %>% 
  filter(grepl("S3.11", gene))
singlem_otu_S3.5 <- singlem_otu_tables %>% 
  filter(grepl("S3.5.ribosomal", gene))


# calc_taxonomic_conflicts <- function(singlem_otu_table = singlem_otu_S3.11, DistinctTally = TRUE, remove_NAs = TRUE, nocollapse = FALSE) {
#   tax_table <- singlem_otu_table %>% 
#     group_by(sequence) %>%
#     mutate(OTU = paste0(gene, "_OTU_", cur_group_id())) %>%
#     separate(taxonomy, 
#              into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
#              sep = "; ") %>%
#     select(OTU, Root, Domain:Species, everything())
#   
#   # Handle taxonomic disagreements - by default we collapse the taxonomy
#   if(nocollapse == FALSE) {
#     # Replace conflicting taxonomy levels with NAs for each OTU
#     # Essentially, taxonomy is the highest level common across all taxa
#     writeLines("Collapsing taxonomy by least common method")
#     tax_table <- tax_table %>%
#       group_by(OTU, Root) %>%
#       mutate(DomainDistinct = n_distinct(Domain, na.rm = remove_NAs),
#              Domain = ifelse(DomainDistinct != 1, NA, Domain)) %>% 
#       mutate(across(Domain:Species, ~ifelse(is.na(Domain), NA, .x))) %>%
#       ungroup() %>%
#       group_by(OTU, Domain) %>%
#       mutate(PhylumDistinct = n_distinct(Phylum, na.rm = remove_NAs),
#              Phylum = ifelse(PhylumDistinct != 1, NA, Phylum)) %>% # if there are mismatches in phylum within the same OTU, make it NA
#       mutate(across(Phylum:Species, ~ifelse(is.na(Phylum), NA, .x))) %>% # Make any otu that is NA in Phylum, NA across it's taxonomy string below phylum
#       ungroup() %>%
#       group_by(OTU, Phylum) %>%
#       mutate(ClassDistinct = n_distinct(Class, na.rm = remove_NAs),
#              Class = ifelse(ClassDistinct != 1, NA, Class)) %>%
#       mutate(across(Class:Species, ~ifelse(is.na(Class), NA, .x))) %>%
#       ungroup() %>%
#       group_by(OTU, Class) %>%
#       mutate(OrderDistinct = n_distinct(Order, na.rm = remove_NAs),
#              Order = ifelse(OrderDistinct != 1, NA, Order)) %>%
#       mutate(across(Order:Species, ~ifelse(is.na(Order), NA, .x))) %>%
#       ungroup() %>%
#       group_by(OTU, Order) %>%
#       mutate(FamilyDistinct = n_distinct(Family, na.rm = remove_NAs),
#              Family = ifelse(FamilyDistinct != 1, NA, Family)) %>%
#       mutate(across(Family:Species, ~ifelse(is.na(Family), NA, .x))) %>%
#       ungroup() %>%
#       group_by(OTU, Family) %>%
#       mutate(GenusDistinct = n_distinct(Genus, na.rm = remove_NAs),
#              Genus = ifelse(GenusDistinct != 1, NA, Genus)) %>%
#       mutate(across(Genus:Species, ~ifelse(is.na(Genus), NA, .x))) %>%
#       ungroup() %>%
#       group_by(OTU, Genus) %>%
#       mutate(SpeciesDistinct = n_distinct(Species, na.rm = remove_NAs),
#              Species = ifelse(SpeciesDistinct != 1, NA, Species)) %>%
#       ungroup() %>%
#       group_by(OTU) %>%
#       distinct()
#   } else {
#     # Identify conflicting taxonomy levels but do not collapse with
#     # NAs
#     writeLines("Identifying taxonomic conflicts")
#     tax_table <- tax_table %>%
#       group_by(OTU, Root) %>%
#       mutate(DomainDistinct = n_distinct(Domain, na.rm = remove_NAs),
#              Domain = ifelse(DomainDistinct != 1, NA, Domain)) %>% 
#       ungroup() %>%
#       group_by(OTU, Domain) %>%
#       mutate(PhylumDistinct = n_distinct(Phylum, na.rm = remove_NAs),
#              Phylum = ifelse(PhylumDistinct != 1, NA, Phylum)) %>% # if there are mismatches in phylum within the same OTU, make it NA
#       ungroup() %>%
#       group_by(OTU, Phylum) %>%
#       mutate(ClassDistinct = n_distinct(Class, na.rm = remove_NAs),
#              Class = ifelse(ClassDistinct != 1, NA, Class)) %>%
#       ungroup() %>%
#       group_by(OTU, Class) %>%
#       mutate(OrderDistinct = n_distinct(Order, na.rm = remove_NAs),
#              Order = ifelse(OrderDistinct != 1, NA, Order)) %>%
#       ungroup() %>%
#       group_by(OTU, Order) %>%
#       mutate(FamilyDistinct = n_distinct(Family, na.rm = remove_NAs),
#              Family = ifelse(FamilyDistinct != 1, NA, Family)) %>%
#       ungroup() %>%
#       group_by(OTU, Family) %>%
#       mutate(GenusDistinct = n_distinct(Genus, na.rm = remove_NAs),
#              Genus = ifelse(GenusDistinct != 1, NA, Genus)) %>%
#       ungroup() %>%
#       group_by(OTU, Genus) %>%
#       mutate(SpeciesDistinct = n_distinct(Species, na.rm = remove_NAs),
#              Species = ifelse(SpeciesDistinct != 1, NA, Species)) %>%
#       ungroup() %>%
#       group_by(OTU) %>%
#       distinct()
#   }
#   
#   # Remove pesky zeros
#   if(remove_NAs == T) {
#     tax_table <- tax_table %>%
#       mutate(across(contains("Distinct"), ~ifelse(.x == 0, NA, .x)))
#   }
#   
#   if(DistinctTally == F) {
#     tax_table <- tax_table %>%
#       select(!contains("Distinct"))
#   }
#   
#   return(tax_table)
# }
# 
# 
# test <- calc_taxonomic_conflicts(singlem_otu_table = singlem_otu_S3.11, remove_NAs = TRUE,
#                                  nocollapse = T)
# 
# test1 <- test %>%
#   select(OTU:Species, sequence, contains("Distinct")) %>%
#   ungroup() %>% group_by(OTU) %>% distinct() %>%
#   mutate(DistinctTally = sum(c_across(contains("Distinct")), na.rm = TRUE)) %>%
#   ungroup() %>%
#   group_by(OTU) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(rows = map_dbl(data, nrow))
# 
# DistinctTally.S3.11 <- test %>%
#   ungroup() %>%
#   select(OTU, ends_with("Distinct"),sequence) %>%
#   distinct() %>%
#   group_by(OTU) %>%
#   add_tally(name = "NOTU") %>% 
#   ungroup() %>%
#   unite("TaxStrLength", contains("Distinct"),sep = "",na.rm = TRUE, remove = FALSE) %>% 
#   group_by(OTU) %>%
#   mutate(maxstrlength = max(str_length(TaxStrLength))) %>%
#   ungroup() %>%
#   mutate(TaxStrIsLongest = ifelse(maxstrlength == str_length(TaxStrLength), T, F)) %>%
#   mutate(keep = ifelse(NOTU == 1, "yes",
#                        ifelse(TaxStrIsLongest == TRUE, "yes", "no"))) %>%
#   filter(keep == "yes") %>%
#   mutate(OTUDistinct = 1)
# 
# Sum <- colSums(DistinctTally.S3.11[,3:9]>1, na.rm = T)
# ncount <- colSums(!is.na(DistinctTally.S3.11[,3:9]))
# 
# (Sum)/ncount
# 
# DistinctTally.S3.11_seq <- DistinctTally.S3.11 %>%
#   select(OTU, contains("Distinct"), sequence)
# 
# 
# test <- calc_taxonomic_conflicts(singlem_otu_table = singlem_otu_S3.11, remove_NAs = FALSE)
# 
# DistinctTally.S3.11_na_kept <- test %>%
#   ungroup() %>%
#   select(OTU, ends_with("Distinct")) %>%
#   distinct() %>%
#   mutate(OTUDistinct = 1) %>%
#   summarize(across(ends_with("Distinct"), ~sum(.))) %>%
#   mutate(across(ends_with("Distinct"), ~.x-OTUDistinct, .names = "NumNonDist_{.col}")) %>%
#   mutate(across(!contains("NumNon"), ~100*(.x-OTUDistinct)/OTUDistinct, 
#                 .names = "PercNonDist_{.col}")) %>%
#   pivot_longer(cols = contains("_")) %>%
#   separate(name, into = c("NonDist", "Dist")) %>%
#   select(!ends_with("Distinct")) %>%
#   pivot_wider(names_from = "Dist", values_from = "value")  
# Sum <- colSums(DistinctTally.S3.11_na_kept[,2:8]>1, na.rm = T)
# ncount <- colSums(!is.na(DistinctTally.S3.11_na_kept[,2:8]))
# 
# 100*(Sum)/ncount
# 
# 
# test <- calc_taxonomic_conflicts(singlem_otu_table = singlem_otu_S3.5, remove_NAs = TRUE)
# 
# DistinctTally.S3.5 <- test %>%
#   ungroup() %>%
#   select(OTU, ends_with("Distinct"), sequence) %>%
#   distinct() %>%
#   group_by(OTU) %>%
#   add_tally(name = "NOTU") %>% 
#   ungroup() %>%
#   unite("TaxStrLength", contains("Distinct"),sep = "",na.rm = TRUE, remove = FALSE) %>% 
#   group_by(OTU) %>%
#   mutate(maxstrlength = max(str_length(TaxStrLength))) %>%
#   ungroup() %>%
#   mutate(TaxStrIsLongest = ifelse(maxstrlength == str_length(TaxStrLength), T, F)) %>%
#   mutate(keep = ifelse(NOTU == 1, "yes",
#                        ifelse(TaxStrIsLongest == TRUE, "yes", "no"))) %>%
#   filter(keep == "yes") %>%
#   mutate(OTUDistinct = 1)
# 
# Sum <- colSums(DistinctTally.S3.5[,3:9]>1, na.rm = T)
# ncount <- colSums(!is.na(DistinctTally.S3.5[,3:9]))
# 
# 100*(Sum)/ncount
# 
# DistinctTally.S3.5_seq <- DistinctTally.S3.5 %>%
#   select(OTU, contains("Distinct"), sequence)
# test <- calc_taxonomic_conflicts(singlem_otu_table = singlem_otu_S3.5, remove_NAs = TRUE)
# 
# DistinctTally.S3.5 <- test %>%
#   ungroup() %>%
#   select(OTU, ends_with("Distinct")) %>%
#   distinct() %>%
#   mutate(OTUDistinct = 1) %>%
#   summarize(across(ends_with("Distinct"), ~sum(.))) %>%
#   mutate(across(ends_with("Distinct"), ~.x-OTUDistinct, .names = "NumNonDist_{.col}")) %>%
#   mutate(across(!contains("NumNon"), ~100*(.x-OTUDistinct)/OTUDistinct, 
#                 .names = "PercNonDist_{.col}")) %>%
#   pivot_longer(cols = contains("_")) %>%
#   separate(name, into = c("NonDist", "Dist")) %>%
#   select(!ends_with("Distinct")) %>%
#   pivot_wider(names_from = "Dist", values_from = "value")  


singlem_otu_S3.11

tax_table <- tax_table %>%
  group_by(sequence) %>%
  mutate(OTU = paste0(gene, "_OTU_", 1:n()))

#### =================================================================================================== ####
# Convert singlm OTU tables to wide format
#### =================================================================================================== ####
# Load the MAG abundances
convert_singlem_otu <- function(otu_table.lf = NULL) {
  taxa_table <- otu_table.lf %>%
    select(gene, sequence) %>%
    distinct() %>%
    mutate(OTU = paste0(gene, "_OTU_",1:n())) %>%
    left_join(otu_table.lf %>% select(sequence, taxonomy),
              by = "sequence") %>%
    select(OTU, gene, sequence, taxonomy)
  
  
  d <- otu_table.lf %>%
    select(temporal_sample_id, gene, sequence, taxonomy, coverage, num_hits) %>%
    left_join(taxa_table %>% 
                select(OTU, sequence) %>%
                distinct(), by = c("sequence")) %>%
    #pivot_longer(cols = !OTU, names_to = "temporal_sample_id", values_to = "coverage") %>%
    filter(coverage > 0) %>%
    group_by(temporal_sample_id) %>%
    mutate(relabund_of_recovered = coverage / sum(coverage)) %>%
    select(temporal_sample_id, OTU, coverage, relabund_of_recovered)
  
  rel_abund <- d %>%
    select(OTU, temporal_sample_id, relabund_of_recovered) %>%
    pivot_wider(names_from = "temporal_sample_id", values_from = "relabund_of_recovered",
                values_fill = 0)
  coverage <- d %>%
    select(OTU, temporal_sample_id, coverage) %>%
    pivot_wider(names_from = temporal_sample_id, values_from = coverage,
                values_fill = 0)
  
  singlem_otu <- list(trimmed_mean = d, 
                      rel_abund = rel_abund, 
                      coverage = coverage,
                      taxa_table = taxa_table)
  
  return(singlem_otu)
}

singlem_otu_S3.11.wf <- convert_singlem_otu(otu_table.lf = singlem_otu_S3.11)

singlem_otu_S3.5.wf <- convert_singlem_otu(otu_table.lf = singlem_otu_S3.5)
#### =================================================================================================== ####

# Create fasta for tree building
#### =================================================================================================== ####
# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the taxonomy dataframe for the writeRepSetFasta function
fasta_ids_S3.11 <- singlem_otu_S3.11.wf$taxa_table %>%
  select(OTU, sequence) %>%
  distinct() %>%
  rename(seq = sequence, 
         name = OTU)

fasta_ids_S3.5 <- singlem_otu_S3.5.wf$taxa_table %>%
  select(OTU, sequence) %>%
  distinct() %>%
  rename(seq = sequence, 
         name = OTU)


# write fasta file
writeRepSetFasta(fasta_ids_S3.11, paste0(outputs.fp, "/S3.11.repset.fasta"))
writeRepSetFasta(fasta_ids_S3.5, paste0(outputs.fp, "/S3.5.repset.fasta"))
#### =================================================================================================== ####

# Save outputs
#### =================================================================================================== ####
# # Data
# # singlem otu tables
saveRDS(singlem_otu_S3.11.wf, paste0(outputs.fp, "/singlem_otu_S3.11.wf.RDS"))
write.csv(singlem_otu_S3.11.wf$rel_abund,
          paste0(outputs.fp, "/singlem_otu_S3.11.rel_abund.csv"),
          quote=FALSE, row.names = FALSE)
write.csv(singlem_otu_S3.11.wf$coverage,
          paste0(outputs.fp, "/singlem_otu_S3.11.count.csv"),
          quote=FALSE, row.names = FALSE)
write.csv(singlem_otu_S3.11.wf$taxa_table %>% distinct(),
          paste0(outputs.fp, "/singlem_otu_S3.11.taxonomy.csv"),
          quote=FALSE, row.names = FALSE)

saveRDS(singlem_otu_S3.5.wf, paste0(outputs.fp, "/singlem_otu_S3.5.wf.RDS"))
write.csv(singlem_otu_S3.5.wf$rel_abund,
          paste0(outputs.fp, "/singlem_otu_S3.5.rel_abund.csv"),
          quote=FALSE, row.names = FALSE)
write.csv(singlem_otu_S3.5.wf$coverage,
          paste0(outputs.fp, "/singlem_otu_S3.5.count.csv"),
          quote=FALSE, row.names = FALSE)
write.csv(singlem_otu_S3.5.wf$taxa_table %>% distinct(),
          paste0(outputs.fp, "/singlem_otu_S3.5.taxonomy.csv"),
          quote=FALSE, row.names = FALSE)
