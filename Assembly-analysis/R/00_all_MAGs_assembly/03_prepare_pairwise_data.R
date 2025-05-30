#' ## Prepare data for pairwise modeling (make betanull.lf.diff)
#' This step takes the betanull.lf and creates columns showing pairwise differences
#' in abiotic and biotic factors

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
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(GGally)
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
#### ====================================================================== ####
#### Read in data 
#### ====================================================================== ####
# Read in OTU table and env data
input <- input_ra
input$otu_table <- input$otu_table[-1]

# Turn Habitat into a factor
input$sample_metadata <- input$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, 
                            levels = c("Palsa",  "Bog", "Fen")))
# Set assembly levels
assembly_levels <- c("Homogenous selection", "Heterogenous selection",
                     "Homogenizing dispersal", "Dispersal limitation and drift",
                     "Drift")
# Read in long format
betanull.lf <- read_csv(file = paste0(outputs.fp, "/betanull.lf.csv"))

# Read in Methanogens
# methanogens <- read_tsv(file = here("data", "Methanogens.txt"), 
#                         col_names = "methanogens") %>%
#   filter(methanogens %in% input_ra$taxonomy$genome)

# manual methanogen calls prepared in setup.R
methanogens <- manual_methanogen_calls %>%
  filter(call) %>%
  filter(genome %in% input_ra$taxonomy$genome) %>%
  select(genome) %>%
  distinct() %>%
  rename(methanogens = genome)

#### ====================================================================== ####

# Prepare and clean data
#### ====================================================================== ####

# pull out taxa of interest and add their differential abundances to betnull.lf
MAGsOfInterest <- c(c("PLGY01", # m. stor
                      "PMEE01", # m.crilli
                      "20170700_S25_26", # A.stor
                      "20100900_E1D_10", # A. sp 2
                      "20120700_S1X_2", # Acidimicrobiales
                      "20120800_E3X_27", # terracidiphilus
                      "20120800_E1D_14", # terracidiphilus #2
                      "20150700_E34_11", # desulfobacterota sp
                      "20120800_S3X_7"), # Granulicella  
                    methanogens$methanogens)

# Remove mags that aren't seen very commonly: 20130700_E3Z_4; 3300037089_14;3300037104_21; 201100800_E3D_10; 201100800_E3D_36; 201200800_E2X_5
MAGsOfInterest <- setdiff(MAGsOfInterest, c("20130700_E3Z_4", 
                                            "3300037089_14",
                                            "3300037104_21", 
                                            "201100800_E3D_10", 
                                            "201100800_E3D_36", 
                                            "201200800_E2X_5"))

#MAGsOfInterest <- input_ra$taxonomy$genome

MAGsOfInterest <- c(vip_members$MAG, MAGsOfInterest)

MAGsOfInterest <- unique(MAGsOfInterest)

MAG_diff_Abund.ls <- lapply(MAGsOfInterest, function(x) {
  otu_diff <- as.data.frame(t(input_ra$otu_table[-1])) %>%
    rownames_to_column(var = "SampleID") %>%
    select(SampleID, all_of(x)) %>%
    column_to_rownames(var = "SampleID") %>% dist() %>%
    as.matrix() %>% as.data.frame() %>% 
    rownames_to_column(var = "Site1") %>%
    pivot_longer(-Site1, names_to = "Site2", values_to = paste0(x,"_mag_diff_abund")) %>%
    filter(Site1 != Site2)
  return(otu_diff)
  
})

MAG_diff_Abund <- reduce(MAG_diff_Abund.ls, left_join, by = c("Site1", "Site2"))
#### ====================================================================== ####

# Prepare Reaction data (code from Sam)
#### ====================================================================== ####
# path_abund <- trimmed_mean$rel_abund %>%
#   select(genome, all_of(sample_metadata$temporal_sample_id)) %>%
#   pivot_longer(-genome, names_to = "sample", values_to = "rel_abund") %>%
#   left_join(emerge_product_refined, by = "genome")

# cumu_path_abund <- path_abund %>%
#   filter(call) %>%
#   group_by(subpathway, sample) %>%
#   summarise(abund = sum(rel_abund)) %>%
#   pivot_wider(names_from = "subpathway", values_from = "abund")

# pull out pathways of interest and add their differential abundances to betnull.lf
if(exists("pathway_abundance_cumu")) {
  print("TRUE")
  PathsOfInterest <- names(pathway_abundance_cumu[-1])
  
  Path_diff_Abund.ls <- lapply(PathsOfInterest, function(x) {
    path_diff <- pathway_abundance_cumu %>%
      select(sample, all_of(x)) %>%
      column_to_rownames(var = "sample") %>% dist() %>%
      as.matrix() %>% as.data.frame() %>% 
      rownames_to_column(var = "Site1") %>%
      pivot_longer(-Site1, names_to = "Site2", values_to = paste0(x,"_path_diff_abund")) %>%
      filter(Site1 != Site2)
    return(path_diff)
    
  })
  
  Path_diff_Abund <- reduce(Path_diff_Abund.ls, left_join, by = c("Site1", "Site2"))
} else { # if the pathways don't exist, just use columns site1 and site2
  Path_diff_Abund <- betanull.lf %>%
    select(Site1, Site2)
}


#### ====================================================================== ####
# Create pairwise calculation of alpha diversity
# based off of Ben's suggestion to see if assembly processes show relationship
# with alpha diversity
#### ====================================================================== ####
alpha_div <- data.frame(Shannon = vegan::diversity(input$otu_table, index = "shannon", MARGIN = 2), 
                        Richness = apply(input$otu_table, 2, function(x) length(x[x > 0])),
                        SampleID = names(input$otu_table))
avg_alpha_div <- betanull.lf %>% 
  select(Site1, Site2) %>%
  left_join(alpha_div, by = c("Site1" = "SampleID")) %>%
  rename(Shannon.Site1 = Shannon, Richness.Site1 = Richness) %>%
  left_join(alpha_div, by = c("Site2" = "SampleID")) %>%
  rename(Shannon.Site2 = Shannon, Richness.Site2 = Richness) %>%
  mutate(AvgShannon = (Shannon.Site1 + Shannon.Site2)/2,
         AvgRichness = (Richness.Site1 + Richness.Site2)/2) %>%
  select(Site1, Site2, AvgShannon, AvgRichness)

#### ====================================================================== ####

# Create a version of betanull.lf that has pairwise differences for various categories
# Including dram annotations
#### ====================================================================== ####
betanull.lf.diff <- betanull.lf %>%
  # make a selectiondispersaldrift category
  mutate(SelectionDispersalDrift = ifelse(grepl("selection", Assembly_Process),
                                          "Selection",
                                          ifelse(grepl("ispersal", Assembly_Process),
                                                 "Dispersal", "Drift"))) %>%
  # calculate pairwise variables
  mutate(Year_diff = abs(Year__.Site1 - Year__.Site2),
         Depth_diff = abs(DepthAvg__.Site1 - DepthAvg__.Site2), 
         habcomb = paste0(Habitat__.Site1,"_", Habitat__.Site2),
         Habitat_comp = ifelse(grepl("Bog_Palsa", habcomb), "Palsa_Bog", 
                               ifelse(grepl("Fen_Palsa", habcomb), "Palsa_Fen",
                                      ifelse(grepl("Fen_Bog", habcomb), "Bog_Fen", habcomb))),
         pH_porewater_diff = abs(pH_porewater.Site1 - pH_porewater.Site2),
         alphaC_diff = abs(alphaC.Site1 - alphaC.Site2),
         T_air_diff = abs(T_air.deg_C.Site1 - T_air.deg_C.Site2),
         T_soil_diff = abs(T_soil.deg_C.Site1 - T_soil.deg_C.Site2),
         T_air_7d_diff = abs(mean_AirTemperature_7d.Site1 - mean_AirTemperature_7d.Site2),
         T_air_14d_diff = abs(mean_AirTemperature_14d.Site1 - mean_AirTemperature_14d.Site2),
         T_air_28d_diff = abs(mean_AirTemperature_28d.Site1 - mean_AirTemperature_28d.Site2),
         T_air_7d_sd_diff = abs(sd_AirTemperature_7d.Site1 - sd_AirTemperature_7d.Site2),
         T_air_14d_sd_diff = abs(sd_AirTemperature_14d.Site1 - sd_AirTemperature_14d.Site2),
         T_air_28d_sd_diff = abs(sd_AirTemperature_28d.Site1 - sd_AirTemperature_28d.Site2),
         pct_time_below_WTD_21d_diff = abs(pct_time_below_WTD_21d.Site1 - pct_time_below_WTD_21d.Site2),
         pct_time_below_WTD_28d_diff = abs(pct_time_below_WTD_28d.Site1 - pct_time_below_WTD_28d.Site2),
         Samp_dist_WT_diff = abs(Samp_dist_WT.Site1 - Samp_dist_WT.Site2),
         Min_PF_dist_from_sample_diff = abs(Min_PF_dist_from_sample.Site1 - Min_PF_dist_from_sample.Site2)) %>%
  # make within/between variables
  mutate(Year_wbtn = ifelse(Year_diff == 0, "within", "between"),
         Depth_wbtn = ifelse(DepthLumping.Site1 == DepthLumping.Site2, "within", "between"),
         Habitat_wbtn = ifelse(Habitat__.Site1 == Habitat__.Site2, "within", "between")) %>%
  # Add in taxa features of interest
  left_join(MAG_diff_Abund, by = c("Site1", "Site2")) %>%
  # Add in the Cumulative pathway abundance differences
  left_join(Path_diff_Abund, by = c("Site1", "Site2")) %>%
  # Add avg alpha diversity
  left_join(avg_alpha_div, by = c("Site1", "Site2")) %>%
  # Scale numeric predictors
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE)),
                .names = "scl_{.col}")) %>%
  # Reformat what needs to be reformatted
  mutate(across(where(is.character), ~factor(.x))) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels)) %>%
  select(contains("Assembly"), StochasticDeterministic, BetaNTI, RCBC.nona, RCBC, contains("diff"), contains("comp"),
         contains("wbtn"), contains("Shannon"), contains("Richness"), Site1, Site2) %>%
  filter(Site1 != "714_P1_30-34") %>%
  filter(Site2 != "714_P1_30-34")


names(betanull.lf.diff)

# how many of each type of assembly process do we have?
# betanull.lf.diff %>%
#   group_by(Assembly_Process) %>%
#   tally() %>%
#   mutate(Percent = 100* n/sum(n)) %>% # not horrendous
#   ggplot(aes(y = Percent, fill = Assembly_Process, x = 1,
#              group = Assembly_Process)) +
#   barchart_plotting_layers

betanull.lf.diff %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Percent = 100* n/sum(n)) # Very uneven

betanull.lf.diff %>%
  group_by(StochasticDeterministic) %>%
  tally() %>%
  mutate(Percent = 100* n/sum(n))

# What do the numeric predictors look like scaled and unscaled?
# ggplot(data = betanull.lf.diff %>%
#          pivot_longer(cols = where(is.numeric), names_to = "Variables", 
#                       values_to = "Values") %>%
#          mutate(scaled = ifelse(grepl("scl", Variables), "scaled", "not_scaled")),
#        aes(x = Values)) +
#   geom_histogram() +
#   facet_wrap(scaled~Variables, scale = "free_x", nrow = 2)

#### ====================================================================== ####

# Save outputs
#### ====================================================================== ####
# # Data
# assembly results in long format
saveRDS(betanull.lf.diff, paste0(outputs.fp, "/betanull.lf.diff.RDS"))
write.csv(betanull.lf.diff,
          paste0(outputs.fp, "/betanull.lf.diff.csv"),
          quote=TRUE, row.names = FALSE)
