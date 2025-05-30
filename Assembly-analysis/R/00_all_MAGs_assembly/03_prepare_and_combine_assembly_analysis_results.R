#' ## Prepare and combine Assembly Analysis results
#' This step takes the output from 01_calc_betaNTI.R and 02_calc_rcbc.R and runs
#' transforms them into useable formats for downstream analysis

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
                            levels = c("Palsa", "Collapsed Palsa", "Bog", "Fen")))

# Read in betaNTI results
betaNTI <- read_csv(file = paste0(outputs.fp, "/weighted_bNTI.csv")) %>%
  rename(SampleID = 1) %>%
  column_to_rownames(var = "SampleID")

# Create a long-form dummy data frame to combine results
samp.df <- data.frame(Site1 = colnames(betaNTI), Site2 = rownames(betaNTI))


betaNTI[upper.tri(betaNTI)] <- 999 

# Read in RCBC results
RCBC.raw <- read_csv(file = paste0(outputs.fp, "/rcbc_matrix.csv"))

RCBC <- RCBC.raw[c(as.character(1:ncol(input$otu_table))),
                 c(as.character(1:ncol(input$otu_table)))]
RCBC <- as.matrix(RCBC)
colnames(RCBC) <- names(input$otu_table)
rownames(RCBC) <- names(input$otu_table)

# Set the upper corner of RCBC to 999 to distinguish self-comparison NAs from 
# NAs generated during RCBC calculation; 
# This is important if you set speedup = TRUE because there will be NAs in the 
# matrix where no comparison was performed
RCBC[upper.tri(RCBC)] <- 999 

# Read in long format of rcbc
RCBC.lf <- read_csv(file = paste0(outputs.fp, "/rcbc_long_form.csv"))

#### ====================================================================== ####

# Transform betaNTI and RCBC results to long format
#### ====================================================================== ####
betaNTI.lf <- betaNTI %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site1") %>% 
  pivot_longer(cols = !matches("Site1"), names_to = "Site2",
               values_to = "BetaNTI") %>%
  full_join(samp.df, by = c("Site1", "Site2")) %>%
  # remove duplicated comparisons
  filter(Site1 != Site2) %>% # self-comparisons
  filter(BetaNTI != 999) %>% # duplicate comparisons from upper triangle
  # remove samples that are not in the metadata
  filter(Site1 %in% input$sample_metadata$temporal_sample_id) %>%
  filter(Site2 %in% input$sample_metadata$temporal_sample_id) %>%
  # Add column of BetaNTI interpretation
  mutate(Assembly_Process = ifelse(BetaNTI < -2, "Homogenous selection",
                                   ifelse(BetaNTI > 2, "Heterogenous selection", 
                                          "Stochastic"))) %>%
  # Join metadata for site 1 and 2
  left_join(input$sample_metadata, by = c("Site1" = "temporal_sample_id"))%>%
  left_join(input$sample_metadata, by = c("Site2" = "temporal_sample_id"),
            suffix = c(".Site1", ".Site2"))

# Sanity check, do we have the right number of rows?
# There are 16653 combinations (without self-comparisons) of 183 samples
nrow(betaNTI.lf) == 16653 # TRUE == good (combn(183, 2))

# Get list of "stochastic" comparisons, by setting any
# non-stochastic comparisons to NA
rcbc <- as.matrix(betaNTI)
rcbc[betaNTI < -2 | betaNTI > 2] <- NA
rcbc[which(!is.na(rcbc))]
image(as.matrix(rcbc))
image(as.matrix(RCBC))
image(as.matrix(betaNTI))

# RCBC
## -2 < betaNTI < 2 and RCBC < -0.95 indicates homogenizing dispersal
## -2 < betaNTI < 2 and RCBC > 0.95 indicates dispersal limitation and drift
## -2 < betaNTI < 2 and -0.95 < RCBC < 0.95 indicates drift
RCBC.lf <- RCBC %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site1") %>% 
  pivot_longer(cols = !Site1, names_to = "Site2",
               values_to = "RCBC") %>%
  # remove duplicated comparisons
  filter(Site1 != Site2) %>% # self-comparisons
  filter(RCBC != 999) %>% # duplicate comparisons from upper triangle
  # remove samples that are not in the metadata
  filter(Site1 %in% input$sample_metadata$temporal_sample_id) %>%
  filter(Site2 %in% input$sample_metadata$temporal_sample_id) %>%
  # Add column of BetaNTI interpretation
  mutate(Assembly_Process = ifelse(RCBC < -0.95, "Homogenizing dispersal",
                                   ifelse(RCBC <= 0.95, "Drift",
                                          ifelse(RCBC > 0.95, "Dispersal limitation and drift",
                                                 NA))))

# Sanity check, do we have the right number of rows?
# There are 16653 combinations (without self-comparisons) of 183 samples
# ncol(combn(xxxsamples, 2))
nrow(RCBC.lf) == 16653
#### ====================================================================== ####

#' Create a table with betaNTI and rcbc results combined
#### ====================================================================== ####
betanull.lf <- full_join(betaNTI.lf %>%
                           mutate(Assembly_Process = ifelse(Assembly_Process == 
                                                              "Stochastic",
                                                            NA, Assembly_Process)),
                         RCBC.lf, by = c("Site1", "Site2")) %>% 
  mutate(Assembly_Process = coalesce(Assembly_Process.x, Assembly_Process.y),
         StochasticDeterministic = ifelse(grepl("selection", Assembly_Process), 
                                          "Deterministic",
                                          "Stochastic"),
         Assembly_Process = factor(Assembly_Process, 
                                   levels = c("Homogenous selection",
                                              "Heterogenous selection",
                                              "Homogenizing dispersal", 
                                              "Dispersal limitation and drift", 
                                              "Drift"))) %>%
  rename(Assembly_Process.RCBC = Assembly_Process.y,
         Assembly_Process.BetaNTI = Assembly_Process.x) %>%
  mutate(Assembly_Process.RCBC = ifelse(grepl("selection", 
                                              Assembly_Process.BetaNTI), 
                                        NA, Assembly_Process.RCBC),
         RCBC.nona = RCBC,
         RCBC = ifelse(grepl("selection", Assembly_Process.BetaNTI), 
                       NA, RCBC)) %>%
  
  dplyr::select(Site1, Site2, BetaNTI, RCBC, Assembly_Process.BetaNTI, 
         Assembly_Process.RCBC, Assembly_Process, everything())

# Sanity check, do we have the right number of rows?
# There are 16653 combinations (without self-comparisons) of 183 samples
nrow(betanull.lf) == 16653


# Transform combined dataframe to wide format
betanull.wf <- betanull.lf %>%
  select(Site1, Site2, Assembly_Process, BetaNTI, RCBC) %>%
  # recode Assembly processes
  mutate(Assembly_Process_Code = factor(Assembly_Process, 
                                        levels = c("Heterogenous selection",
                                                   "Dispersal limitation and drift",
                                                   "Drift",
                                                   "Homogenizing dispersal",
                                                   "Homogenous selection")),
         Assembly_Process_level = as.numeric(Assembly_Process_Code)) %>%
  dplyr::select(Site1, Site2, Assembly_Process_level) %>%
  pivot_wider(names_from = Site1, values_from = Assembly_Process_level) %>% 
  column_to_rownames(var = "Site2")

# Fill in upper diagonal
diag(betanull.wf) <- NA
betanull.mat <- as.matrix(betanull.wf)

#### ====================================================================== ####

# Save outputs
#### ====================================================================== ####
# # Data
# assembly results in long format
saveRDS(betanull.lf, paste0(outputs.fp, "/betanull.lf.RDS"))
write.csv(betanull.lf,
          paste0(outputs.fp, "/betanull.lf.csv"),
          quote=TRUE, row.names = FALSE)

saveRDS(betanull.wf, paste0(outputs.fp, "/betanull.wf.RDS"))
write.csv(betanull.wf,
          paste0(outputs.fp, "/betanull.wf.csv"),
          quote=FALSE, row.names = FALSE)

