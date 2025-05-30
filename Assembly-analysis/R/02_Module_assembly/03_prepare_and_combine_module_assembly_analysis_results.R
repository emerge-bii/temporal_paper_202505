#' ## Prepare and combine Assembly Analysis results
#' This step takes the output from 01_calc_betaNTI.R and 02_calc_rcbc.R and runs
#' initial figures and analyses on them.

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
#### Helpful Functions
#### ====================================================================== ####
read_module_input_data <- function(module_input_data_folder = input_data_folder,
                                   module_assembly_folder = mod_assemb_fold,
                                   module_habitat = "bog") {
  # Read input data
  input_files <- list.files(module_input_data_folder, pattern = ".RDS")
  mod_info_list_names <- gsub(paste0(module_habitat, "_module_"),"", input_files)
  mod_info_list_names <- gsub("\\.RDS","", mod_info_list_names)
  mod_info_list <- lapply(input_files, function(x) {
    readRDS(file = paste0(module_input_data_folder, x))
  })
  names(mod_info_list) <- mod_info_list_names
  
  # Read betaNTI data
  habitat_module <- paste0(module_habitat, "_", mod_info_list_names)
  mod_info_list <- lapply(mod_info_list_names, function(x) {
    bNTI <- read_csv(file = paste0(module_assembly_folder, "/", 
                                                      module_habitat, "_", x, "_weighted_bNTI.csv")) %>%
      rename(SampleID = 1) %>%
      column_to_rownames(var = "SampleID")
    
    mod_info_list[[which(mod_info_list_names == x)]][["bNTI"]] <- bNTI
    return(mod_info_list[[which(mod_info_list_names == x)]])
  })
  names(mod_info_list) <- mod_info_list_names
  
  # Read RCBC data
  habitat_module <- paste0(module_habitat, "_", mod_info_list_names)
  mod_info_list <- lapply(mod_info_list_names, function(x) {
    RCBC_file <- paste0(module_assembly_folder, "/", 
                        module_habitat, "_", x, "_rcbc_matrix.csv")
    if(file.exists(RCBC_file)) {
      RCBC.raw <- read_csv(file = RCBC_file)
      
      nsamples <- ncol(mod_info_list[[which(mod_info_list_names == x)]][["otu_table"]]) - 1
      
      RCBC <- RCBC.raw[c(as.character(1:nsamples)),
                       c(as.character(1:nsamples))]
      RCBC <- as.matrix(RCBC)
      colnames(RCBC) <- names(mod_info_list[[which(mod_info_list_names == x)]][["otu_table"]])[2:(nsamples+1)] # first column is "genome"
      rownames(RCBC) <- names(mod_info_list[[which(mod_info_list_names == x)]][["otu_table"]])[2:(nsamples+1)]
      
      RCBC[upper.tri(RCBC)] <- 999 
      
      mod_info_list[[which(mod_info_list_names == x)]][["RCBC"]] <- RCBC
      return(mod_info_list[[which(mod_info_list_names == x)]])
    } else {
      return(mod_info_list[[which(mod_info_list_names == x)]])
    }
    
  })
  names(mod_info_list) <- mod_info_list_names
  
  return(mod_info_list)
}


# Convert to long format
convert_to_lf <- function(module_info_list = input_info_list, 
                          module = "red") {
  # Get BNTI and RCBC wide format
  betaNTI <- module_info_list[[module]][["bNTI"]]
  RCBC <- module_info_list[[module]][["RCBC"]]
  sample_metadata <- module_info_list[[module]][["sample_metadata"]]
  
  
  # Convert BetaNTI to long format
  bnti.lf <- betaNTI %>%
    as.data.frame() %>%
    rownames_to_column(var = "Site1") %>%
    pivot_longer(cols = !matches("Site1"), names_to = "Site2",
                 values_to = "BetaNTI") %>% 
    # remove duplicated comparisons
    filter(!is.na(BetaNTI)) %>% 
    # Add column of BetaNTI interpretation
    mutate(Assembly_Process = ifelse(BetaNTI < -2, "Homogenous selection",
                                     ifelse(BetaNTI > 2, "Heterogenous selection", 
                                            "Stochastic"))) %>%
    # Join metadata for site 1 and 2
    left_join(sample_metadata, by = c("Site1" = "temporal_sample_id"))%>%
    left_join(sample_metadata, by = c("Site2" = "temporal_sample_id"),
              suffix = c(".Site1", ".Site2")) %>%
    mutate(Module_name = module)
  
  # Convert RCBC to long format
  RCBC.lf <- RCBC %>%
    as.data.frame() %>%
    rownames_to_column(var = "Site1") %>% 
    pivot_longer(cols = !Site1, names_to = "Site2",
                 values_to = "RCBC") %>%
    # remove duplicated comparisons
    filter(Site1 != Site2) %>% # self-comparisons
    filter(RCBC != 999) %>% # duplicate comparisons from upper triangle
    # Add column of BetaNTI interpretation
    mutate(Assembly_Process = ifelse(RCBC < -0.95, "Homogenizing dispersal",
                                     ifelse(RCBC <= 0.95, "Drift",
                                            ifelse(RCBC > 0.95, "Dispersal limitation and drift",
                                                   NA)))) %>%
    mutate(Module_name = module)

  
  # Join the two formats
  betanull.lf <- full_join(bnti.lf %>%
                             mutate(Assembly_Process = ifelse(Assembly_Process == 
                                                                "Stochastic",
                                                              NA, Assembly_Process)),
                           RCBC.lf, by = c("Site1", "Site2", "Module_name")) %>% 
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
    
    select(Site1, Site2, BetaNTI, RCBC, Assembly_Process.BetaNTI, 
           Assembly_Process.RCBC, Assembly_Process, everything())
  
  
  return(betanull.lf)
  
}

# Function for combining module assembly stats
comb_mod_assembly_data <- function(input_info_list = all_info_list,
                                   habitat = "all") {
  input_mod_names <- names(input_info_list)
  
  # Convert the assemlby results to long format
  input_info_list <- lapply(names(input_info_list), function(x) {
    print(x)
    betanull.lf.tmp <- convert_to_lf(module_info_list = input_info_list,
                                     module = x)
    input_info_list[[x]][["betanull.lf"]] <- betanull.lf.tmp
    
    return(input_info_list[[x]])
  } )
  
  names(input_info_list) <- input_mod_names # rename the list
  
  # Extract one betanull.lf for all modules combined
  input_mod_betanull <- lapply(seq_along(input_info_list), function(x) {
    return(input_info_list[[x]][["betanull.lf"]])
  })
  names(input_mod_betanull) <- input_mod_names # rename
  
  betanull.lf_mod <- reduce(input_mod_betanull, bind_rows)
  
  # combine with original assembly analysis for comparison
  if(habitat == "all") {
    betanull.lf_mod <- bind_rows(betanull.lf_mod, 
                                 betanull.lf %>%
                                   mutate(Module_name = "all") %>%
                                   filter(Habitat__.Site1 == Habitat__.Site2)) %>% # removes cross-habitat compairons for "all"
      filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
      # Remove the grey module
      filter(Module_name != "grey")
  } else {
    betanull.lf_mod <- bind_rows(betanull.lf_mod, 
                                 betanull.lf %>%
                                   mutate(Module_name = "all") %>%
                                   filter(Habitat__.Site1 == habitat & Habitat__.Site2 == habitat)) %>%
      filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
      # Remove the grey module
      filter(Module_name != "grey")
  }
   
  # Combine two sets of results
  return(list(betanull.lf_mod = betanull.lf_mod,
              input_info_list = input_info_list))
}







#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# set up input folders
all_data_folder <- paste0(outputs.fp, "/all_module_otus/")
palsa_data_folder <- paste0(outputs.fp, "/palsa_module_otus/")
bog_data_folder <- paste0(outputs.fp, "/bog_module_otus/")
fen_data_folder <- paste0(outputs.fp, "/fen_module_otus/")


mod_assemb_fold <- paste0(outputs.fp, "/module_assembly/")

# Read in each habitat module info
all_info_list <- read_module_input_data(module_input_data_folder = all_data_folder,
                       module_assembly_folder = mod_assemb_fold,
                       module_habitat = "all")


palsa_info_list <- read_module_input_data(module_input_data_folder = palsa_data_folder,
                                          module_assembly_folder = mod_assemb_fold,
                                          module_habitat = "palsa")

bog_info_list <- read_module_input_data(module_input_data_folder = bog_data_folder,
                                          module_assembly_folder = mod_assemb_fold,
                                          module_habitat = "bog")

fen_info_list <- read_module_input_data(module_input_data_folder = fen_data_folder,
                                          module_assembly_folder = mod_assemb_fold,
                                          module_habitat = "fen")



# Read in full betaNTI results for comparison
betaNTI <- read_csv(file = paste0(outputs.fp, "/weighted_bNTI.csv")) %>%
  rename(SampleID = 1) %>%
  column_to_rownames(var = "SampleID")
betaNTI.all <- betaNTI %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site1") %>%
  pivot_longer(cols = !matches("Site1"), names_to = "Site2",
               values_to = "BetaNTI") %>% 
  # remove duplicated comparisons
  filter(!is.na(BetaNTI)) %>% 
  # Add column of BetaNTI interpretation
  mutate(Assembly_Process = ifelse(BetaNTI < -2, "Homogenous selection",
                                   ifelse(BetaNTI > 2, "Heterogenous selection", 
                                          "Stochastic"))) %>%
  # Join metadata for site 1 and 2
  left_join(input_ra$sample_metadata, by = c("Site1" = "temporal_sample_id"))%>%
  left_join(input_ra$sample_metadata, by = c("Site2" = "temporal_sample_id"),
            suffix = c(".Site1", ".Site2")) %>%
  mutate(Module_name = "all")

betanull.lf <- read_csv(file = paste0(outputs.fp, "/betanull.lf.csv"))

#### ====================================================================== ####

# Convert to long format and prepare for plotting
#### ====================================================================== ####
all_comb <- comb_mod_assembly_data(input_info_list = all_info_list,
                                     habitat = "all")
all_mod_betanull.lf <- all_comb[[1]]
all_info_list <- all_comb[[2]]

palsa_comb <- comb_mod_assembly_data(input_info_list = palsa_info_list,
                                     habitat = "Palsa")
palsa_mod_betanull.lf <- palsa_comb[[1]]
palsa_info_list <- palsa_comb[[2]]

bog_comb <- comb_mod_assembly_data(input_info_list = bog_info_list,
                                     habitat = "Bog")
bog_mod_betanull.lf <- bog_comb[[1]]
bog_info_list <- bog_comb[[2]]

fen_comb <- comb_mod_assembly_data(input_info_list = fen_info_list,
                       habitat = "Fen")
fen_mod_betanull.lf <- fen_comb[[1]]
fen_info_list <- fen_comb[[2]]



#### ====================================================================== ####

# Plot All
#### ====================================================================== ####
all_colors <- unique(all_mod_betanull.lf$Module_name)
names(all_colors) <- unique(all_mod_betanull.lf$Module_name)
all_colors["all"] <- "white"

intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              DistanceMeasure = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
all_betanull.lf_mod.plot <- all_mod_betanull.lf %>%
  pivot_longer(all_of(c("BetaNTI", "RCBC")), 
               names_to = "DistanceMeasure", values_to = "Distance") %>% 
  select(Distance, DistanceMeasure, Assembly_Process, Site1, Site2, Module_name, everything())

all_dot_plot <- ggplot(all_betanull.lf_mod.plot, 
                         aes(y = Distance, x = Module_name)) +
  geom_point(position = "jitter", alpha = 0.5) + 
  geom_violin(aes(color = Module_name, fill = Module_name), alpha = 0.3, weight = 2) +
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
  geom_hline(data = intercept_guide,
             aes(yintercept = Distance),
             linetype = "dashed") +
  facet_wrap(~DistanceMeasure, scales = "free_y")
all_dot_plot

# Proportions
all_mod_assembly_plot <- all_mod_betanull.lf %>%
  group_by(Module_name, Assembly_Process) %>%
  # filter out cross-depth comparisons
  tally() %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = c("Homogenous selection",
                                                                "Heterogenous selection",
                                                                "Homogenizing dispersal",
                                                                "Dispersal limitation and drift",
                                                                "Drift"))) %>%
  ggplot(aes(x = Module_name, y = Percent)) +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", position = "dodge") +
  scale_fill_viridis(name = "Assembly Process", discrete = TRUE) +
  theme_bw() +
  geom_vline(xintercept = 1.5)
all_mod_assembly_plot

#### ====================================================================== ####

# Plot Palsa
#### ====================================================================== ####
palsa_colors <- unique(palsa_mod_betanull.lf$Module_name)
names(palsa_colors) <- unique(palsa_mod_betanull.lf$Module_name)
palsa_colors["all"] <- "white"

intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              DistanceMeasure = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
palsa_betanull.lf_mod.plot <- palsa_mod_betanull.lf %>%
  pivot_longer(all_of(c("BetaNTI", "RCBC")), 
               names_to = "DistanceMeasure", values_to = "Distance") %>% 
  select(Distance, DistanceMeasure, Assembly_Process, Site1, Site2, Module_name, everything())

palsa_dot_plot <- ggplot(palsa_betanull.lf_mod.plot, 
                         aes(y = Distance, x = Module_name)) +
  geom_point(position = "jitter", alpha = 0.5) + 
  geom_violin(aes(color = Module_name, fill = Module_name), alpha = 0.3, weight = 2) +
  scale_color_manual(values = palsa_colors) +
  scale_fill_manual(values = palsa_colors) +
  geom_hline(data = intercept_guide,
             aes(yintercept = Distance),
             linetype = "dashed") +
  facet_wrap(~DistanceMeasure, scales = "free_y")
palsa_dot_plot

# Proportions
palsa_mod_assembly_plot <- palsa_mod_betanull.lf %>%
  group_by(Module_name, Assembly_Process) %>%
  # filter out cross-depth comparisons
  tally() %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = c("Homogenous selection",
                                                                "Heterogenous selection",
                                                                "Homogenizing dispersal",
                                                                "Dispersal limitation and drift",
                                                                "Drift"))) %>%
  ggplot(aes(x = Module_name, y = Percent)) +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", position = "dodge") +
  scale_fill_viridis(name = "Assembly Process", discrete = TRUE) +
  theme_bw() +
  geom_vline(xintercept = 1.5)
palsa_mod_assembly_plot

#### ====================================================================== ####

# Plot Bog
#### ====================================================================== ####
bog_colors <- unique(bog_mod_betanull.lf$Module_name)
names(bog_colors) <- unique(bog_mod_betanull.lf$Module_name)
bog_colors["all"] <- "white"

intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              DistanceMeasure = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
bog_betanull.lf_mod.plot <- bog_mod_betanull.lf %>%
  pivot_longer(all_of(c("BetaNTI", "RCBC")), 
               names_to = "DistanceMeasure", values_to = "Distance") %>% 
  select(Distance, DistanceMeasure, Assembly_Process, Site1, Site2, Module_name, everything())

bog_dot_plot <- ggplot(bog_betanull.lf_mod.plot, 
                         aes(y = Distance, x = Module_name)) +
  geom_point(position = "jitter", alpha = 0.5) + 
  geom_violin(aes(color = Module_name, fill = Module_name), alpha = 0.3, weight = 2) +
  scale_color_manual(values = bog_colors) +
  scale_fill_manual(values = bog_colors) +
  geom_hline(data = intercept_guide,
             aes(yintercept = Distance),
             linetype = "dashed") +
  facet_wrap(~DistanceMeasure, scales = "free_y")
bog_dot_plot

# Proportions
bog_mod_assembly_plot <- bog_mod_betanull.lf %>%
  group_by(Module_name, Assembly_Process) %>%
  # filter out cross-depth comparisons
  tally() %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = c("Homogenous selection",
                                                                "Heterogenous selection",
                                                                "Homogenizing dispersal",
                                                                "Dispersal limitation and drift",
                                                                "Drift"))) %>%
  ggplot(aes(x = Module_name, y = Percent)) +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", position = "dodge") +
  scale_fill_viridis(name = "Assembly Process", discrete = TRUE) +
  theme_bw() +
  geom_vline(xintercept = 1.5)
bog_mod_assembly_plot

#### ====================================================================== ####

# Plot Fen
#### ====================================================================== ####
fen_colors <- unique(fen_mod_betanull.lf$Module_name)
names(fen_colors) <- unique(fen_mod_betanull.lf$Module_name)
fen_colors["all"] <- "white"

intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              DistanceMeasure = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
fen_betanull.lf_mod.plot <- fen_mod_betanull.lf %>%
  pivot_longer(all_of(c("BetaNTI", "RCBC")), 
               names_to = "DistanceMeasure", values_to = "Distance") %>% 
  select(Distance, DistanceMeasure, Assembly_Process, Site1, Site2, Module_name, everything())

fen_dot_plot <- ggplot(fen_betanull.lf_mod.plot, 
                         aes(y = Distance, x = Module_name)) +
  geom_point(position = "jitter", alpha = 0.5) + 
  geom_violin(aes(color = Module_name, fill = Module_name), alpha = 0.3, weight = 2) +
  scale_color_manual(values = fen_colors) +
  scale_fill_manual(values = fen_colors) +
  geom_hline(data = intercept_guide,
             aes(yintercept = Distance),
             linetype = "dashed") +
  facet_wrap(~DistanceMeasure, scales = "free_y")
fen_dot_plot

# Proportions
fen_mod_assembly_plot <- fen_mod_betanull.lf %>%
  group_by(Module_name, Assembly_Process) %>%
  # filter out cross-depth comparisons
  tally() %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = c("Homogenous selection",
                                                                "Heterogenous selection",
                                                                "Homogenizing dispersal",
                                                                "Dispersal limitation and drift",
                                                                "Drift"))) %>%
  ggplot(aes(x = Module_name, y = Percent)) +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", position = "dodge") +
  scale_fill_viridis(name = "Assembly Process", discrete = TRUE) +
  theme_bw() +
  geom_vline(xintercept = 1.5)
fen_mod_assembly_plot

#### ====================================================================== ####

# Save outputs
#### ====================================================================== ####
# # Data
# assembly results in long format
write.csv(all_mod_betanull.lf,
          paste0(outputs.fp, "/all_mod_betanull.lf.csv"),
          quote=TRUE, row.names = FALSE)
ggsave(all_mod_assembly_plot, device = "png",
       filename = paste0(figures.fp, "/all_module_assembly_prop.png"))

write.csv(palsa_mod_betanull.lf,
          paste0(outputs.fp, "/palsa_mod_betanull.lf.csv"),
          quote=TRUE, row.names = FALSE)
ggsave(palsa_mod_assembly_plot, device = "png",
       filename = paste0(figures.fp, "/palsa_module_assembly_prop.png"))

write.csv(bog_mod_betanull.lf,
          paste0(outputs.fp, "/bog_mod_betanull.lf.csv"),
          quote=TRUE, row.names = FALSE)
ggsave(bog_mod_assembly_plot, device = "png",
       filename = paste0(figures.fp, "/bog_module_assembly_prop.png"))

write.csv(fen_mod_betanull.lf,
          paste0(outputs.fp, "/fen_mod_betanull.lf.csv"),
          quote=TRUE, row.names = FALSE)
ggsave(fen_mod_assembly_plot, device = "png",
       filename = paste0(figures.fp, "/fen_module_assembly_prop.png"))
