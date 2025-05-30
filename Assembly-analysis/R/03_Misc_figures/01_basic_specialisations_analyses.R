#' ## Interpret Assembly Analysis results
#' This step takes the output from 
#' 03_prepare_and_combine_assembly_analysis_results.R and runs initial figures 
#' and analyses on them.

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
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(GUniFrac)
library(cowplot)
library(multcompView)
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
# Rename input
input <- input_ra

#### ====================================================================== ####

# Setup plotting parameters
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(15, 16, 17)
fill_habitat   <- colour_habitat
depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth
year_levels <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017")
year_labels <- year_levels
fill_year   <- RColorBrewer::brewer.pal(7, "OrRd")
colour_year <- fill_year

special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
special_fill <- special_colour
n_special_levels <- c("nitrogen_fixer", "nitrate_reducer", "denitrifier", "nitrite_oxidiser", NA)
n_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-5:-8]
n_special_fill <- n_special_colour
s_special_levels <- c("sulfur_oxidiser", "sulfate_reducer", NA)
s_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-3:-8]
s_special_fill <- s_special_colour

all_special_levels <- c(na.exclude(special_levels), na.exclude(n_special_levels), na.exclude(s_special_levels), NA)
all_special_labels <- gsub("_", " ", all_special_levels)
all_special_colour <- c(RColorBrewer::brewer.pal(7, "OrRd")[c(7,1,6,2,5,3,4)],RColorBrewer::brewer.pal(5, "Greens")[c(2,4,1,3)],
                        RColorBrewer::brewer.pal(9, "Purples")[c(5,8)])
all_special_fill <- all_special_colour

#### ====================================================================== ####

# Define some helper functions
#### ====================================================================== ####

#### ====================================================================== ####

# Specialisations
#### ====================================================================== ####
# remove nitric oxide oxidisers since they're probably just detoxifying
specialisations <- specialisations %>% 
  filter(nitrogen_specialisation != "nitric_oxide_oxidiser")

vip_df <- left_join(vip_members, input$taxonomy, by = c("MAG" = "genome")) %>%
  left_join(specialisations, by = c("MAG" = "genome")) %>%
  mutate(PhylumGlom = ifelse(!is.na(Phylum) & !(Phylum %in% phylum_levels), "other", Phylum)) %>% 
  mutate(PhylumLabel = ifelse(!is.na(PhylumGlom), PhylumGlom, "Environmental Predictor")) %>%
  mutate(Variable = fct_reorder(MAG, `VIP Score`, .desc = T))

ggplot(vip_df, aes(x = MAG, y = `VIP Score`)) +
  geom_bar(stat = "identity", aes(fill = PhylumLabel, group = `Environmental Variable`),
           position = "dodge",
           color = "grey20") +
  facet_wrap(~`Data Origin`) +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  ylab("VIP Score") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white"))

vip_specialisations.plot.df <- vip_df %>% 
  filter(`Data Origin` != "All") %>%
  mutate(`Data Origin` = factor(`Data Origin`, levels = c("All", habitat_levels))) %>%
  group_by(MAG, `Data Origin`) %>%
  mutate(mean_VIP_Score = mean(`VIP Score`)) %>% # Sum of MAG's VIP values across all modules within habitat
  ungroup() %>% 
  # select(-`Module Color`, -`VIP Score`, -`Environmental Variable`, -Pathways) %>% distinct() %>%
  mutate(MAG = fct_reorder(MAG, `mean_VIP_Score`, .desc = TRUE)) %>%
  ungroup() %>%
  rename(carbon_specialisation = specialisation) %>%
  pivot_longer(contains("specialisation"), names_to = "specialisation_type",
               values_to = "specialisation") %>% 
  mutate(specialisation_type = gsub("_", " ", specialisation_type))

vip_specialisations.plot <- ggplot(vip_specialisations.plot.df, 
                                   aes(x = MAG, y = `VIP Score`)) +
  geom_bar(stat = "identity", aes(fill = specialisation, group = `Environmental Variable`),
           position = position_dodge2(width = 0.5, preserve = "single", ),
           color = "grey20") +
  scale_fill_manual("Specialisation", values = all_special_colour, breaks = all_special_levels, labels = all_special_labels) +
  facet_grid(specialisation_type ~`Data Origin`, scales = "free_x", space = "free_x", switch = "y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  ylab("VIP Score") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.1)),
        strip.placement = "outside",
        plot.background = element_rect(fill = "white"))

vip_specialisations.plot


ggsave("~/Downloads/vip_specialisations_plot.png", plot = vip_specialisations.plot, device = "png", dpi = 300,
       height = 10, width = 25)
#### ====================================================================== ####

# Specialisations across habitat
#### ====================================================================== ####
habitat_specialisation_df <- trimmed_mean$trimmed_mean %>%
  right_join(input$sample_metadata %>% select(temporal_sample_id, Habitat__, DepthAvg__, 
                                              Year__), 
             by = "temporal_sample_id") %>%
  filter(genome %in% input$taxonomy$genome) %>%
  left_join(specialisations, by = "genome") %>%
  rename(carbon_specialisation = specialisation) %>%
  pivot_longer(contains("specialisation"), names_to = "specialisation_type",
               values_to = "specialisation") %>%
  group_by(specialisation, temporal_sample_id) %>%
  mutate(sum_Relative_Abundance = sum(relabund_of_recovered, na.rm = TRUE)) %>%
  ungroup() %>%
  select(temporal_sample_id, sum_Relative_Abundance, Habitat__, contains("specialisation")) %>%
  distinct() %>%
  filter(!is.na(specialisation)) %>%
  group_by(Habitat__, specialisation) %>%
  mutate(mean_specialisation = mean(sum_Relative_Abundance), 
         sd_specialisation = sd(sum_Relative_Abundance)) %>% 
  select(Habitat__, temporal_sample_id, contains("specialisation"), sum_Relative_Abundance) %>%
  mutate(specialisation_type = gsub("_", " ", specialisation_type)) %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels))

habitat_specialisations_plot <- ggplot(habitat_specialisation_df %>%
         select(-temporal_sample_id, -sum_Relative_Abundance) %>%
         distinct(), 
       aes(x = fct_reorder(specialisation, mean_specialisation, .desc = T),
             y = mean_specialisation)) +
  geom_bar(stat = "identity", aes(fill = specialisation),
           position = "dodge",
           color = "grey20") +
  scale_fill_manual("Specialisation", values = all_special_colour, breaks = all_special_levels[-15], labels = all_special_labels[-15]) +
  geom_errorbar(aes(ymin = mean_specialisation - sd_specialisation, 
                    ymax = mean_specialisation + sd_specialisation),
                width = 0.5) +
  # geom_point(data = habitat_specialisation_df, 
  #            aes(y = sum_Relative_Abundance),
  #            alpha = 0.5, shape = 21,
  #            position = "jitter") +
  facet_wrap(specialisation_type ~ Habitat__) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(limits = all_special_levels[-14], labels = all_special_labels[-15]) +
  ylab("Summed Relative Abundance of MAGs with Specialisation") + xlab("Specialisation") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),
        plot.background = element_rect(fill = "white"))

ggsave("~/Downloads/Habitat_specialisations_plot.png", 
       plot = habitat_specialisations_plot, device = "png", dpi = 300,
       height = 10, width = 20)
### ====================================================================== ####

# Specialisation turnover
#### ====================================================================== ####
# get subsetted otu tables
specialisation_lf <- specialisations %>% 
  rename(carbon_specialisation = specialisation) %>%
  pivot_longer(contains("specialisation"), names_to = "specialisation_type",
               values_to = "specialisation")
all_specialists <- specialisation_lf %>% 
  filter(specialisation != "generalist")

all_specialists_otu_table <- input_ra$otu_table %>%
  filter(genome %in% all_specialists$genome)

carbon_generalists <- specialisation_lf %>% 
  filter(specialisation == "generalist")

carbon_generalists_otu_table <- input_ra$otu_table %>%
  filter(genome %in% carbon_generalists$genome)

carbon_specialists <- specialisation_lf %>% 
  filter(specialisation_type %in% "carbon_specialisation") %>%
  filter(specialisation != "generalist")

carbon_specialists_otu_table <- input_ra$otu_table %>%
  filter(genome %in% carbon_specialists$genome)

nitrogen_specialists <- specialisation_lf %>% 
  filter(specialisation_type %in% "nitrogen_specialisation") %>%
  filter(!is.na(specialisation))

nitrogen_specialists_otu_table <- input_ra$otu_table %>%
  filter(genome %in% nitrogen_specialists$genome)

sulfur_specialists <- specialisation_lf %>% 
  filter(specialisation_type %in% "sulfur_specialisation") %>%
  filter(!is.na(specialisation))

sulfur_specialists_otu_table <- input_ra$otu_table %>%
  filter(genome %in% sulfur_specialists$genome)


# Dist Matrices
all_specialists_bc <- all_specialists_otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")
carbon_generalists_bc <- carbon_generalists_otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")
carbon_specialists_bc <- carbon_specialists_otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")
nitrogen_specialists_bc <- nitrogen_specialists_otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")
sulfur_specialists_bc <- sulfur_specialists_otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")
whole_community_bc <- input$otu_table[,-1] %>%
  t() %>%
  vegdist(., method = "bray")


# Convert to long format for plotting
all_dm_bc.lf <- as.matrix(whole_community_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  #mutate(SampleID_B = gsub("^X", "", SampleID_B)) %>%
  mutate(SpecialistAll = "Whole")
Cgen_dm_bc.lf <- as.matrix(carbon_generalists_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  mutate(SpecialistAll = "Carbon Generalists")
spec_dm_bc.lf <- as.matrix(all_specialists_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  mutate(SpecialistAll = "All Specialists")

Cspec_dm_bc.lf <- as.matrix(carbon_specialists_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  mutate(SpecialistAll = "Carbon Specialists")

Nspec_dm_bc.lf <- as.matrix(nitrogen_specialists_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  mutate(SpecialistAll = "Nitrogen Specialists")

Sspec_dm_bc.lf <- as.matrix(sulfur_specialists_bc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_A") %>%
  pivot_longer(names_to = "SampleID_B", 
               values_to = "Bray-Curtis", -SampleID_A) %>%
  mutate(SpecialistAll = "Sulfur Specialists")


## Plot

Specialisation_bc <- bind_rows(spec_dm_bc.lf, all_dm_bc.lf, Cgen_dm_bc.lf,
          Cspec_dm_bc.lf, Nspec_dm_bc.lf, Sspec_dm_bc.lf) %>%
  mutate(SpecialistAll = factor(SpecialistAll, 
                                levels = c("Whole", "Carbon Generalists", "All Specialists",
                                           "Carbon Specialists", "Nitrogen Specialists", 
                                           "Sulfur Specialists"))) %>%
  filter(SampleID_A != SampleID_B) %>%
  left_join(input$sample_metadata %>% 
              select(temporal_sample_id, Habitat__, Year__, DepthLumping),
            by = c("SampleID_A" = "temporal_sample_id" )) %>%
  rename(Habitat__A = Habitat__,
         Year__A = Year__,
         DepthLumpingA = DepthLumping) %>%
  left_join(input$sample_metadata %>% 
              select(temporal_sample_id, Habitat__, Year__, DepthLumping),
            by = c("SampleID_B" = "temporal_sample_id" )) %>%
  rename(Habitat__B = Habitat__,
         Year__B = Year__,
         DepthLumpingB = DepthLumping)

# Plot by habitat
specialisation_bc_hab.plot <- Specialisation_bc %>%
  filter(Habitat__A == Habitat__B) %>%
  mutate(Habitat__A = factor(Habitat__A, levels = habitat_levels)) %>%
  ggplot(aes(x = SpecialistAll, y = `Bray-Curtis`)) +
  geom_violin() +
  geom_point(aes(color = Habitat__A),position = "jitter", alpha = 0.05) +
  scale_color_manual("Habitat", values = colour_habitat, 
                     breaks = habitat_levels, labels = habitat_levels) +
  facet_wrap(~Habitat__A) + 
  xlab("Community") + ylab("Community Turnover Rate (Bray-Curtis)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
specialisation_bc_hab.plot

ggsave("~/Downloads/BC_specialisations_hab.png", 
       plot = specialisation_bc_hab.plot, 
       device = "png", dpi = 300,
       height = 10, width = 20)

# Plot by habitat and year
specialisation_bc_hab_year.plot <- Specialisation_bc %>%
  filter(Habitat__A == Habitat__B,
         Year__A == Year__B) %>%
  mutate(Habitat__A = factor(Habitat__A, levels = habitat_levels)) %>%
  ggplot(aes(x = SpecialistAll, y = `Bray-Curtis`)) +
  geom_violin() +
  geom_point(aes(color = Habitat__A),position = "jitter", alpha = 0.05) +
  scale_color_manual("Habitat", values = colour_habitat, 
                    breaks = habitat_levels, labels = habitat_levels) +
  facet_grid(Habitat__A ~Year__A, switch = "y") + 
  xlab("Community") + ylab("Community Turnover Rate (Bray-Curtis)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
specialisation_bc_hab_year.plot

ggsave("~/Downloads/BC_specialisations_hab_year.png", 
       plot = specialisation_bc_hab_year.plot, 
       device = "png", dpi = 300,
       height = 10, width = 20)

# Plot by habitat and Depth
specialisation_bc_hab_depth.plot <- Specialisation_bc %>%
  filter(Habitat__A == Habitat__B,
         DepthLumpingA == DepthLumpingB) %>%
  mutate(Habitat__A = factor(Habitat__A, levels = habitat_levels)) %>%
  ggplot(aes(x = SpecialistAll, y = `Bray-Curtis`)) +
  geom_violin() +
  geom_point(aes(color = Habitat__A),position = "jitter", alpha = 0.05) +
  scale_color_manual("Habitat", values = colour_habitat, 
                     breaks = habitat_levels, labels = habitat_levels) +
  facet_grid(Habitat__A ~DepthLumpingA, switch = "y") + 
  xlab("Community") + ylab("Community Turnover Rate (Bray-Curtis)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
specialisation_bc_hab_depth.plot

ggsave("~/Downloads/BC_specialisations_hab_depth.png", 
       plot = specialisation_bc_hab_depth.plot, 
       device = "png", dpi = 300,
       height = 10, width = 20)
#### ====================================================================== ####

