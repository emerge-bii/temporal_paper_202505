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
# Read in assembly analysis results
betanull.lf <- read_csv(file = paste0(outputs.fp, "/betanull.lf.csv"))

# Change variables into a factors
betanull.lf <- betanull.lf %>%
  mutate(across(starts_with("Habitat__"), ~ factor(.x, 
                            levels = c("Palsa", "Collapsed Palsa", 
                                       "Bog", "Fen")))) %>%
  mutate(across(starts_with("DepthLumping"),
                ~ factor(.x, levels =
                           c("0-9", "10-19", "20-29", "30-39", 
                             "40-49")))) %>%
  mutate(across(starts_with("DepthCodeMod."), 
                ~ factor(.x, levels = 
                           c("S","M", "D", "X0", "X1", "X2", 
                             "X3", "X4", "X5", "X6")))) %>%
  mutate(Assembly_Process = factor(Assembly_Process, 
                                        levels = c("Homogenous selection",
                                                   "Heterogenous selection",
                                                   "Homogenizing dispersal",
                                                   "Dispersal limitation and drift",
                                                   "Drift")))
# Read in betanull.lf.diff
betanull.lf.diff <- read_csv(file = paste0(outputs.fp, "/betanull.lf.diff.csv"))

# Read in Dylan's within-between analysis
bc_dist <- read_csv(file = here("data", "within-between-year-comp.csv"))

# Filter MainAutochamber.201107_E_1_1to4 out of metadata as it wasn't around when the assembly analysis was run

# Setup plotting parameters
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(15, 16, 17)
fill_habitat   <- colour_habitat
depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
habdepth_levels <- paste0(rep(c("P", "B", "F"), times = 4), " ", rep(depth_levels, each = 3))
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth
year_levels <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017")
year_labels <- year_levels
fill_year   <- RColorBrewer::brewer.pal(7, "OrRd")
colour_year <- fill_year
assembly_levels <- c("Homogenous selection", "Heterogenous selection",
                     "Homogenizing dispersal", "Dispersal limitation and drift",
                     "Drift")
assembly_labels <- assembly_levels
#colour_assembly <- c("#133253", "#738CA6","#1C4A00","#80BA5D", "#B1DF95") # blue/green
#colour_assembly <- c("#3F002E", "#BE7FAD","#871200","#DA3015", "#FF816D") # purple/red
#colour_assembly <- c("#6C358D", "#AB81C4","#085A4F","#2D8478", "#79BDB4") # purple/teal
colour_assembly <- c("#521168", "#8e318f","#BD4D0C","#FF8000", "#FADBAC") # purple/orange
fill_assembly <- colour_assembly
phylum_levels <- c(
  # Bacteria
  "p__Acidobacteriota",
  "p__Actinobacteriota",
  "p__Bacteroidota",
  "p__Chloroflexota",
  "p__Desulfobacterota",
  "p__Eremiobacterota",
  "p__Proteobacteria",
  "p__Verrucomicrobiota",
  # Archaea
  "p__Halobacteriota",
  "p__Methanobacteriota",
  "p__Thermoplasmatota",
  "p__Thermoproteota",
  "other",
  "Environmental Predictor"
)
phylum_labels <- c(
  # Bacteria
  "Acidobacteriota",
  "Actinobacteriota",
  "Bacteroidota",
  "Chloroflexota",
  "Desulfobacterota",
  "Eremiobacterota",
  "Proteobacteria",
  "Verrucomicrobiota",
  # Archaea
  "Halobacteriota",
  "Methanobacteriota",
  "Thermoplasmatota",
  "Thermoproteota",
  "other",
  "Environmental Predictor"
)
colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white", "grey30")
fill_phylum <- colour_phylum
metapathway_levels <- metapathway_groups %>% arrange(metapath) %>% `$`(metapath) %>% unique()
colour_metapathway <- RColorBrewer::brewer.pal(10, "Set3")
fill_metapathway <- colour_metapathway


#### ====================================================================== ####

  
# Define helper functions for plotting and subsetting dataframes
#### ====================================================================== ####
# This function translates habitat transitions so that order doesn't matter
habitat_transition_translate <- function(transitions) {
  case_when(
    transitions == "Palsa vs. Palsa" ~ "Palsa",
    transitions == "Palsa vs. Bog" ~ "Palsa vs. Bog",
    transitions == "Bog vs. Palsa" ~ "Palsa vs. Bog",
    transitions == "Bog vs. Bog" ~ "Bog",
    transitions == "Bog vs. Fen" ~ "Bog vs. Fen",
    transitions == "Fen vs. Bog" ~ "Bog vs. Fen",
    transitions == "Fen vs. Fen" ~ "Fen",
    transitions == "Palsa vs. Fen" ~ "Fen vs. Palsa",
    transitions == "Fen vs. Palsa" ~ "Fen vs. Palsa",
    TRUE ~ NA_character_
  )
}

#### ====================================================================== ####

# Calculate proportions of pairwise comparisons overall
#### ====================================================================== ####
# Proportions across all samples
betanull.lf %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total) %>%
  arrange(desc(Percent))

#### ====================================================================== ####

# Summary plot by habitat
#### ====================================================================== ####
intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              DistanceMeasure = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
# Cross-depth comparisons allowed
# betanull.lf.hab.plot <- betanull.lf %>%
#   select(Site1, Site2, Assembly_Process, 
#          starts_with("Habitat__"), starts_with("Depth"),  BetaNTI, RCBC) %>%
#   filter(Habitat__.Site1 != "Collapsed Palsa",
#          Habitat__.Site2 != "Collapsed Palsa") %>%
#   mutate(Habitat__.Site1  = fct_drop(Habitat__.Site1),
#          Habitat__.Site2  = fct_drop(Habitat__.Site2)) %>%
#   # Filter so that we don't keep comparisons between different depths
#   filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
#   # Filter so that the only transitions we're considering are
#   # thaw state transitions
#   filter(as.numeric(Habitat__.Site1) == 
#            as.numeric(fct_shift(Habitat__.Site2, n = 1)) |
#            as.numeric(Habitat__.Site2) == 
#            as.numeric(fct_shift(Habitat__.Site1, n = -1)) |
#            Habitat__.Site1 == Habitat__.Site2) %>%
#   unite(Habitat__.Site1, Habitat__.Site2, col = "HabitatTransition", sep = " vs. ",
#         remove = FALSE) %>%
#   mutate(HabitatTransition = ifelse(Habitat__.Site1 == Habitat__.Site2,
#                                     as.character(Habitat__.Site1), HabitatTransition)) %>%
#   group_by(HabitatTransition) %>%
#   mutate(HabitatTransition = ifelse(HabitatTransition == "Fen vs. Palsa", 
#                                     "Palsa vs. Fen", HabitatTransition)) %>%
#   mutate(HabitatTransition = factor(HabitatTransition, 
#                                     levels = c("Palsa",
#                                                "Palsa vs. Bog", "Bog",
#                                                "Bog vs. Fen", "Fen", "Palsa vs. Fen"))) %>%
#   filter(!grepl("vs.", HabitatTransition)) %>%
#   pivot_longer(cols = matches(c("BetaNTI", "RCBC")),
#                names_to = "DistanceMeasure", values_to = "Distance")
# 
# hab.plot <- betanull.lf.hab.plot %>% 
#   ggplot(aes(x = HabitatTransition, y = Distance)) +
#   geom_point(aes(color = Assembly_Process),
#              position = "jitter", alpha = 0.3) +
#   geom_violin(fill = NA) +
#   scale_color_manual(name = "Assembly Process", values = colour_assembly, 
#                      breaks = assembly_levels, 
#                      labels = assembly_labels) +
#   geom_hline(data = intercept_guide,
#              aes(yintercept = Distance),
#              linetype = "dashed") +
#   facet_wrap(~DistanceMeasure, scales = "free_y") +
#   xlab("") + ylab("Turnover Distance") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none")#+
#   #ggtitle("Assembly Processes underlying turnover between samples")
# hab.plot
# 
# ggsave(hab.plot, 
#        filename = paste0(figures.fp, "/Habitat_fullbetanull_violin_plot_nocrossdepth.png"),
#        width = 10, height = 5, dpi = 400)
# 
# hab.plot.blank <- betanull.lf.hab.plot %>% 
#   ggplot(aes(x = HabitatTransition, y = Distance)) +
#   geom_point(aes(color = Assembly_Process),
#              position = "jitter", alpha = 0) +
#   #geom_violin(fill = NA) +
#   #geom_boxplot(fill = NA) +
# #  scale_color_viridis(name = "Assembly Process", discrete = TRUE) +
#   scale_color_manual(name = "Assembly Process", values = colour_assembly, 
#                      breaks = assembly_levels, 
#                      labels = assembly_labels) +
#   geom_hline(data = intercept_guide,
#              aes(yintercept = Distance),
#              linetype = "dashed") +
#   facet_wrap(~DistanceMeasure, scales = "free_y") +
#   xlab("") + ylab("Turnover Distance") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none")#+
# #ggtitle("Assembly Processes underlying turnover between samples")
# hab.plot.blank
# 
# ggsave(hab.plot.blank, 
#        filename = paste0(figures.fp, "/Habitat_fullbetanull_blank_plot.png"),
#        width = 10, height = 5, dpi = 400)

# Restrict to depth self-comparisions only
hab.plot.prop.df <- betanull.lf %>%
  filter(Habitat__.Site1 != "Collapsed Palsa",
         Habitat__.Site2 != "Collapsed Palsa") %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>% # Depth self comparisons only
  mutate(Habitat__.Site1  = fct_drop(Habitat__.Site1),
         Habitat__.Site2  = fct_drop(Habitat__.Site2)) %>%
  unite(Habitat__.Site1, Habitat__.Site2, col = "HabitatTransition", sep = " vs. ",
        remove = FALSE) %>% 
  mutate(HabitatTransition = habitat_transition_translate(HabitatTransition)) %>% 
  group_by(HabitatTransition) %>% 
  mutate(HabitatTransition = factor(HabitatTransition, 
                                    levels = c("Palsa",
                                               "Palsa vs. Bog", "Bog",
                                               "Bog vs. Fen", "Fen", "Fen vs. Palsa"))) %>%
  select(Site1, Site2, Assembly_Process, DepthLumping.Site1, DepthLumping.Site2, HabitatTransition) %>%
  group_by(HabitatTransition, Assembly_Process) %>%
  tally() %>% 
  ungroup() %>% group_by(HabitatTransition) %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total)
  
hab.plot.prop <- hab.plot.prop.df %>%
  filter(!grepl("vs.", HabitatTransition)) %>% # remove cross-habitat comparisons
  ggplot(aes(x= HabitatTransition, y = Percent, 
             fill = Assembly_Process)) +
  facet_wrap(~HabitatTransition, shrink = TRUE, drop = TRUE, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")
hab.plot.prop

hab.plot.prop.legend <- get_legend(hab.plot.prop)
hab.plot.prop <- hab.plot.prop + theme(legend.position = "none")

ggsave(hab.plot.prop, 
       filename = paste0(figures.fp, "/Habitat_prop_no_crossdepthcomp.png"),
       width = 15, height = 3, dpi = 400)

comb_plot <- plot_grid(hab.plot.prop, NULL, hab.plot, hab.plot.prop.legend,
                       nrow = 4,
                       labels = c("A", NULL, "B", NULL),
                       label_size = 12, 
                       rel_heights = c(0.9,0.1,1, 0.2),
                       align = "hv", 
                       axis = "lb")
comb_plot
ggsave(comb_plot, 
       filename = paste0(figures.fp, "/Habitat_fullbetanull_comb_plot_nocrossdepth.png"),
       width = 10, height = 8, dpi = 400)

# With habitat cross-comparisons
hab.plot.prop.trans <- hab.plot.prop.df %>%
  ggplot(aes(x= HabitatTransition, y = Percent, 
             fill = Assembly_Process)) +
  facet_wrap(~HabitatTransition, shrink = TRUE, drop = TRUE, scales = "free_x",ncol = 2) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")
hab.plot.prop.trans

ggsave(hab.plot.prop.trans, 
       filename = paste0(figures.fp, "/Habitat_prop_no_crossdepthcomp_habitat_transitions.png"),
       width = 15, height = 10, dpi = 400)

#### ====================================================================== ####

# By Year
#### ====================================================================== ####
# Over the years
# Version with distance in years
year.plot.prop.df <- betanull.lf.diff %>%
  filter(Depth_diff < 10) %>%
  group_by(Habitat_comp, Year_diff, Assembly_Process) %>%
  filter(Habitat_comp %in% c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")) %>%
  tally() %>%
  mutate(Percent = 100* n/sum(n),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels),
         Habitat_comp = gsub("(.{1,})(_)(.{1,})", "\\1", Habitat_comp),
         Habitat_comp = factor(Habitat_comp, levels = habitat_levels))

year.plot.prop <-  year.plot.prop.df %>%
  ggplot(aes(x = Year_diff, y = Percent, fill = Assembly_Process)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Habitat_comp, strip.position = "left", ncol = 1) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  xlab("Number of years between samples") + ylab("Percent of Pairwise Comparisons") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        plot.margin = margin(2,2,2,2, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.placement = "outside",
        strip.text.x = element_text(size = rel(1)),
        strip.text.y = element_text(face = "bold", size = rel(1)),
        legend.position = "bottom")
year.plot.prop 

ggsave(year.plot.prop, 
       filename = paste0(figures.fp, "/year_prop_by_habitat_plot_nocrossdepthcomp.png"),
       width = 16, height = 9, dpi = 600)


# Cross Depth comparisons allowed
# We do this with depth cross-comparisons allowed because otherwise many observation points would have no data
yr_betanull <- betanull.lf %>%
  #filter(abs(DepthAvg__.Site1 - DepthAvg__.Site2) < 10) %>%
  select(Site1, Site2, Assembly_Process, Habitat__.Site1, Habitat__.Site2,
         Year__.Site1, Year__.Site2) %>%
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  mutate(SameYrDiff = ifelse(Year__.Site1 == Year__.Site2, "Same", "Different")) %>%
  mutate(Year2011 = ifelse(Year__.Site1 == 2011 | Year__.Site2 == 2011, 2011, NA),
         Year2012 = ifelse(Year__.Site1 == 2012 | Year__.Site2 == 2012, 2012, NA),
         Year2013 = ifelse(Year__.Site1 == 2013 | Year__.Site2 == 2013, 2013, NA),
         Year2014 = ifelse(Year__.Site1 == 2014 | Year__.Site2 == 2014, 2014, NA),
         Year2015 = ifelse(Year__.Site1 == 2015 | Year__.Site2 == 2015, 2015, NA),
         Year2016 = ifelse(Year__.Site1 == 2016 | Year__.Site2 == 2016, 2016, NA),
         Year2017 = ifelse(Year__.Site1 == 2017 | Year__.Site2 == 2017, 2017, NA))

year_year_comparisons <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)

  filt_year <- yr_betanull %>%
    filter(SameYrDiff == "Different") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x)

  return(filt_year)
})

names(year_year_comparisons) <- as.character(c(2011:2017))



same_year_comparisons <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)

  filt_year <- yr_betanull %>%
    filter(SameYrDiff == "Same") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x,
           YrComp = as.character(x))

  return(filt_year)
})


# For rectangle box around self-comparisons
rectangle_box <- year_year_comparisons %>%
  reduce(bind_rows) %>%
  ungroup() %>%
  select(Habitat__.Site1, RefYr, YrComp) %>%
  distinct() %>%
  mutate(xmin = factor(RefYr, levels = as.character(c(2011:2017))),
         xmax = factor(RefYr + 1, levels = as.character(c(2011:2018)))) %>%
  mutate(xmin = as.numeric(xmin)-0.5,
         xmax = as.numeric(xmax)-0.5) %>%
  select(-YrComp) %>% distinct() %>%
  mutate(Habitat__.Site1 = factor(Habitat__.Site1, levels = habitat_levels))


# Combine dataframes into 1 for plotting
year_to_year_comparison.df <- year_year_comparisons %>%
  reduce(bind_rows) %>%
  bind_rows(same_year_comparisons %>%
              reduce(bind_rows))  %>%
  select(-n) %>%
  #  filter(Habitat__.Site1 == "Palsa") %>%
  #  pivot_wider(names_from = Assembly_Process, values_from = prop, values_fill = 0) %>%
  #  group_by(Habitat__.Site1, RefYr) %>%
  #  nest() %>%
  #  pluck(3,1) %>%
  mutate(RefYr = factor(RefYr),
         YrComp = factor(YrComp),
         Habitat__.Site1 = factor(Habitat__.Site1, levels = habitat_levels),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels)) 


year_yearcomp.plot.prop <- year_to_year_comparison.df %>%
  ggplot(aes(x = YrComp, y = prop)) +
  facet_grid(Habitat__.Site1~RefYr, switch = "x") +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", width = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  # cover up self-comparisons
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 1,
            fill = "grey20", color = "grey20", linewidth = 0.8,
            data = rectangle_box,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly,
                    breaks = assembly_levels,
                    labels = assembly_labels) +
  xlab("Year of Comparison") + ylab("Proportion of Assembly Process") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.placement = "outside",
        #aspect.ratio=1,
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.75)))
year_yearcomp.plot.prop

ggsave(year_yearcomp.plot.prop,
       filename = paste0(figures.fp, "/year_yearcomp_prop_by_habitat_plot_crossdepthcomp.png"),
       width = 10, height = 8, dpi = 600)



# Line plots
year_year_comparisons %>%
  reduce(bind_rows) %>%
  bind_rows(same_year_comparisons %>%
              reduce(bind_rows))  %>%
  select(-n) %>%
  mutate(RefYr = factor(RefYr),
         YrComp = factor(YrComp)) %>%
  ggplot(aes(x = YrComp, y = prop)) +
  geom_point(aes(color = Assembly_Process)) +
  facet_grid(Habitat__.Site1~RefYr) +
  geom_line(aes(group = Assembly_Process, color = Assembly_Process), stat = "identity", width = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 1,
            fill = NA, color = "black", linewidth = 0.8,
            data = rectangle_box,
            inherit.aes = FALSE) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly,
                    breaks = assembly_levels,
                    labels = assembly_labels) +
  xlab("Year of Comparison") + ylab("Proportion of Assembly Process") +
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background = element_rect(fill = "white", color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))


year_yearcomp.plot.prop_box <- year_year_comparisons %>%
  reduce(bind_rows) %>%
  mutate(RefYr = as.numeric(as.character(RefYr))) %>%
  ggplot(aes(x = RefYr, y = YrComp, alpha = prop)) +
  geom_tile(aes(fill = Assembly_Process)) +
  facet_grid(Habitat__.Site1 ~ Assembly_Process) +
  theme_bw() +
  scale_x_continuous(breaks = c(2011:2017), expand = c(0,0)) +
  # scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly,
                    breaks = assembly_levels,
                    labels = assembly_labels) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly,
                     breaks = assembly_levels,
                     labels = assembly_labels) +
  theme(aspect.ratio=1,
        strip.background = element_rect(fill = "white", color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(year_yearcomp.plot.prop_box,
       filename = paste0(figures.fp, "/year_yearcomp_prop_by_habitat_box_plot_crossdepthcomp.png"),
       width = 20, height = 9, dpi = 600)


# Note, when cross-depth comparisons are disallowed, not enough data in each group
# to reliably calculate year-to-year transisions (some groups only  have one oberservation at each depth)
# or don't have observations at depths (ex. 2014 palsa)
#### ====================================================================== ####

#### ====================================================================== ####
# Cross Depth comparisons excluded
# We do this with depth cross-comparisons allowed because otherwise many observation points would have no data
yr_betanull_nocd <- betanull.lf %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>% # may not be the right way to do this...
  #filter(abs(DepthAvg__.Site1 - DepthAvg__.Site2) < 10) %>%
  select(Site1, Site2, Assembly_Process, Habitat__.Site1, Habitat__.Site2,
         Year__.Site1, Year__.Site2, DepthLumping.Site1) %>%
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  mutate(SameYrDiff = ifelse(Year__.Site1 == Year__.Site2, "Same", "Different")) %>%
  mutate(Year2011 = ifelse(Year__.Site1 == 2011 | Year__.Site2 == 2011, 2011, NA),
         Year2012 = ifelse(Year__.Site1 == 2012 | Year__.Site2 == 2012, 2012, NA),
         Year2013 = ifelse(Year__.Site1 == 2013 | Year__.Site2 == 2013, 2013, NA),
         Year2014 = ifelse(Year__.Site1 == 2014 | Year__.Site2 == 2014, 2014, NA),
         Year2015 = ifelse(Year__.Site1 == 2015 | Year__.Site2 == 2015, 2015, NA),
         Year2016 = ifelse(Year__.Site1 == 2016 | Year__.Site2 == 2016, 2016, NA),
         Year2017 = ifelse(Year__.Site1 == 2017 | Year__.Site2 == 2017, 2017, NA))

year_year_comparisons_nocd <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)
  
  filt_year <- yr_betanull_nocd %>%
    filter(SameYrDiff == "Different") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x)
  
  return(filt_year)
})

names(year_year_comparisons_nocd) <- as.character(c(2011:2017))



same_year_comparisons_nocd <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)
  
  filt_year <- yr_betanull_nocd %>%
    filter(SameYrDiff == "Same") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x,
           YrComp = as.character(x))
  
  return(filt_year)
})


rectangle_box <- year_year_comparisons_nocd %>%
  reduce(bind_rows) %>%
  ungroup() %>%
  select(Habitat__.Site1, RefYr, YrComp) %>%
  distinct() %>%
  mutate(xmin = factor(RefYr, levels = as.character(c(2011:2017))),
         xmax = factor(RefYr + 1, levels = as.character(c(2011:2018)))) %>%
  mutate(xmin = as.numeric(xmin)-0.5,
         xmax = as.numeric(xmax)-0.5) %>%
  select(-YrComp) %>% distinct() %>%
  mutate(Habitat__.Site1 = factor(Habitat__.Site1, levels = habitat_levels)) 




year_to_year_comparison_nocd.df <- year_year_comparisons_nocd %>%
  reduce(bind_rows) %>%
  bind_rows(same_year_comparisons_nocd %>%
              reduce(bind_rows))  %>%
  select(-n) %>%
  #  filter(Habitat__.Site1 == "Palsa") %>%
  #  pivot_wider(names_from = Assembly_Process, values_from = prop, values_fill = 0) %>%
  #  group_by(Habitat__.Site1, RefYr) %>%
  #  nest() %>%
  #  pluck(3,1) %>%
  mutate(RefYr = factor(RefYr),
         YrComp = factor(YrComp),
         Habitat__.Site1 = factor(Habitat__.Site1, levels = habitat_levels),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels)) 
year_yearcomp.plot_nocd.prop <- year_to_year_comparison_nocd.df %>%
  ggplot(aes(x = YrComp, y = prop)) +
  facet_grid(Habitat__.Site1~RefYr) +
  geom_bar(aes(fill = Assembly_Process), stat = "identity", width = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  # cover up self-comparisons
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 1,
            fill = "grey20", color = "grey20", linewidth = 0.8,
            data = rectangle_box,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly,
                    breaks = assembly_levels,
                    labels = assembly_labels) +
  xlab("Year of Comparison") + ylab("Proportion of Assembly Process") +
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background = element_rect(fill = "white", color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))
year_yearcomp.plot_nocd.prop

ggsave(year_yearcomp.plot.prop,
       filename = paste0(figures.fp, "/year_yearcomp_prop_by_habitat_plot_nocrossdepthcomp.png"),
       width = 9, height = 16, dpi = 600)




# Line plots
year_year_comparisons_nocd <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)
  
  filt_year <- yr_betanull_nocd %>%
    filter(SameYrDiff == "Different") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, DepthLumping.Site1, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x)
  
  return(filt_year)
})

names(year_year_comparisons_nocd) <- as.character(c(2011:2017))



same_year_comparisons_nocd <- lapply(c(2011:2017), function(x) {
  yearcolumn = paste0("Year", x)
  
  filt_year <- yr_betanull_nocd %>%
    filter(SameYrDiff == "Same") %>%
    filter(.data[[yearcolumn]] == x) %>%
    select(!matches(yearcolumn)) %>%
    unite(starts_with("Year2"), col = "YrComp", na.rm = T) %>%
    group_by(Habitat__.Site1, YrComp, DepthLumping.Site1, Assembly_Process) %>% tally() %>%
    mutate(prop = n/sum(n)) %>%
    mutate(RefYr = x,
           YrComp = as.character(x))
  
  return(filt_year)
})


Habitat_depth_levels <- c(paste(rep("P", times = 4), c("0-9", "10-19", "20-29", "30-39")),
                         paste(rep("B", times = 4), c("0-9", "10-19", "20-29", "30-39")),
                         paste(rep("F", times = 4), c("0-9", "10-19", "20-29", "30-39")))
habitat_background <- data.frame(Habitat__.Site1 = habitat_levels,
                                 DepthLumping.Site1 = rep(depth_levels, 3),
                                 xmin = rep("2011", 12),
                                 xmax = rep("2018", 12)) %>%
  mutate(HabitatDepth = factor(paste0(substr(Habitat__.Site1,1,1), " ", DepthLumping.Site1),
                                        levels = Habitat_depth_levels))






year_year_comparisons_no_cross_depth <- year_year_comparisons_nocd %>%
  reduce(bind_rows) %>%
  bind_rows(same_year_comparisons_nocd %>%
              reduce(bind_rows))  %>%
  select(-n) %>%
  mutate(RefYr = factor(RefYr, levels = as.character(c(2011:2017))),
         YrComp = factor(YrComp, levels = as.character(c(2011:2017))),
         HabitatDepth = factor(paste0(substr(Habitat__.Site1,1,1), " ", DepthLumping.Site1),
                               levels = Habitat_depth_levels)) %>%
  ggplot(aes(x = YrComp, y = prop)) +
  geom_point(aes(color = Assembly_Process)) +
  facet_grid(HabitatDepth~RefYr) +
  geom_rect(aes(xmin = xmin, xmax = xmax, fill = Habitat__.Site1), ymin = 0, ymax = 1,
            alpha = 0.2,
            data = habitat_background,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "Habitat", values = colour_habitat, 
                    breaks = habitat_levels, labels = habitat_levels) +
  geom_line(aes(group = Assembly_Process, color = Assembly_Process), stat = "identity", width = 1) +
  scale_x_discrete(expand = c(0,0), breaks = c(2011, 2013, 2015, 2017)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 1,
            fill = NA, color = "black", linewidth = 0.8,
            data = rectangle_box,
            inherit.aes = FALSE) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly,
                     breaks = assembly_levels,
                     labels = assembly_labels) +
  xlab("Year of Comparison") + ylab("Proportion of Assembly Process") +
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.y = element_text(angle = 0,hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(year_year_comparisons_no_cross_depth,
       filename = paste0(figures.fp, "/year_yearcomp_prop_by_habitat_plot_NOcrossdepthcomp.png"),
       width = 20, height = 10, dpi = 600)
#### ====================================================================== ####


# Depth
#### ====================================================================== ####
# Discrete depth distances by 5cm
depth.plot.prop.df <- betanull.lf.diff %>%
  mutate(Depth_diff_lumped = cut(Depth_diff, breaks = seq(from = 0, to = 35, by = 5), include.lowest = TRUE,
                                 labels = c("0",  "5", "10", "15", "20", "25", "30"))) %>%
  group_by(Habitat_comp, Depth_diff_lumped, Assembly_Process) %>%
  filter(Habitat_comp %in% c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")) %>%
  tally() %>%
  mutate(Percent = 100* n/sum(n),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels),
         Habitat_comp = gsub("(.{1,})(_)(.{1,})", "\\1", Habitat_comp),
         Habitat_comp = factor(Habitat_comp, levels = habitat_levels),
         Depth_diff_lumped = as.numeric(as.character(Depth_diff_lumped))) %>% 
  ungroup() %>%
  group_by(Habitat_comp, Depth_diff_lumped) %>%
  mutate(total_depthlump = sum(n)) %>%
  ungroup() %>%
  mutate(Depth_diff_lumped_label = paste0(Depth_diff_lumped, " (n = ", total_depthlump, ")"))
depth.plot.prop <-   depth.plot.prop.df %>% 
  ggplot(aes(x = Depth_diff_lumped, y = Percent, fill = Assembly_Process)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Habitat_comp, strip.position = "left", ncol = 1) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  xlab("Number of cm between samples") + ylab("Percent of Pairwise Comparisons") +
  scale_x_continuous(expand = c(0,0), n.breaks = 7) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12), 
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = rel(1)),
        strip.text.y = element_text(face = "bold", size = rel(1)),
        legend.position = "bottom")
depth.plot.prop 


# Depth discrete
depthlumping.plot.prop.df <- betanull.lf %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  group_by(Habitat__.Site1, DepthLumping.Site1, Assembly_Process) %>%
  tally() %>%
  mutate(Percent = 100* n/sum(n),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels)) %>% 
  ungroup() %>%
  group_by(Habitat__.Site1, DepthLumping.Site1) %>%
  mutate(total_depthlump = sum(n)) %>%
  ungroup() %>%
  mutate(Depthlumped_label = ifelse(Assembly_Process == "Homogenizing dispersal", paste0("(n = ", total_depthlump, ")"), NA)) %>%
  mutate(Habitat__.Site1 = factor(Habitat__.Site1, levels = habitat_levels))


depthlumping.plot.prop <-   depthlumping.plot.prop.df %>% 
  ggplot(aes(x = DepthLumping.Site1, y = Percent, fill = Assembly_Process)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Habitat__.Site1, strip.position = "left", ncol = 1) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  #geom_text(aes(label = Depthlumped_label), y = 50) +
  xlab("Depth Interval") + ylab("Percent of Pairwise Comparisons") +
  ggtitle("Depth Comparisons Assembly Processes") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title = element_text(size = rel(3)),
        axis.text = element_text(size = rel(2)),
        panel.grid.major = element_blank(),panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(3,3,3,3), "lines"),
        strip.placement = "outside",
        panel.spacing = unit(2, "lines"),
        strip.text.x = element_text(size = rel(3)),
        strip.text.y = element_text(face = "bold", size = rel(3)),
        plot.background = element_rect(color = "white"),
        legend.position = "none")
depthlumping.plot.prop
ggsave(depthlumping.plot.prop, 
              filename = paste0(figures.fp, "/depth_prop_by_habitat_plot.png"),
              width = 20, height = 11, dpi = 600)


#### ====================================================================== ####

# Mantel - BetaNTI
#### ====================================================================== ####
betaNTI <- betanull.lf %>% select(Site1, Site2, BetaNTI) %>%
  pivot_wider(names_from = Site2, values_from = BetaNTI) %>%
  column_to_rownames("Site1") %>% as.matrix()

(rownames(betaNTI)) == (input_ra$sample_metadata$temporal_sample_id)

RCBC <- betanull.lf %>% select(Site1, Site2, RCBC.nona) %>%
  pivot_wider(names_from = Site2, values_from = RCBC.nona) %>%
  column_to_rownames("Site1") %>% as.matrix()

(rownames(RCBC)) == (input_ra$sample_metadata$temporal_sample_id)



library(ecodist)

metadata_filt <- input_ra$sample_metadata %>% 
  select(temporal_sample_id, Habitat__, T_soil.deg_C, 
         #T_air.deg_C,
         samplingdate_mean_AirTemperature, 
         DepthLumping, Year__) %>%
  filter(Habitat__ == "Palsa") %>%
  na.omit() %>%
  column_to_rownames("temporal_sample_id")

euc_soilT <- dist(metadata_filt %>% select(T_soil.deg_C))
euc_airT <- dist(metadata_filt %>% select(matches("[Aa]ir")))
euc_DepthLumping <- as.matrix(dist(model.matrix(~0 + as.character(metadata_filt$DepthLumping)), method = "binary"))
euc_Year <- as.matrix(dist(model.matrix(~0 + as.character(metadata_filt$Year__)), method = "binary"))

betaNTI_var <- betaNTI[match(rownames(as.matrix(euc_airT)), rownames(betaNTI)), 
                       match(colnames(as.matrix(euc_airT)), colnames(betaNTI))]

palsa.MRM <- MRM(as.dist(betaNTI_var) ~ as.dist(euc_Year) + as.dist(euc_DepthLumping) + euc_airT + euc_soilT)
palsa.MRM


# Create a function for all samples
habitat <- "Palsa"
depth <- "0-9"
test_time_correlation <- function(habitat, depth){
  metadata_filt <- input_ra$sample_metadata %>% 
    select(temporal_sample_id, Habitat__, 
           DepthLumping, Year__) %>%
    filter(Habitat__ == habitat) %>%
    filter(DepthLumping == depth) %>%
    na.omit() %>%
    # Add in a Time since 2011 column
    mutate(TimeSince2011 = Year__-2011) %>%
    column_to_rownames("temporal_sample_id")

  # Variants on year to test
  # euc_Year: Simple euclidean distance in years
  euc_Year <- dist(metadata_filt %>% select(Year__), method = "euclidean")
  # bin_Year: Whether two samples are from the same (0) or different years (1), no continuous distance 
  bin_Year <- as.matrix(dist(model.matrix(~0 + as.character(metadata_filt$Year__)), method = "binary"))
  rownames(bin_Year) <- rownames(metadata_filt)
  colnames(bin_Year) <- rownames(metadata_filt)
  # euc_2011_baseline: Whether 2 samples are the same distance from the year 2011 or not; Samples that are a similar distance from
  # 2011 are closer to 0, samples that are dissimilarly distant from 2011 (eg 2012 and 2017) are considered far apart. These results
  # should be functionally similar to euc_Year except they use an arbitrary 3rd year as the baseline
  euc_2011_baseline <- dist(metadata_filt %>% select(TimeSince2011), method = "euclidean")
  
  # Add back in reciprocal comparisons
  betaNTI_var <- betaNTI[match(rownames(as.matrix(euc_Year)), rownames(betaNTI)), 
                         match(colnames(as.matrix(euc_Year)), colnames(betaNTI))]
  
  # Test if Year distances and betaNTI are correlated
  euc_Year_mantel <- ecodist::mantel(as.dist(betaNTI_var) ~ euc_Year, mrank = T)
  bin_Year_mantel <- ecodist::mantel(as.dist(betaNTI_var) ~ as.dist(bin_Year), mrank = T)
  #euc_2011_baseline_mantel <- ecodist::mantel(as.dist(betaNTI_var) ~ euc_2011_baseline)
  
  df_bnti_return <- data.frame(Habitat__ = habitat, DepthLumping = depth,
                          MantelX = c("EucYear", "BinYear"),
                          MantelY = "BetaNTI",
    bind_rows(euc_Year_mantel, bin_Year_mantel))
  
  
  # Add back in reciprocal comparisons
  RCBC_var <- RCBC[match(rownames(as.matrix(euc_Year)), rownames(RCBC)), 
                   match(colnames(as.matrix(euc_Year)), colnames(RCBC))]
  
  # Test if Year distances and betaNTI are correlated
  euc_Year_mantel <- ecodist::mantel(as.dist(RCBC_var) ~ euc_Year, mrank = T)
  bin_Year_mantel <- ecodist::mantel(as.dist(RCBC_var) ~ as.dist(bin_Year), mrank = T)
  #euc_2011_baseline_mantel <- ecodist::mantel(as.dist(RCBC_var) ~ euc_2011_baseline)
  
  df_RCBC_return <- data.frame(Habitat__ = habitat, DepthLumping = depth,
                               MantelX = c("EucYear", "BinYear"),
                               MantelY = "RCBC",
                               bind_rows(euc_Year_mantel, bin_Year_mantel))
  df_return <- bind_rows(df_bnti_return, df_RCBC_return)
  return(df_return)
}

# Create mantels with time. 
# Note:
# mantelr	
# pval1	
# one-tailed p-value (null hypothesis: r <= 0).
# 
# pval2	
# one-tailed p-value (null hypothesis: r >= 0).
# 
# pval3	
# two-tailed p-value (null hypothesis: r = 0).
# 
# llim	
# lower confidence limit.
# 
# ulim	
# upper confidence limit.
# eukyear = distances continuous over time (directional)
#binyear = distances year is same as comparison or different. (non-directional)
mantel_with_time <- input_ra$sample_metadata %>% 
  select(Habitat__, DepthLumping) %>%
  distinct() %>%
  rename(habitat = Habitat__, depth = DepthLumping) %>%
  pmap(., test_time_correlation) %>%
  reduce(bind_rows) %>%
#  filter(MantelX != "BinYear") %>%
  mutate(padj = p.adjust(pval3, method = "fdr"),
         sig = ifelse(padj < 0.05, "significant", "not significant"),
         sig_bin = ifelse(padj < 0.05, 1, 0)) %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))

mantel_with_time_plotdf <- mantel_with_time %>%
  filter(MantelX != "BinYear") %>%
  select(Habitat__, DepthLumping, MantelY, mantelr, sig_bin) %>%
  pivot_longer(matches(c("mantelr", "sig_bin")), names_to = "MeasureType", values_to = "value") %>%
  mutate(StochDet = ifelse(MantelY == "BetaNTI", "Deterministic", "Stochastic"),
         y = ifelse(MeasureType == "sig_bin", 1, 2)) %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))
value_labels <- mantel_with_time_plotdf %>%
  distinct(MeasureType, y) %>%
  mutate(MeasureType = ifelse(MeasureType == "mantelr", "Mantel Corr",
                              "Significance\n(<0.05 = blue)")) %>%
  mutate(Habitat__ = "Palsa", DepthLumping = "0-9", StochDet = "Deterministic") %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))

mantel_with_time.plot <- mantel_with_time_plotdf %>%
  ggplot(aes(x = StochDet, y = y)) +
  geom_tile(aes(fill = value)) +
  #geom_tile(aes(fill = mantelr)) +
  geom_text(data = value_labels, aes(label = MeasureType), size = rel(3), color = "grey20") +
  scale_fill_gradient2(name = "Mantel Correlation", 
                       limits = c(-1, 1)) +
  facet_grid(DepthLumping ~ Habitat__, switch = "y") +
  xlab("Stochastic Or Deterministic Correlation with Time") +
  ylab("") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank() 
        )

ggsave(mantel_with_time.plot, 
       filename = paste0(figures.fp, "/mantel_with_time.png"),
       width = 10, height = 10, dpi = 300)

# Write table of mantels out for sharing data
write.table(mantel_with_time, 
            file = paste0(outputs.fp, "/mantel_with_time.txt"),
            sep = "\t", row.names = F)
#### ====================================================================== ####
habitat = "Fen"
depth = "10-19"
temperature_metric = "T_air.deg_C"
test_temperature_correlation <- function(habitat, depth, temperature_metric){
  metadata_filt <- input_ra$sample_metadata %>% 
    select(temporal_sample_id, Habitat__, 
           DepthLumping, Year__, matches(temperature_metric)) %>%
    filter(Habitat__ == habitat) %>%
    filter(DepthLumping == depth) %>%
    na.omit() %>%
    # Add in a Time since 2011 column
    mutate(TimeSince2011 = Year__-2011) %>%
    column_to_rownames("temporal_sample_id")
  
  # Variants on year to test
  # euc_Year: Simple euclidean distance of temperature
  euc_temp <- dist(metadata_filt %>% select(matches(temperature_metric)), method = "euclidean")

  # Add back in reciprocal comparisons
  betaNTI_var <- betaNTI[match(rownames(as.matrix(euc_temp)), rownames(betaNTI)), 
                         match(colnames(as.matrix(euc_temp)), colnames(betaNTI))]
  
  # Test if Year distances and betaNTI are correlated
  euc_temp_mantel <- ecodist::mantel(as.dist(betaNTI_var) ~ euc_temp, mrank = T)
  #euc_2011_baseline_mantel <- ecodist::mantel(as.dist(betaNTI_var) ~ euc_2011_baseline)
  
  df_bnti_return <- data.frame(Habitat__ = habitat, DepthLumping = depth,
                               temperature_metric = temperature_metric,
                               MantelX = c("EucTemp"),
                               MantelY = "BetaNTI",
                               bind_rows(euc_temp_mantel))
  
  
  # Add back in reciprocal comparisons
  RCBC_var <- RCBC[match(rownames(as.matrix(euc_temp)), rownames(RCBC)), 
                   match(colnames(as.matrix(euc_temp)), colnames(RCBC))]
  if(sd(RCBC_var, na.rm = T) == 0) {
    print("standard deviation of matrix is 0 not computing")
    df_RCBC_return <- NULL
  } else {
    # Test if Year distances and betaNTI are correlated
    euc_temp_mantel <- ecodist::mantel(as.dist(RCBC_var) ~ euc_temp, mrank = T)
    
    df_RCBC_return <- data.frame(Habitat__ = habitat, DepthLumping = depth,
                                 temperature_metric = temperature_metric,
                                 MantelX = c("EucTemp"),
                                 MantelY = "RCBC",
                                 bind_rows(euc_temp_mantel))
  }
  
  
  df_return <- bind_rows(df_bnti_return, df_RCBC_return)
  return(df_return)
}

# Create mantels with temperature 
# Note:
# mantelr	
# pval1	
# one-tailed p-value (null hypothesis: r <= 0).
# 
# pval2	
# one-tailed p-value (null hypothesis: r >= 0).
# 
# pval3	
# two-tailed p-value (null hypothesis: r = 0).
# 
# llim	
# lower confidence limit.
# 
# ulim	
# upper confidence limit.
# eukyear = distances continuous over time (directional)
#binyear = distances year is same as comparison or different. (non-directional)
mantel_with_temperature <- input_ra$sample_metadata %>% 
  select(Habitat__, DepthLumping, T_air.deg_C, mean_AirTemperature_28d, T_soil.deg_C) %>%
  pivot_longer(cols = all_of(c("T_air.deg_C", "mean_AirTemperature_28d", "T_soil.deg_C")),
               names_to = "TemperatureMetric", values_to = "Temperature") %>%
  select(-Temperature) %>%
  distinct() %>%
  rename(habitat = Habitat__, depth = DepthLumping, temperature_metric = TemperatureMetric) %>%
  pmap(., test_temperature_correlation) %>%
  reduce(bind_rows) %>%
  #  filter(MantelX != "BinYear") %>%
  mutate(padj = p.adjust(pval3, method = "fdr"),
         sig = ifelse(padj < 0.05, "significant", "not significant"),
         sig_bin = ifelse(padj < 0.05, 1, 0)) %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))

mantel_with_time_plotdf <- mantel_with_time %>%
  filter(MantelX != "BinYear") %>%
  select(Habitat__, DepthLumping, MantelY, mantelr, sig_bin) %>%
  pivot_longer(matches(c("mantelr", "sig_bin")), names_to = "MeasureType", values_to = "value") %>%
  mutate(StochDet = ifelse(MantelY == "BetaNTI", "Deterministic", "Stochastic"),
         y = ifelse(MeasureType == "sig_bin", 1, 2)) %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))
value_labels <- mantel_with_time_plotdf %>%
  distinct(MeasureType, y) %>%
  mutate(MeasureType = ifelse(MeasureType == "mantelr", "Mantel Corr",
                              "Significance\n(<0.05 = blue)")) %>%
  mutate(Habitat__ = "Palsa", DepthLumping = "0-9", StochDet = "Deterministic") %>%
  mutate(HabitatDepth = paste0(Habitat__, " ", DepthLumping),
         Habitat__ = factor(Habitat__, levels = habitat_levels),
         DepthLumping = factor(DepthLumping, levels = depth_levels))

mantel_with_time.plot <- mantel_with_time_plotdf %>%
  ggplot(aes(x = StochDet, y = y)) +
  geom_tile(aes(fill = value)) +
  #geom_tile(aes(fill = mantelr)) +
  geom_text(data = value_labels, aes(label = MeasureType), size = rel(3), color = "grey20") +
  scale_fill_gradient2(name = "Mantel Correlation", 
                       limits = c(-1, 1)) +
  facet_grid(DepthLumping ~ Habitat__, switch = "y") +
  xlab("Stochastic Or Deterministic Correlation with Time") +
  ylab("") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank() 
  )

ggsave(mantel_with_time.plot, 
       filename = paste0(figures.fp, "/mantel_with_time.png"),
       width = 10, height = 10, dpi = 300)

# Write table of mantels out for sharing data
write.table(mantel_with_time, 
            file = paste0(outputs.fp, "/mantel_with_time.txt"),
            sep = "\t", row.names = F)



#### ====================================================================== ####

# Change with time BetaNTI
#### ====================================================================== ####

test <- betaNTI %>%
  as.matrix() %>%
  as_tibble(rownames = "SampleID") %>%
  pivot_longer(-SampleID) %>%
  #filter(SampleID < name) %>%
  left_join(input_ra$sample_metadata, by = c("SampleID" = "temporal_sample_id")) %>%
  left_join(input_ra$sample_metadata, by = c("name" = "temporal_sample_id")) %>%
  select(SampleID, name, value, starts_with("Habitat__"), starts_with("DepthLumping"), starts_with("Year__")) %>%
  mutate(Year_diff = abs(Year__.x - Year__.y),
         MaxYear = ifelse(Year__.y > Year__.x, Year__.y, Year__.x)) %>%
  filter(Habitat__.x == Habitat__.y) %>% 
  filter(DepthLumping.x == DepthLumping.y) %>%
  filter(Habitat__.x == "Palsa") %>%
#  select(starts_with("Year")) %>% 
  #distinct() %>% 
  filter(DepthLumping.x == "30-39") %>%
  View()


Time_change_assembly <- betanull.lf %>% 
  select(Assembly_Process, BetaNTI, RCBC.nona, 
         Habitat__.Site1, Habitat__.Site2, DepthLumping.Site1, DepthLumping.Site2,
         Year__.Site1, Year__.Site2) %>%
  mutate(Year_diff = abs(Year__.Site1 - Year__.Site2),
         MaxYear = ifelse(Year__.Site2 > Year__.Site1, Year__.Site2, Year__.Site1)) %>%
  # Only select comparisons within a HabitatDepth combination
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
  # Only select neighboring years
  filter(Year_diff == 1) %>%
  mutate(HabDepth = paste0(substr(Habitat__.Site1, 1,1), " ", DepthLumping.Site1), 
         HabDepth = factor(HabDepth, levels = habdepth_levels))

intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              AssemblyMetric = c(rep("BetaNTI", 2),
                                                 rep("RCBC", 2)))

Time_change_assembly %>%
  ggplot(aes(x = MaxYear, y = BetaNTI, group = HabDepth)) +
#  geom_line(aes(color = Habitat__.Site1)) +
  geom_point(aes(color = Habitat__.Site1), shape = 21, size = 2) +
  geom_hline(data = intercept_guide %>% filter(AssemblyMetric == "BetaNTI"),
             aes(yintercept = Distance),
             linetype = "dashed", color = "black") +
  facet_wrap(~HabDepth, nrow = 4, strip.position = "right") +
  geom_smooth(color = "black", linewidth = 2, se = F, alpha = 0.3) +
  theme_bw()




#### ====================================================================== ####

# Plot betanull with alpha diversity
#### ====================================================================== ####
intercept_guide <- data.frame(Distance = c(2, -2, 0.95, -0.95),
                              AssemblyMetric = c(rep("BetaNTI", 2),
                                                  rep("RCBC", 2)))
div_betanull.df <- betanull.lf.diff %>%
  filter(Depth_diff < 10) %>%
  select(Site1, Site2, Habitat_comp, Assembly_Process, BetaNTI, RCBC, AvgRichness,
         AvgShannon) %>%
  group_by(Habitat_comp, Assembly_Process) %>%
  add_tally() %>%
  mutate(Percent = 100* n/sum(n),
         Assembly_Process = factor(Assembly_Process, levels = assembly_levels),
         #Habitat_comp = gsub("(.{1,})(_)(.{1,})", "\\1 vs. \\2", Habitat_comp),
         Habitat_comp = factor(Habitat_comp, levels = c("Palsa_Palsa",
                                                        "Palsa_Bog", "Bog_Bog",
                                                        "Bog_Fen", "Fen_Fen", "Palsa_Fen"))) %>%  
  filter(Habitat_comp %in% c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")) %>%
  pivot_longer(contains("Avg"), names_to = "DiversityMetric", values_to = "Diversity") %>%
  pivot_longer(all_of(c("BetaNTI", "RCBC")), names_to = "AssemblyMetric", values_to = "DistanceMeasure")



div_betanull.df %>%
  mutate(DistanceMeasureBin = ifelse(AssemblyMetric == "RCBC" & Assembly_Process == "Homogenizing dispersal",
                                     0, 
                              ifelse(AssemblyMetric == "RCBC" & Assembly_Process == "Dispersal limitation and drift",
                                     1, NA))) %>%
  group_by(Habitat_comp, AssemblyMetric, DiversityMetric) %>%
  nest() %>%
  #pluck(4,2) %>%
  # lm(DistanceMeasure ~ Diversity, data = .) %>% summary()
  # Running LMs
  mutate(lm.model = purrr::map(data, ~lm(DistanceMeasure~Diversity, data = .)),
         lm.tidied = purrr::map(lm.model, broom::tidy),
         lm.summaried = purrr::map(lm.model, summary),
         lm.model.rsqr = map_dbl(lm.model, ~broom::glance(.)$r.squared),
         lm.model.pval = map_dbl(lm.model, ~broom::glance(.)$p.value)) %>%
  unnest(lm.tidied) #%>% View()
  # Running logistic regression
#   unnest(lm.tidied, names_sep = ".") %>%
#   # mutate(cor.model = map(data, ~cor.test(~CH4_Flux + mean_AirTemp, data = .)),
#   #        cor.tidied = map(cor.model, tidy)) %>%
#   # unnest(cor.tidied, names_sep = ".") %>% 
#   unnest(data)
#   View()

# Notes:
#' Diversity is a poor, albeit significant predictor of alpha diversity. Although
#' All models showed a significant linear relationship between alpha diversity 
#' and assembly metrics, the rsquared was quite low ranging from 1.2e-2  at the lowest
#' significant model ("Bog/Shannon/BetaNTI) to 0.25 (or a 25% of variation explained) for highest significant
#' model  (Palsa/Shannon/RCBC). We note that for the RCBC models, linear models are not really appropriate,
#' given the bimodal nature of the data. See below for logistic regression results for these.
#' 
#' Slopes were not always positive within each habitat, although the general trend across 
#' all habitats was to show a somewhat positive relationship in betaNTI at least. 
#' The strongest positive and negative slopes were for shannon diversity and betaNTI,
#' and found in the Fen (slope = 2.18 H'/BetaNTI) and Palsa (-0.76 H'/BetaNTI), respectively.
#relationship betweeen diversity and assembly metrics NOT significant in Bog/Richness/BetaNTI

div_betanull.df2 <- div_betanull.df %>%
  mutate(DistanceMeasure = ifelse(AssemblyMetric == "RCBC" & DistanceMeasure < -0.95, 0, 
                           ifelse(AssemblyMetric == "RCBC" & DistanceMeasure >0.95, 1, 
                           ifelse(AssemblyMetric == "BetaNTI", DistanceMeasure, NA )))) %>%
  mutate(AssemblyMetric = ifelse(AssemblyMetric == "RCBC", "RCBCBin", "BetaNTI"))

bnti_rcbc_alpha_div <- div_betanull.df %>%
  # Manually add column for significance
  #mutate(linsig = )
  ggplot(aes(x = Diversity, y = DistanceMeasure)) +
  geom_point(aes(fill = Assembly_Process), shape = 21, size = 2, color = "white") +
  facet_wrap(AssemblyMetric ~ DiversityMetric, scales = "free") +
  # geom_hline(data = intercept_guide,
  #            aes(yintercept = Distance),
  #            linetype = "dashed", color = "black") + 
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  # scale_color_manual(name = "Assembly Process", values = colour_assembly, 
  #                    breaks = assembly_levels, 
  #                    labels = assembly_labels) +
  stat_smooth(method = "lm", color = "black") +
  stat_smooth(data = div_betanull.df, 
              aes(color = Habitat_comp, group = Habitat_comp),
              method = "lm") +
  stat_ellipse(aes(color = Habitat_comp, group = Habitat_comp), 
               type = "norm", linetype = 2) +
  scale_color_manual(name = "Habitat", values = colour_habitat,
                     breaks = c("Palsa_Palsa", "Bog_Bog", "Fen_Fen"),
                     labels = habitat_levels) +
  xlab("Average alpha diversity between samples") + ylab("Assembly Distance Measure") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        plot.margin = margin(2,2,2,2, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.placement = "outside",
        strip.text.x = element_text(size = rel(1)),
        strip.text.y = element_text(face = "bold", size = rel(1)),
        legend.position = "bottom")
ggsave(bnti_rcbc_alpha_div, 
       filename = paste0(figures.fp, "/Diversity_linear_assembly.png"),
       width = 10, height = 5, dpi = 400)

# Multinomial version of this
# Prepare df for multinomial calcs
div_mn_pred.prep <- betanull.lf.diff %>%
  filter(Depth_diff < 10) %>%
  select(Site1, Site2, Habitat_comp, Assembly_Process, AvgShannon, AvgRichness) %>%
  filter(Habitat_comp %in% c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels),
         Assembly_Process = fct_rev(Assembly_Process), # reverse the levels so drift is the reference
         Habitat_comp = factor(Habitat_comp, levels = c("Palsa_Palsa",
                                                        "Palsa_Bog", "Bog_Bog",
                                                        "Bog_Fen", "Fen_Fen", "Palsa_Fen"))) %>%
  filter(Habitat_comp %in% c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")) %>%
  mutate(Habitat_comp = case_when(
    Habitat_comp == "Palsa_Palsa" ~ "Palsa",
    Habitat_comp == "Bog_Bog" ~ "Bog",
    Habitat_comp == "Fen_Fen" ~ "Fen",
    TRUE ~ NA_character_
    ),
    Habitat_comp = factor(Habitat_comp, levels = habitat_levels)) %>%
  pivot_longer(contains("Avg"), names_to = "DiversityMetric", values_to = "Diversity") %>%
  mutate(DiversityMetric = fct_recode(DiversityMetric, "Shannon" = "AvgShannon", 
                                   "Richness"="AvgRichness"))

# Version 1 with habitats separate
div_mn_pred <- div_mn_pred.prep %>%
  group_by(Habitat_comp, DiversityMetric) %>%
  nest() %>%
  # Running MNMs
  mutate(data = purrr::map(data, ~ mutate(., Assembly_Process = droplevels(Assembly_Process)))) %>%
  mutate(mn.model = purrr::map(data, ~ nnet::multinom(Assembly_Process~Diversity, data = ., Hess = T)), # run multinomial model
         mn.tidied = purrr::map(mn.model, ~ broom::tidy(., exponentiate = T)), # extract results
         mn.summaried = purrr::map(mn.model, ~ summary(., Wald.ratios = T)), # extract results
         mn.walds = purrr::map(mn.summaried, ~ .[["Wald.ratios"]]),
         mn.pvals =  purrr::map(mn.summaried, ~ data.frame(2 * pnorm(abs(.$Wald.ratios), lower.tail = FALSE)) %>% rownames_to_column())) %>% # get mapped p-values
  mutate(effterm = purrr::map_chr(data, ~paste0("Diversity[",floor(min(.$Diversity)),":",ceiling(max(.$Diversity)), ",by=0.05]")),
           eff = purrr::map(mn.model, ~ggeffects::ggeffect(.,terms = effterm)))

# Plot Multinomial
# extract fitted values
div_mn_pred.fitted <- div_mn_pred %>%
  select(Habitat_comp, DiversityMetric, data, mn.model) %>%
  mutate(mn.fitted = purrr::map(mn.model, ~data.frame(fitted(.)))) %>% 
  select(Habitat_comp, DiversityMetric, data, mn.fitted) %>% #pluck("data", 4)
  unnest(c(data, mn.fitted))

# Plot
div_mn_pred.plot <- div_mn_pred %>%
  unnest(eff)  %>%
  left_join(div_mn_pred %>% 
              select(Habitat_comp, DiversityMetric, mn.pvals) %>% distinct() %>%
              unnest(mn.pvals) %>% distinct() %>%
              rename(DivPval = Diversity) %>%
              mutate(sig = ifelse(DivPval < 0.05, T, F),
                     rowname = gsub(" ", ".", rowname)), 
            by = c("Habitat_comp", "DiversityMetric", "response.level" = "rowname")) %>% 
ggplot(aes(x = x, y = predicted, 
                fill = response.level, color = response.level)) +
  geom_line(aes(linetype = sig)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 1/3) +
  # Fitted values for points
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Drift), color = colour_assembly[5], 
             fill = colour_assembly[5], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Dispersal.limitation.and.drift), color = colour_assembly[4], 
             fill = colour_assembly[4], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Homogenizing.dispersal), color = colour_assembly[3], 
             fill = colour_assembly[3], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Heterogenous.selection), color = colour_assembly[2], 
             fill = colour_assembly[2], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Homogenous.selection), color = colour_assembly[1], 
             fill = colour_assembly[1], 
             alpha = 0.3, shape = 21) +
  # Scale Aesthetics
  scale_linetype_manual(name = "Significantly diff from 0",
                          breaks = c(TRUE, FALSE),
                          values = c("solid", "dashed")) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = gsub(" ", ".", assembly_levels), 
                     labels = assembly_labels) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = gsub(" ", ".", assembly_levels), 
                     labels = assembly_labels) +
  facet_wrap(DiversityMetric ~ Habitat_comp, scales = "free_x") +
  labs(x = 'Average Alpha Diversity of Comparison', y = 'Modeled Probability') +
  theme_bw()
div_mn_pred.plot

ggsave(div_mn_pred.plot, 
       filename = paste0(figures.fp, "/Diversity_multinom_assembly_by_habitat.png"),
       width = 10, height = 5, dpi = 400)

# All habitats together
div_mn_pred_all <- div_mn_pred.prep %>% 
  select(-Site1, -Site2, -Habitat_comp) %>%
  group_by(DiversityMetric) %>%
  nest() %>%
  # Running MNMs
  mutate(data = purrr::map(data, ~ mutate(., Assembly_Process = droplevels(Assembly_Process)))) %>%
  mutate(mn.model = purrr::map(data, ~ nnet::multinom(Assembly_Process~Diversity, data = ., Hess = T)), # run multinomial model
         mn.tidied = purrr::map(mn.model, ~ broom::tidy(., exponentiate = T)), # extract results
         mn.summaried = purrr::map(mn.model, ~ summary(., Wald.ratios = T)), # extract results
         mn.walds = purrr::map(mn.summaried, ~ .[["Wald.ratios"]]),
         mn.pvals =  purrr::map(mn.summaried, ~ data.frame(2 * pnorm(abs(.$Wald.ratios), lower.tail = FALSE)) %>% rownames_to_column())) %>% # get mapped p-values
  mutate(effterm = purrr::map_chr(data, ~paste0("Diversity[",floor(min(.$Diversity)),":",ceiling(max(.$Diversity)), ",by=0.05]")),
         eff = purrr::map(mn.model, ~ggeffects::ggeffect(.,terms = effterm)))

div_mn_pred_all.fitted <- div_mn_pred_all %>%
  select(DiversityMetric, data, mn.model) %>%
  mutate(mn.fitted = purrr::map(mn.model, ~data.frame(fitted(.)))) %>% 
  select(DiversityMetric, data, mn.fitted) %>% #pluck("data", 4)
  unnest(c(data, mn.fitted))

div_mn_pred_all.plot <- div_mn_pred_all %>%
  select(DiversityMetric, eff) %>%
  unnest(eff) %>%
  left_join(div_mn_pred_all %>%
              select(DiversityMetric, mn.pvals) %>% distinct() %>%
              unnest(mn.pvals) %>% distinct() %>%
              rename(DivPval = Diversity) %>%
              mutate(sig = ifelse(DivPval < 0.05, T, F),
                     rowname = gsub(" ", ".", rowname)) %>%
              select(DiversityMetric, rowname, sig),
            by = c("DiversityMetric", "response.level" = "rowname")) %>%
  ggplot(aes(x = x, y = predicted, 
             fill = response.level, color = response.level)) +
  geom_line(aes(linetype = sig)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 1/3) +
  # Fitted values for points
  geom_point(data = div_mn_pred_all.fitted,
             aes(x = Diversity, y = Drift), color = colour_assembly[5], 
             fill = colour_assembly[5], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred_all.fitted,
             aes(x = Diversity, y = Dispersal.limitation.and.drift), color = colour_assembly[4], 
             fill = colour_assembly[4], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred_all.fitted,
             aes(x = Diversity, y = Homogenizing.dispersal), color = colour_assembly[3], 
             fill = colour_assembly[3], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred_all.fitted,
             aes(x = Diversity, y = Heterogenous.selection), color = colour_assembly[2], 
             fill = colour_assembly[2], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred_all.fitted,
             aes(x = Diversity, y = Homogenous.selection), color = colour_assembly[1], 
             fill = colour_assembly[1], 
             alpha = 0.3, shape = 21) +
  # Scale Aesthetics
  scale_linetype_manual(name = "Significantly diff from 0",
                        breaks = c(TRUE, FALSE),
                        values = c("solid", "dashed")) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = gsub(" ", ".", assembly_levels), 
                     labels = assembly_labels) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = gsub(" ", ".", assembly_levels), 
                    labels = assembly_labels) +
  facet_wrap(DiversityMetric ~ ., scales = "free_x") +
  labs(x = 'Average Alpha Diversity of Comparison', y = 'Modeled Probability') +
  theme_bw()
div_mn_pred_all.plot
ggsave(div_mn_pred_all.plot, 
       filename = paste0(figures.fp, "/Diversity_multinom_assembly_allhabs.png"),
       width = 10, height = 5, dpi = 400)
#### ====================================================================== ####




# See if there's any pattern with Dylan's bray-curtis diff data
#### ====================================================================== ####
test <- bc_dist %>%
  # correct the names in Var1 and Var2 to match temporal_sample_id
  mutate(across(all_of(c("Var1", "Var2")), ~gsub("^X", "", .))) %>%
  mutate(across(all_of(c("Var1", "Var2")), ~gsub("[.]", "-", .))) %>%
  full_join(betanull.lf, by = c("Var1" = "Site1", "Var2" = "Site2"))
  

bc.mat <- vegdist(t(input_ra$otu_table[,-1]), method = "bray")

bc.mat_long <- bc.mat %>% as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site1") %>%
  pivot_longer(-Site1, names_to = "Site2", values_to = "value") 

bc_assemb <- inner_join(betanull.lf, bc.mat_long, by = c("Site1", "Site2")) %>%
  select(Site1, Site2, Assembly_Process, Habitat__.Site1, Habitat__.Site2, Year__.Site1, Year__.Site2, 
         DepthLumping.Site1, DepthLumping.Site2, BetaNTI, RCBC,RCBC.nona,
         value) %>% 
  # Filter for matching habitats
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
  filter(Site1 != Site2) %>%
  mutate(year_comparison = ifelse(Year__.Site1 == Year__.Site2, "within a year", "between years"))


ggplot(bc_assemb, aes(y = value, x = year_comparison)) +
  geom_boxplot() +
  geom_point(position = "jitter", aes(fill = Assembly_Process), shape = 21, size = 2) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  facet_grid(Habitat__.Site1 ~ Assembly_Process) +
  ylab("Bray-Curtis Distance") + xlab("") +
  theme_bw()



ggplot(bc_assemb, aes(y = value, x = BetaNTI)) +
  geom_point(aes(fill = year_comparison), shape = 21, size = 2) +
  # scale_color_manual(name = "Assembly Process", values = colour_assembly, 
  #                    breaks = assembly_levels, 
  #                    labels = assembly_labels) +
  # scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
  #                   breaks = assembly_levels, 
  #                   labels = assembly_labels) +
  facet_grid(~Habitat__.Site1) +
  theme_bw()


ggplot(bc_assemb, aes(y = value, x = RCBC)) +
  geom_point(aes(fill = year_comparison), shape = 21, size = 2) +
  # scale_color_manual(name = "Assembly Process", values = colour_assembly, 
  #                    breaks = assembly_levels, 
  #                    labels = assembly_labels) +
  # scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
  #                   breaks = assembly_levels, 
  #                   labels = assembly_labels) +
  facet_grid(~Habitat__.Site1) +
  theme_bw()




#betanull.lf

#### ====================================================================== ####
# multinomial and bray-curtis
#### ====================================================================== ####
# Multinomial version of this
# Prepare df for multinomial calcs
div_mn_pred.prep <- inner_join(betanull.lf, bc.mat_long, by = c("Site1", "Site2")) %>%
  select(Site1, Site2, Assembly_Process, Habitat__.Site1, Habitat__.Site2, Year__.Site1, Year__.Site2, 
         DepthLumping.Site1, DepthLumping.Site2, BetaNTI, RCBC,RCBC.nona,
         value) %>% 
  rename(BrayCurtis = value) %>%
  # Filter for matching habitats
  filter(Habitat__.Site1 == Habitat__.Site2) %>%
  filter(DepthLumping.Site1 == DepthLumping.Site2) %>%
  filter(Site1 != Site2) %>%
  select(Site1, Site2, Habitat__.Site1, Assembly_Process, BrayCurtis) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels),
         Assembly_Process = fct_rev(Assembly_Process), # reverse the levels so drift is the reference
         Habitat = factor(Habitat__.Site1, levels = c("Palsa", "Bog",
                                                        "Fen"))) %>%
  pivot_longer(contains("Avg"), names_to = "DiversityMetric", values_to = "Diversity") %>%
  mutate(DiversityMetric = fct_recode(DiversityMetric, "Shannon" = "AvgShannon", 
                                      "Richness"="AvgRichness"))

# Version 1 with habitats separate
div_mn_pred <- div_mn_pred.prep %>%
  group_by(Habitat_comp, DiversityMetric) %>%
  nest() %>%
  # Running MNMs
  mutate(data = purrr::map(data, ~ mutate(., Assembly_Process = droplevels(Assembly_Process)))) %>%
  mutate(mn.model = purrr::map(data, ~ nnet::multinom(Assembly_Process~Diversity, data = ., Hess = T)), # run multinomial model
         mn.tidied = purrr::map(mn.model, ~ broom::tidy(., exponentiate = T)), # extract results
         mn.summaried = purrr::map(mn.model, ~ summary(., Wald.ratios = T)), # extract results
         mn.walds = purrr::map(mn.summaried, ~ .[["Wald.ratios"]]),
         mn.pvals =  purrr::map(mn.summaried, ~ data.frame(2 * pnorm(abs(.$Wald.ratios), lower.tail = FALSE)) %>% rownames_to_column())) %>% # get mapped p-values
  mutate(effterm = purrr::map_chr(data, ~paste0("Diversity[",floor(min(.$Diversity)),":",ceiling(max(.$Diversity)), ",by=0.05]")),
         eff = purrr::map(mn.model, ~ggeffects::ggeffect(.,terms = effterm)))

# Plot Multinomial
# extract fitted values
div_mn_pred.fitted <- div_mn_pred %>%
  select(Habitat_comp, DiversityMetric, data, mn.model) %>%
  mutate(mn.fitted = purrr::map(mn.model, ~data.frame(fitted(.)))) %>% 
  select(Habitat_comp, DiversityMetric, data, mn.fitted) %>% #pluck("data", 4)
  unnest(c(data, mn.fitted))

# Plot
div_mn_pred.plot <- div_mn_pred %>%
  unnest(eff)  %>%
  left_join(div_mn_pred %>% 
              select(Habitat_comp, DiversityMetric, mn.pvals) %>% distinct() %>%
              unnest(mn.pvals) %>% distinct() %>%
              rename(DivPval = Diversity) %>%
              mutate(sig = ifelse(DivPval < 0.05, T, F),
                     rowname = gsub(" ", ".", rowname)), 
            by = c("Habitat_comp", "DiversityMetric", "response.level" = "rowname")) %>% 
  ggplot(aes(x = x, y = predicted, 
             fill = response.level, color = response.level)) +
  geom_line(aes(linetype = sig)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 1/3) +
  # Fitted values for points
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Drift), color = colour_assembly[5], 
             fill = colour_assembly[5], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Dispersal.limitation.and.drift), color = colour_assembly[4], 
             fill = colour_assembly[4], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Homogenizing.dispersal), color = colour_assembly[3], 
             fill = colour_assembly[3], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Heterogenous.selection), color = colour_assembly[2], 
             fill = colour_assembly[2], 
             alpha = 0.3, shape = 21) +
  geom_point(data = div_mn_pred.fitted,
             aes(x = Diversity, y = Homogenous.selection), color = colour_assembly[1], 
             fill = colour_assembly[1], 
             alpha = 0.3, shape = 21) +
  # Scale Aesthetics
  scale_linetype_manual(name = "Significantly diff from 0",
                        breaks = c(TRUE, FALSE),
                        values = c("solid", "dashed")) +
  scale_color_manual(name = "Assembly Process", values = colour_assembly, 
                     breaks = gsub(" ", ".", assembly_levels), 
                     labels = assembly_labels) +
  scale_fill_manual(name = "Assembly Process", values = colour_assembly, 
                    breaks = gsub(" ", ".", assembly_levels), 
                    labels = assembly_labels) +
  facet_wrap(DiversityMetric ~ Habitat_comp, scales = "free_x") +
  labs(x = 'Average Alpha Diversity of Comparison', y = 'Modeled Probability') +
  theme_bw()
div_mn_pred.plot

ggsave(div_mn_pred.plot, 
       filename = paste0(figures.fp, "/Diversity_multinom_assembly_by_habitat.png"),
       width = 10, height = 5, dpi = 400)


# # Plot Null and empirical distributions
#### ====================================================================== ####
# dist_plot.df <- map_df(null_bMNTD, ~{as.data.frame(.x)}, .id = "null_id") %>% 
#   mutate(null_id = paste0("null_", null_id)) %>% 
#   bind_rows(.,beta.mntd.weighted %>%
#               as_tibble() %>%
#               mutate(null_id = "empirical")) %>%
#   #group_by(null_id) %>%
#   #nest()
#   mutate(null_id = as.factor(null_id)) %>%
#   pivot_longer(cols = !matches("null_id"),
#                names_to = "rep", values_to = "betaMNTD")
# 
# dist_plot <- dist_plot.df %>%
#   filter(null_id != "empirical") %>%
#   ggplot(aes(x = betaMNTD, fill = null_id)) +
#   geom_density(alpha = 0.3) +
#   scale_fill_grey() +
#   geom_density(data = dist_plot.df %>%
#                  filter(null_id == "empirical"),
#                aes(x = betaMNTD),
#                fill = "red", alpha = 0.3) +
#   theme(legend.position = "none")
# 
#### ====================================================================== ####


#### ====================================================================== ####
betanull.lf %>%
  select(Site1, Site2, Assembly_Process,BetaNTI, RCBC.nona,
         Habitat__.Site1, Habitat__.Site2,
         Year__.Site1, Year__.Site2,
         DepthAvg__.Site1, DepthAvg__.Site2,
         DepthLumping.Site1, DepthLumping.Site2) %>%
  select(Site1, Site2, BetaNTI) %>%
  arrange(Site1, Site2) 
  pivot_wider(names_from = "Site2", values_from = BetaNTI) %>%
  column_to_rownames("Site1") %>% as.dist()

#### ====================================================================== ####


#### Save Data and Figures
#### ====================================================================== ####

# Data
# BetaNTI long format
saveRDS(betaNTI.lf, paste0(outputs.fp, "/betaNTI.lf.RDS"))
