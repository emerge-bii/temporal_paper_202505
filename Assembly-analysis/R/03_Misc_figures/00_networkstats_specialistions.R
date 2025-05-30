

library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)
library(gganimate)

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

#### Define plotting colors/levels
#### ====================================================================== ####
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
  "other"
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
  "other"
)
colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white", "grey30")
fill_phylum <- colour_phylum

#### ====================================================================== ####
#### Read in data 
#### ====================================================================== ####
# Read in input
input <- input_counts

#### ====================================================================== ####
# Add Specialisations, taxonomy, and VIP status to network stats
#### ====================================================================== ####
network_stats_spec <- network_stats %>%
  group_by(genome_id, Habitat) %>%
  left_join(specialisations, by = c("genome_id" = "genome")) %>% 
  left_join(input_counts$taxonomy, by = c("genome_id" = "genome")) %>%
  left_join(vip_members %>% select(MAG) %>% distinct() %>% mutate(VIP = TRUE),
            by = c("genome_id" = "MAG")) %>%
  mutate(VIP = ifelse(is.na(VIP), FALSE, VIP))
#### ====================================================================== ####
# Plotting
#### ====================================================================== ####
network_stats %>%
  pivot_longer(contains("centrality"), names_to = "centrality_metric", values_to = "centrality") %>%
  group_by(Habitat, centrality_metric, genome_id) %>%
  summarize(across(all_of("centrality"), list(mean = ~mean(., na.rm = T), 
                                              sd = ~sd(., na.rm = T),
                                              se = ~sd(.,na.rm = T)/sqrt(n())))) %>%
  mutate(lower = centrality_mean - centrality_se,
         upper = centrality_mean + centrality_se,
         Habitat = factor(Habitat, levels = habitat_levels)) %>%
  left_join(specialisations, by = c("genome_id" = "genome")) %>%
  left_join(input_counts$taxonomy, by = c("genome_id" = "genome")) %>%
  mutate(PhylumLabel = ifelse(!is.na(Phylum) & !(Phylum %in% phylum_levels), "other", Phylum)) %>% 
  ggplot(aes(x = centrality_mean, y = specialisation)) + 
#  geom_boxplot() +
  geom_pointrange(aes(xmin = lower,
                      xmax = upper,
                      fill = PhylumLabel), position = "jitter", shape = 21, alpha = 0.5) +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  facet_wrap(Habitat~centrality_metric, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


network_stats %>%
  ggplot(aes(x = eigenvector_centrality)) +
  geom_density(aes(group = genome_id)) +
  facet_wrap(~Habitat)

# Of note: no monomer degraders present in network VIPs
network_stats_spec.plot <- network_stats_spec %>%
  mutate(jitter = runif(1, min = -0.3, max = 0.3),
         spec_fact = ifelse(is.na(specialisation), "Not Assigned",
                               specialisation),
         spec_fact = factor(spec_fact, levels = c(special_levels[1:6], "Not Assigned")),
         numeric_spec = as.numeric(spec_fact) + jitter) # add consistent jitter for smooth plotting
cent_time <- ggplot(network_stats_spec.plot, 
       aes(x = eigenvector_centrality, y = spec_fact)) +
  geom_blank() +
  geom_violin(aes(fill = spec_fact), alpha = 0.3) +
  geom_point(aes(y = numeric_spec, fill = spec_fact), shape = 21, size = 2, alpha = 0.4) +
  scale_fill_manual(name = "Carbon Specialisation",
                    values = c(special_colour[1:6], "grey80"),
                    breaks = c(special_levels[1:6], "Not Assigned"),
                    labels = c(special_levels[1:6], "Not Assigned")) + # omits monomer degraders
  geom_point(data = network_stats_spec.plot %>% filter(VIP == TRUE), 
             aes(y = numeric_spec),
             size = 0.5, color = "black") +
  facet_wrap(Habitat~Year, scale = "free_x", ncol = 7) +
  theme_bw()

ggsave(plot = cent_time, filename = paste0(figures.fp, "/Centrality_Time_Hab_CSpec.png"),
       width = 20, height = 10, dpi = 400)

network_stats_spec %>%
  ungroup() %>% group_by(Habitat) %>%
  nest() %>%
  mutate(KruskalWallis = purrr::map(data, ~rstatix::kruskal_test(eigenvector_centrality ~ specialisation, data = .)),
         Dunn = purrr::map(data, ~rstatix::dunn_test(eigenvector_centrality ~ specialisation, data = ., p.adjust.method = "fdr", detailed = T))) %>% 
  mutate(KW_P = map_dbl(KruskalWallis, ~.$p)) %>%
  unnest(Dunn) %>% 
  filter(p.adj < 0.05) %>%
  filter(if_any(contains("group"), ~grepl("generalist", .))) %>% View()


network_stats_spec <- network_stats %>%
  group_by(genome_id, Habitat) %>%
  left_join(specialisations, by = c("genome_id" = "genome")) %>% 
  left_join(input_counts$taxonomy, by = c("genome_id" = "genome")) %>%
  left_join(vip_members %>% select(MAG) %>% distinct() %>% mutate(VIP = TRUE),
            by = c("genome_id" = "MAG")) %>%
  mutate(VIP = ifelse(is.na(VIP), FALSE, VIP))

ggplot(network_stats_spec, aes(x = eigenvector_centrality, y = specialisation)) +
  geom_violin(aes(fill = specialisation), alpha = 0.3) +
  geom_point(aes(fill = specialisation), position = "jitter", shape = 21, alpha = 0.4) +
  scale_fill_manual(name = "Carbon Specialisation", 
                    values = special_colour[c(1:4,6,8)],
                    breaks = special_levels[c(1:4,6,8)],
                    labels = special_levels[c(1:4,6,8)]) + # omits monomer degraders
  #geom_density(aes(color = specialisation, fill = specialisation), alpha = 0.5) +
  facet_wrap(Habitat~Year, scale = "free_x", ncol = 7) +
  theme_bw()

anim <- ggplot(network_stats_spec %>% #filter(Habitat == "Palsa") %>% 
                 ungroup() %>% group_by(Habitat, Year) %>%
                 slice_max(n = 40, order_by = eigenvector_centrality) %>%
                 mutate(rank = rank(eigenvector_centrality, ties.method = "random")*1.0) %>% # smoother transitions with doubles; have to break ties or bars overlap
                 mutate(YearLabel = ifelse(Year == 2011, "Start: 2011", # make labels for years
                                           ifelse(Year == 2017, "End: 2017",
                                                  ifelse(Year == 2013, "2013*", as.character(Year)))),
                        YearLabel = fct_reorder(YearLabel, Year)) %>%
                 mutate(TaxLabel = paste0(gsub("f__", "", Family), ": ", genome_id),
                        TaxLabel = ifelse(VIP == TRUE, paste0("*", TaxLabel), TaxLabel)), # add asterisk for VIPs 
               aes(y = eigenvector_centrality, 
                   x = rank)) +
  geom_tile(aes(y = eigenvector_centrality/2, # because y demarks the center of tile
                height = eigenvector_centrality,
                fill = specialisation, width = 0.9), alpha = 0.9) +
  scale_fill_manual(name = "Carbon Specialisation", 
                    values = special_colour[c(1:4,6,8)],
                    breaks = special_levels[c(1:4,6,8)],
                    labels = special_levels[c(1:4,6,8)]) + # omits monomer degraders
#  geom_text(aes(y = eigenvector_centrality, label = genome_id), vjust = -0.5) + # labels above bars
  geom_text(aes(y = 0, label = TaxLabel), hjust = 1) + # labels on side
  coord_cartesian(clip = "off", expand = FALSE) + # clip off for bottom labels
  coord_flip(expand = FALSE, clip = "off") +
  facet_wrap(~Habitat, scales = "free", ncol = 1) +
  labs(title = '{closest_state}', x = "", y = "Eigenvector centrality") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 22),
        axis.ticks.y = element_blank(), ## axis.ticks.y shows the ticks on the flipped x-axis (the now metric), and hides the ticks from the geog layer
        axis.text.y = element_blank(),
        plot.margin = margin(10, 10, 10, 450))  +
  transition_states(states = YearLabel, transition_length = 2, state_length = 1) +
  ease_aes('cubic-in-out') +
  enter_fade() +
  exit_fade()
animate(anim, nframes = 420, fps = 20, height = 24, width = 10, units = "in", 
        res = 100, start_pause = 10, end_pause = 10)
anim_save(filename = "~/Downloads/centrals.gif")
animate(anim, nframes = 420, fps = 20, height = 24, width = 10, units = "in", 
        res = 100, start_pause = 10, end_pause = 10, renderer = av_renderer())
anim_save(filename = "~/Downloads/centrals.mp4")


# Animated plot settings:
animate_rank_plot_settings <- list(
  geom_tile(aes(y = eigenvector_centrality/2, # because y demarks the center of tile
                height = eigenvector_centrality, width = 0.9), alpha = 0.9),
  geom_text(aes(y = 0, label = genome_id), hjust = 1), # labels on side
  coord_flip(expand = FALSE, clip = "off"), # clip off for bottom labels
  facet_wrap(~Habitat, scales = "free", ncol = 1),
  labs(title = '{closest_state}', x = "", y = "Eigenvector centrality"),
  theme_bw(),
  theme(plot.title = element_text(hjust = 1, size = 22),
        axis.ticks.y = element_blank(), ## axis.ticks.y shows the ticks on the flipped x-axis (the now metric), and hides the ticks from the geog layer
        axis.text.y = element_blank(),
        plot.margin = margin(10, 10, 10, 350)),
  # animation settings
  transition_states(states = YearLabel, transition_length = 2, state_length = 1),
  ease_aes('cubic-in-out'),
  enter_fade(),
  exit_fade()
)
anim <- ggplot(network_stats_spec %>% #filter(Habitat == "Palsa") %>% 
                 ungroup() %>% group_by(Habitat, Year) %>%
                 slice_max(n = 40, order_by = eigenvector_centrality) %>%
                 mutate(rank = rank(eigenvector_centrality, ties.method = "random")*1.0) %>% # smoother transitions with doubles; have to break ties or bars overlap
                 mutate(YearLabel = ifelse(Year == 2011, "Start: 2011", # make labels for years
                                           ifelse(Year == 2017, "End: 2017",
                                                  ifelse(Year == 2013, "2013*", as.character(Year)))),
                        YearLabel = fct_reorder(YearLabel, Year)) %>%
                 mutate(TaxLabel = paste0(gsub("f__", "", Family), ": ", genome_id),
                        TaxLabel = ifelse(VIP == TRUE, paste0("*", TaxLabel), TaxLabel)), # add asterisk for VIPs 
               aes(y = degree_centrality, 
                   x = rank)) +
  animate_rank_plot_settings +
  ylab("Degree Centrality") +
  scale_fill_manual(name = "Carbon Specialisation", 
                    values = special_colour[c(1:4,6,8)],
                    breaks = special_levels[c(1:4,6,8)],
                    labels = special_levels[c(1:4,6,8)]) # omits monomer degraders
  
animate(anim, nframes = 420, fps = 20, height = 24, width = 10, units = "in", 
        res = 100, start_pause = 10, end_pause = 10, renderer = av_renderer())
anim_save(filename = "~/Downloads/centrals_decent.mp4")


#### ====================================================================== ####

# Centrality over time
#### ====================================================================== ####
cent_time <- network_stats_spec %>% ungroup() %>%
  group_by(Habitat, genome_id) %>%
  mutate(avg_eigenvector_centrality = mean(eigenvector_centrality, na.rm = T),
         genome_rank = rank(avg_eigenvector_centrality, ties.method = "random")*1.0) %>%
  ungroup() %>% group_by(Habitat) %>%
  slice_sample(n = 300) %>%
  ggplot(aes(x = Year, y = eigenvector_centrality)) +
  geom_line(aes(group = genome_id, color = genome_rank)) +
  geom_point() +
  theme(legend.position = "none")
ggsave(plot = cent_time, filename = paste0(figures.fp, "/Centrality_Time_Hab_CSpec.png"),
       width = 10, height = 5, dpi = 400)
#### ====================================================================== ####
