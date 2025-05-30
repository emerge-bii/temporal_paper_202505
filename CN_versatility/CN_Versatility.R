# CN Versatility igure

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("CN_versatility", "outputs")
figures.fp <- here("CN_versatility", "figures")

output_dir <- figures.fp
if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}


#### ====================================================================== ####
# Set colors
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(22, 21, 24)
fill_habitat   <- colour_habitat

depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth

year_levels <- c(2011, 2012, 2013, 2014, 2015, 2016, 2017)
year_labels <- c(1, 2, 3, 4, 5, 6, 7)
fill_year   <- RColorBrewer::brewer.pal(7, "BuPu")
colour_year <- fill_year

special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
special_fill <- special_colour

n_special_levels <- c("nitrite_oxidiser", "nitrogen_fixer", "nitrate_reducer", "denitrifier", NA)
n_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-5:-8]
n_special_fill <- n_special_colour


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
colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white")
fill_phylum <- colour_phylum

set.seed(42)
#### ====================================================================== ####
#### Read in data
#### ====================================================================== ####
calculate_cn_versatility <- function(product_refined, genome_clusters) {
  carbon_metabolism <- function(pathways) {
    case_when(
      grepl("cazymes", pathways) ~ "carbon",
      grepl("fermentation", pathways) ~ "carbon",
      grepl("C fixation", pathways) ~ "carbon",
      grepl("methanotrophy", pathways) ~ "carbon",
      grepl("methanogenesis", pathways) ~ "carbon",
      grepl("degradation", pathways) ~ "carbon",
      TRUE ~ NA_character_
    )
  }
  
  nitrogen_metabolism <- function(pathways) {
    case_when(
      grepl("nitrogen", pathways) ~ "nitrogen",
      grepl("reductive_glycine", pathways) ~ "nitrogen",
      pathways == "CAZy-Chitin" ~ "nitrogen",
      pathways == "urea_degradation-all" ~ "nitrogen",
      TRUE ~ NA_character_
    )
  }
  
  
  CNVers <- metapathway_groups %>%
    mutate(CVersatility = carbon_metabolism(metapath),
           NVersatility = nitrogen_metabolism(subpathway))

  
  pathways <- product_refined %>%
    filter(genome %in% genome_clusters$representative) %>%
    left_join(CNVers, by = "subpathway") %>%
    filter(call)
  
  # Create CN paths for type 1 figure
  cpaths <- pathways %>% filter(CVersatility == "carbon") %>% select(-NVersatility) %>% rename(Versatility = CVersatility)
  npaths <- pathways %>% filter(NVersatility == "nitrogen") %>% select(-CVersatility) %>% rename(Versatility = NVersatility)
  
  cn_paths <- bind_rows(cpaths, npaths)
  
  # Create CN paths for type 2 figure
  cpaths <- pathways %>% filter(CVersatility == "carbon") %>% select(-NVersatility) %>% mutate(CVersatility = subpathway)
  npaths <- pathways %>% filter(NVersatility == "nitrogen") %>% select(-CVersatility) %>% mutate(NVersatility = subpathway)
  
  cn_paths_wide <- left_join(cpaths, npaths, by = "genome", suffix = c(".C", ".N")) %>%
    mutate(NVersatility = case_when(is.na(call.N) ~ "None Detected",
                                    TRUE ~ NVersatility))

  return(list(cn_paths_long = cn_paths, cn_paths_wide = cn_paths_wide))
}


pathways_CN <- calculate_cn_versatility(product_refined, genome_clusters)


#### ====================================================================== ####
#plot figure V1 (Hannah's idea building off of Gin's - plotting first because easier to visualize)
#### ====================================================================== ####
# heatmap, C and N abilities on the y axis ordered by C and N ability
# each vip on x axis ordered by vip score and other organisms

# Plot vip
plot_vip_specials_typeOne <- function(habitat, environmental_var){
  vip_CN_plot_df <- vip_members %>% 
    filter(`Environmental Variable` == environmental_var) %>%
    filter(`Data Origin` == habitat) %>%
    left_join(pathways_CN$cn_paths_long, by = c("MAG" = "genome")) %>%
    left_join(specialisations, by = c("MAG" = "genome")) %>%
    mutate(specialisations = factor(specialisation, levels = special_levels))
  
  
  p <- ggplot(data = vip_CN_plot_df, aes(x = fct_reorder(MAG, `VIP Score`, .desc = TRUE), y = subpathway)) +
    geom_point(aes(fill = specialisation, shape = nitrogen_specialisation), size = 3) +
    scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                      guide = guide_legend(override.aes = list(size = 3, shape = 21))) +
    scale_shape_manual(breaks = n_special_levels, values = c(22,23,24,9,21), na.value = 21) +
#    guides(fill = guide_legend(overide.aes = list(shape = 21, size = 15))) +
    facet_grid(Versatility~., scales = "free_y", space = "free_y", switch = "y") +
    xlab("VIP") + ylab("Functions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
  ggsave(plot = p, str_c("type1", environmental_var, habitat, "VIP_specials.png", 
               sep = "_"), 
              height = 7, width = 10,
              path = output_dir, dpi = 300)
}

vip_vars <- vip_members$`Environmental Variable` %>% unique()


# All
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "All", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typeOne(habitat = "All", i)
  }
  
}

# Palsa
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Palsa", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typeOne(habitat = "Palsa", i)
  }
  
}

# Bog
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Bog", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typeOne(habitat = "Bog", i)
  }
  
}

# Fen
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Fen", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typeOne(habitat = "Fen", i)
  }
  
}

#### ====================================================================== ####
#plot figure V2 (Gin's idea)
#### ====================================================================== ####

# figure needs to have C abilities on the y axis and N abilities on the y
# each VIP will be represented by a small dot jittered so as not to overlap showing the metabolic potential of a vip for some VIP list
# goal is cluster of points representing all the vips with a function and their C/N overlap in abilities
plot_vip_specials_typetwo <- function(habitat, environmental_var) {
  vip_CN_plot_df <- vip_members %>% 
    filter(`Environmental Variable` == environmental_var) %>%
    filter(`Data Origin` == habitat)
  
  # number of MAGs
  MAG_num <- length(unique(vip_CN_plot_df$MAG))
  
  vip_wide <- pathways_CN$cn_paths_wide %>%
    filter(genome %in% vip_CN_plot_df$MAG) %>%
    left_join(specialisations, by = c("genome")) %>%
    group_by(CVersatility, NVersatility, specialisation) %>%
    add_tally() %>%
    rename(VIPCspec_count = n) %>%
    ungroup() %>% group_by(CVersatility, NVersatility, nitrogen_specialisation) %>%
    add_tally() %>%
    rename(VIPNspec_count = n) %>% ungroup()
  
  
  p <- ggplot(data = vip_wide, aes(y = CVersatility, x = NVersatility)) +
    geom_point(aes(group = genome, fill = specialisation, shape = nitrogen_specialisation), 
               size = 3,
               position = position_dodge(width=0.1*MAG_num)) +
    scale_shape_manual(breaks = n_special_levels, values = c(22,23,24,9,21), na.value = 21) +
    scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                      guide = guide_legend(override.aes = list(size = 3, shape = 21))) +
    # Add the faceting to make the C-functions easier to parse
    facet_grid(metapath.C~., space = "free_y", scales = "free_y") +
    #coord_fixed(0.7) + 
    xlab("Nitrogen Flexibility") + ylab("Carbon Flexibility") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(0.1, "lines"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          panel.grid = element_blank(),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, face = "bold"))
  
  ggsave(plot = p, str_c("type2", environmental_var, habitat,"VIP_specials.png", 
                         sep = "_"), 
         height = 7, width = 9,
         path = output_dir, dpi = 300)
}

# plot_vip_specials_typetwo("alphaC")
# plot_vip_specials_typetwo("d13C_CO2")
# plot_vip_specials_typetwo("N.percent")
# plot_vip_specials_typetwo("CO2.mM")

vip_vars <- vip_members$`Environmental Variable` %>% unique()

# All
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "All", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typetwo(habitat = "All", i)
  }
  
}

# Palsa 
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Palsa", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typetwo(habitat = "Palsa", i)
  }
  
}

# Bog
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Bog", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typetwo(habitat = "Bog", i)
  }
  
}


# Fen
for(i in vip_vars) {
  test <- vip_members %>% filter(`Data Origin` == "Fen", 
                                 `Environmental Variable` == i)
  if(nrow(test > 1)) {
    plot_vip_specials_typetwo(habitat = "Fen", i)
  }
  
}

#### ====================================================================== ####
#plot figure V3 (Sarah, please ignore, just playing with other options)
#### ====================================================================== ####
# Plot vip
plot_vip_specials_all_typeOne <- function(habitat, environmental_var){
  vip_CN_plot_df <- vip_members %>% 
    mutate(VIP = "yes") %>%
    #filter(`Environmental Variable` == environmental_var) %>%
    filter(`Data Origin` == habitat) %>%
    left_join(pathways_CN$cn_paths_long, by = c("MAG" = "genome")) %>%
    left_join(specialisations, by = c("MAG" = "genome")) %>%
    mutate(specialisations = factor(specialisation, levels = special_levels)) #%>%
    # mutate(VIP = ifelse(is.na(`VIP Score`),"no", "yes"),
    #        `VIP Score` = ifelse(VIP=="yes", `VIP Score`, 0))
  
  
  notvip <- vip_CN_plot_df %>%
    filter(VIP == "no")
  
  p <- ggplot(data = vip_CN_plot_df, aes(x = fct_reorder(MAG, `VIP Score`, .desc = TRUE), y = subpathway)) +
    geom_point(aes(fill = specialisation, shape = nitrogen_specialisation), size = 3) +
    #geom_point(data = notvip, aes(shape = nitrogen_specialisation), fill = "grey50", alpha = 0.3, size = 3) +
    scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                      guide = guide_legend(override.aes = list(size = 3, shape = 21))) +
    scale_shape_manual(breaks = n_special_levels, values = c(22,23,24,9,21), na.value = 21) +
    #    guides(fill = guide_legend(overide.aes = list(shape = 21, size = 15))) +
    facet_grid(Versatility~`Environmental Variable`, scales = "free", space = "free", switch = "y") +
    xlab("VIP") + ylab("Functions") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank())
  
  p
  ggsave(plot = p, str_c("type1", environmental_var, habitat, "VIP_specials.png", 
                         sep = "_"), 
         height = 7, width = 10,
         path = output_dir, dpi = 300)
}

#### ====================================================================== ####
#plot figure V4
#### ====================================================================== ####
# Plot vip
plot_vip_specials_all_typeOne <- function(habitat, environmental_var){
  vip_CN_plot_df <- vip_members %>% 
    #filter(`Environmental Variable` == environmental_var) %>%
    filter(`Data Origin` == habitat)
  
  # number of MAGs
  MAG_num <- length(unique(vip_CN_plot_df$MAG))
  
  vip_wide <- pathways_CN$cn_paths_wide %>%
    filter(genome %in% vip_CN_plot_df$MAG) %>%
    left_join(vip_CN_plot_df, by = c("genome" = "MAG")) %>%
    left_join(specialisations, by = c("genome")) %>%
    group_by(`Data Origin`, `Environmental Variable`, CVersatility, NVersatility) %>%
    add_tally() %>%
    rename(VIPCNspec_count = n) %>% ungroup()
  
  
  # vip_wide_facets <- unique(vip_wide[,c("metapath.x")])
  # vip_wide_facets$CVersatility <- max(vip_wide$CVersatility)
  # vip_wide_facets$NVersatility <- max(vip_wide$NVersatility)
  
  p <- ggplot(data = vip_wide, aes(y = CVersatility, x = NVersatility)) +
    geom_point(aes(fill = specialisation, shape = nitrogen_specialisation,
                   size = VIPCNspec_count)) +
    scale_shape_manual(breaks = n_special_levels, values = c(22,23,24,9,21), na.value = 21) +
    scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                      guide = guide_legend(override.aes = list(size = 3, shape = 21))) +
    facet_grid(metapath.C~`Environmental Variable`, space = "free_y", scales = "free_y") +
    #coord_fixed(0.7) + 
    xlab("Nitrogen Flexibility") + ylab("Carbon Flexibility") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(0.1, "lines"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          panel.grid = element_blank(),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, face = "bold"))
  
  p
  ggsave(plot = p, str_c("type2", environmental_var, habitat,"VIP_specials.png", 
                         sep = "_"), 
         height = 7, width = 9,
         path = output_dir, dpi = 300)
}
