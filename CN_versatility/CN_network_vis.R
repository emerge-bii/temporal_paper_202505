# CN Versatility igure

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library('ggraph')
library('igraph')
library('ggplot2')
library('tidygraph')
library('dplyr')
library('HiveR')
library('grid')
library("cowplot")
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

metapathway_levels <- metapathway_groups %>% `$`(metapath) %>% unique()
colour_metapathway <- RColorBrewer::brewer.pal(10, "Set3")

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
vip_members <- read_tsv(here("data","vip_members_pathways-08-29-2024.txt")) %>%
  select(-Pathways) %>% distinct()

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


# Read in organisms for naming
named_organisms <- read_csv(here("data", "all_organisms_for_naming_v3.csv"))



named_organisms_right_pathway <- named_organisms %>%
  left_join(product_refined, by = "genome") %>%
  filter(call) %>%
  group_by(genome) %>%
  mutate(Pathways = paste(subpathway, collapse = "; ")) %>%
  select(-subpathway, -call) %>%
  distinct()

write_csv(named_organisms_right_pathway, here("data", "all_organisms_for_naming_v4.csv"))

###############################################################################

# Hive plot
# Calculate correlation of each VIP to alphaC
plot_hive <- function(EnvVar = "C.percent", 
                      DatOrig = "Fen", 
                      check_cor = FALSE,
                      y_title = expression(paste("Correlation to ", alpha[c], " across fen and bog samples")),
                      Significant_only = T,
                      mag_text = F # text labels for MAGs
                      ) {
  
  # y_title = title of y axis in hive plot as an expression (useful for special characters)
  vip_specials <- vip_members %>% 
    filter(`Environmental Variable` == EnvVar & `Data Origin` == DatOrig) %>%
    distinct()  # this is because the d15N_peat blue module (All) is duplicated in vip_members for some reason

  # Set up habitat filtering
  dat_orig <- DatOrig
  if(DatOrig == "All") {
    dat_orig <- c("Palsa", "Bog", "Fen")
    } 
  
  env_var_corr <- trimmed_mean$trimmed_mean %>%
    filter(temporal_sample_id %in% sample_metadata$temporal_sample_id) %>% # filters out outliers
    left_join(sample_metadata %>% select(temporal_sample_id, matches(EnvVar), Habitat__), by = "temporal_sample_id") %>%
    filter(if_any(matches(EnvVar), ~ !is.na(.))) %>% # remove NAs for environmental variable
    # filter by habitat
    filter(Habitat__ %in% dat_orig) %>% # remove NAs for environmental variable
    # rename EnvVar in case it doesn't match
    rename(!!EnvVar := 5) %>%
    group_by(genome) %>%
    nest() %>% 
    filter(genome %in% vip_specials$MAG) %>%
    mutate(cor.model = purrr::map(data, ~cor.test(.x[[EnvVar]], .x[["relabund_of_recovered"]], method = "spearman")),
           cor.tidied = purrr::map(cor.model, broom::tidy)) %>%
    unnest(cor.tidied) %>%
    ungroup() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    left_join(vip_specials, by = c("genome" = "MAG")) %>%
    #select(genome, estimate,`VIP Score`, p.value) %>%
    rename(spearman_correlation = estimate) #%>%
    # ungroup() %>%
    # mutate(type = dplyr::case_when(
    #   p.value <= 0.05 ~ 'CorToAlphaC',
    #   `VIP Score` >= 1 ~ 'VIP',
    #   TRUE ~ 'Both'
    # ))
  
  
  # Side task: Investigate if bog/fen samples are having undue impact on results
  if(check_cor) {
    scatter_fun = function(x, y, df, genome) {
      dat <- left_join(df, sample_metadata %>% 
                         select(temporal_sample_id, Habitat__) %>%
                         filter(Habitat__ %in% dat_orig),
                       by = c("temporal_sample_id", "Habitat__"))
      print(genome)
      p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]]) ) +
        geom_point(aes(color = Habitat__)) +
        scale_color_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels) +
        geom_smooth(method = "loess", se = FALSE, color = "grey74") +
        theme_bw() +
        ggtitle(label = genome)
      p
    }
    
    
    cor_check_plots <- env_var_corr %>%
      mutate(plottitle = paste0(genome," (", DatOrig, " - ", EnvVar,"): ", "rho = ", round(spearman_correlation, digits = 2), "; p.adj = ",
                                round(p.adj, digits = 2)))
    plot_list <- purrr::map2(cor_check_plots$data, cor_check_plots$plottitle, ~scatter_fun(x = "relabund_of_recovered", y = EnvVar, df = .x, genome = .y) )
    
    names(plot_list) <- cor_check_plots$genome
    
  }
  
  # Get ready for Hive plot
  axes_data <- tibble(x1 = 0.9, y1 = min(0, min(env_var_corr$spearman_correlation)), spearman_correlation = 1,
                      `VIP Score` = 1.05*max(env_var_corr$`VIP Score`)) 
  
  env_var_corr.plot <- env_var_corr %>%
    mutate(Significant = ifelse(p.adj < 0.05, "yes", "no")) %>%
    mutate(y1 = axes_data$y1, x1 = 0.9) %>%
    left_join(specialisations, by = "genome")
  if(Significant_only) {
    env_var_corr.plot <- env_var_corr.plot %>%
      filter(Significant == "yes")
  }
  
  p <- ggplot(env_var_corr.plot, aes(color = specialisation)) +
    #  expand_limits(x = 0, y = c(0, 9)) +
    scale_x_continuous(breaks = c(0.9, round(1*max(env_var_corr$`VIP Score`), digits = 1)), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    # add dotted line for 0 correlation
    geom_hline(yintercept =  0, color = "grey80", linetype = "dashed") +
    #y axis arrow
    geom_segment(data = axes_data, aes(x = x1, y = y1, yend = spearman_correlation), xend = 0.9, color = "black", 
                 arrow = arrow(length = unit(.15,"inches"), type = "closed", angle = 15)) +
    # x axis arrow
    geom_segment(data = axes_data, aes(y = y1, xend = `VIP Score`, yend = y1),x = 0.9, color = "black", 
                 arrow = arrow(length = unit(.15,"inches"), type = "closed", angle = 15)) +
    # Hive curves and points
    geom_curve(aes(y = spearman_correlation, yend = y1, x = x1, xend = `VIP Score`, alpha = Significant), curvature = -0.3) +
    # y-points and text
    geom_point(aes(y = spearman_correlation, fill = specialisation, alpha = Significant,
                   shape = Significant), x = 0.9, color = "black", size = 2) +
   
    # x-points
    geom_point(aes(x = `VIP Score`, fill = specialisation, alpha = Significant, 
                   shape = Significant),
               y = axes_data$y1, color = "black", size = 2) +
    scale_color_manual(name = "Specialisation", values = special_colour, breaks = special_levels) +
    scale_fill_manual(name = "Specialisation", values = special_colour, breaks = special_levels) +
    scale_alpha_manual(name = "", values = c(0.2, 1), breaks = c("no", "yes")) +
    scale_shape_manual(name = "", values = c(16, 21), breaks = c("no", "yes")) + 
    guides(alpha = "none", shape = "none", color = guide_legend(override.aes = list(shape = 21))) +
    xlab("Relative VIP Score") +
    ylab(y_title) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  if(mag_text) {
    p <- p +
      # y text
      geom_text_repel(aes(y = spearman_correlation, label = genome, alpha = Significant),
                      x = 0.9, color = "black", size = 2, hjust = -0.3, vjust = 0) +
      geom_text_repel(aes(x = `VIP Score`, label = genome, alpha = Significant),  
                      y = axes_data$y1, color = "black", size = 2, hjust = -0.2, vjust = 0, angle = 45)
  }
  
  if(check_cor) {
    
    env_var_corr <- env_var_corr %>% 
      left_join(input_ra$taxonomy, by = "genome")
    
    return_list <- list(plot = p, cor_check_plots = plot_list, data = env_var_corr)
    return(return_list)
  } else {
    return(p)
  }
}


# AlphaC plus CH4 and CO2 enrichment
p <- plot_hive(EnvVar = "d13C_CH4", DatOrig = "All", check_cor = T, Significant_only = F)
#p$plot

pdf(paste0(output_dir, "/d13C_CH4_corr_plots.pdf"))
p$cor_check_plots 
dev.off()

# Bog AlphaC
p <- plot_hive(EnvVar = "alphaC", DatOrig = "Bog", check_cor = T, Significant_only = T)
#p$plot

pdf(paste0(output_dir, "/bog_alphaC_corr_plots.pdf"))
p$cor_check_plots 
dev.off()

# Fen C percent
p <- plot_hive(EnvVar = "C.percent", DatOrig = "Fen", check_cor = T, Significant_only = T)
#p$plot

pdf(paste0(output_dir, "/Fen_C.percent_corr_plots.pdf"))
p$cor_check_plots 
dev.off()


p <- plot_hive(EnvVar = "d13C_CO2", DatOrig = "All", check_cor = T, Significant_only = F)
#p$plot

pdf(paste0(output_dir, "/d13C_C02_corr_plots.pdf"))
p$cor_check_plots
dev.off()

p <- plot_hive(check_cor = T, Significant_only = F)
#p$plot

pdf(paste0(output_dir, "/alphaC_corr_plots.pdf"))
p$cor_check_plots
dev.off()





# Save all hive plot combos:
save_all_hive <- function(wgcna_combos, significant_only = T){
  DatOrig <- wgcna_combos[[1]]
  if(DatOrig == "All") {hab <- "all habitats"} else {hab <- paste(str_to_lower(DatOrig), "samples")}
  Mod_col <- wgcna_combos[[2]]
  EnvVar <- wgcna_combos[[3]]
  
  y_title_list <- list(
    d15N_peat = paste("Correlation to \u03B4\U00B9\u2075", "N across", hab),
    d13C_CO2 = paste("Correlation to \u03B4\U00B9\u00B3", "C CO\u2082 across", hab),
    d13C_CH4 =  paste("Correlation to \u03B4\U00B9\u00B3","C CH\u2084", " across ", hab),
    d13C_peat =  paste("Correlation to \u03B4\U00B9\u00B3","C peat", " across ", hab),
    alphaC =  paste("Correlation to \u0251C", " across ", hab),
    N.percent =  paste("Correlation to total %N across ", hab),
    C.percent =  paste("Correlation to total %C across ", hab),
    CO2.mM =  paste("Correlation to [CO\u2082](mM) across ", hab),
    TN.mM =  paste("Correlation to total N (mM) across ", hab),
    DOC.mM =  paste("Correlation to Dissolved Organic Carbon (mM) across ", hab)
  )
  y_title <- y_title_list[[EnvVar]]
  
  
  hive_combo_plot <- plot_hive(EnvVar = EnvVar, DatOrig = DatOrig, check_cor = T, Significant_only = significant_only,
                               y_title = y_title)
  if(significant_only) {significant <- "Significant"} else {significant <- "all"}
  ggsave(plot = hive_combo_plot$plot, str_c(DatOrig, Mod_col, EnvVar, significant, "VIP_hive_plot.png", 
                         sep = "_"), 
         height = 10, width = 10,
         path = output_dir, dpi = 300)
  
  # write correlations out to table
  write_tsv(hive_combo_plot$data, 
            file = paste0(output_dir, "/", str_c(DatOrig, Mod_col, EnvVar, significant, "VIP_corr_table.txt", 
                                               sep = "_")))
  
  # write_cor_check plots
  pdf(paste0(output_dir, "/", str_c(DatOrig, Mod_col, EnvVar, significant, "corr_check_plots.pdf", 
                                      sep = "_")))
  print(hive_combo_plot$cor_check_plots)
  dev.off()
}


#save_all_hive(wgcna_combos = wgcna_combos[1,])

# Get list of all module combinations
wgcna_combos <- vip_members %>% select(`Data Origin`, `Module Color`, `Environmental Variable`) %>%
  distinct()

# Significant
wgcna_combos %>%
  #slice_head(n = 1) %>%
  split(1:nrow(.)) %>%
  purrr::map(., ~save_all_hive(.))

wgcna_combos %>%
  split(1:nrow(.)) %>%
  purrr::map(., ~save_all_hive(., significant_only = F))

# Read in correlation tables and combine into one file
tsv_files <- fs::dir_ls(output_dir, regexp = "\\VIP_corr_table.txt$")


#old <- combo_VIP_corr_table
combo_VIP_corr_table <- tsv_files %>% 
  map_dfr(read_tsv) %>%
  distinct() %>% 
  select(-cor.model, -data)


write_tsv(combo_VIP_corr_table, file = paste0(output_dir, "/all_vip_env_var_corr.tsv"))


################################################################################
# CN Capability Plots
################################################################################
NVers <- pathways_CN$cn_paths_wide %>%
  select(genome, NVersatility) %>% distinct()

environmental_var <- "alphaC"
habitat <- "All"


# Hannah filter out correlations below 0.5

plot_CN_met <- function(EnvVar = "alphaC", DatOrig = "Bog", MAGs = NULL,
                        Significant_only = F, sigMags = NULL) {
  # MAGs, string of MAGs to plot CN versatilty 
  if(is.null(MAGs)) {
    
    # # If significant only, check that significant mag list provided
    # # Significant only works only for plotting CN of module members
    # if(is.null(sigMags) & Significant_only) {
    #   print(paste0("Error: You want to plot only significant mags,",
    #                " but no list of significant Mags provided.\n",
    #                "Please specify list in sigMags and try again"))
    #   stop()
    # }
    
    if(Significant_only) {
      vip_members_filt <- vip_members %>%
        left_join(combo_VIP_corr_table %>% select(genome:`Environmental Variable`),
                  by = c("MAG" = "genome", "Environmental Variable", "Data Origin", "VIP Score", "Module Color")) %>%
        #filter(MAG %in% sigMags) %>%
        filter(p.adj < 0.05) %>%
        filter(spearman_correlation >= 0.5) %>%
        filter(`Environmental Variable` == EnvVar) %>%
        filter(`Data Origin` == DatOrig)
      
      if(nrow(vip_members_filt) == 0) {
        return(NULL)
      }
    } else {
      vip_members_filt <- vip_members %>%
        left_join(combo_VIP_corr_table %>% select(genome:`Environmental Variable`),
                  by = c("MAG" = "genome", "Environmental Variable", "Data Origin", "VIP Score", "Module Color")) %>%
        filter(`Environmental Variable` == EnvVar) %>%
        filter(`Data Origin` == DatOrig)
    }

    
    mag_CN_plot_df <- vip_members_filt %>%
      mutate(MAG_Rank = rank(abs(spearman_correlation))) %>%
      left_join(input_counts$taxonomy, by = c("MAG" = "genome")) %>%
      mutate(across(Domain:Species, ~gsub("[dpcofgs]__", "", .))) %>%
      mutate(genome_labels = paste0(Phylum, ": ", Genus, " - ", MAG)) %>%
      # Change names of previously named genomes
      mutate(genome_labels = case_when(MAG == "20170700_S25_26" ~ paste0(Phylum, ": ", "Ca. A. stordalenmirensis"),
                                       MAG == "PLGY01" ~ paste0(Phylum, ": ", "Ca. M. stordalenmirensis"),
                                       MAG == "PMEE01" ~ paste0(Phylum, ": ", "Ca. M. crilli"),
                                       MAG == "20120700_S1X_2" ~ paste0(Phylum, ": ", "Ca. Acididegradans cibantis"),
                                       MAG == "20170700_S25_14" ~ paste0(Phylum, ": ", "Ca. Acidimolecula latum"),
                                       MAG == "3300037175_6" ~ paste0(Phylum, ": ", "Ca. Acidiacetogenum stordalenmirensis"),
                                       MAG == "20170700_S25_9" ~ paste0(Phylum, ": ", "Ca. Acidimolecula amplum"), 
                                       MAG == "20110800_E3D_9" ~ paste0(Phylum, ": ", "Ca. Acididegradans virginiarichiae"),
                                       MAG == "20110800_S1D_14" ~ paste0(Phylum, ": ", "Ca. Abiskobacterium ruthvarnerae"),
                                       MAG == "20160700_S15_2" ~ paste0(Phylum, ": ", "Ca. Acidiflorens viviennsis"),
                                       MAG == "20100900_E1D_10" ~ paste0(Phylum, ": ", "Ca. Acidiflorens anabellensis"),
                                       MAG == "3300036929_5" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum latum"),
                                       MAG == "3300037085_4" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum fermentum"),
                                       MAG %in% c("20110800_E2D_15", "20111000_P1M_1","20111000_P2D_6", "20120500_P15_4",
                                                  "20120500_P23_1", "20120500_P24_7", "20120500_P26_8", "20120500_P36_7",
                                                  "20120700_P2D_10", "20120700_P3D_89", "20130700_C2X_13","20140700_E22_12",
                                                  "20150700_S23_3", "20170700_S25_25", "PKWB01", "PLNE01", "PLYD01", "PLYZ01",
                                                  "PLZA01", "PMOL01", "PNAG01", "PNBE01") ~ paste0(Phylum, ": ", "Ca. Changshengia", " spp."),
                                       TRUE ~ genome_labels),
      #mutate(genome_labels = paste0(genome_labels, ": ", MAG)) %>%
             genome_labels = str_wrap(genome_labels, width = 30)) %>%
      left_join(pathways_CN$cn_paths_long, by = c("MAG" = "genome")) %>%
      left_join(specialisations, by = c("MAG" = "genome")) %>%
      mutate(specialisations = factor(specialisation, levels = special_levels),
             MAG = fct_reorder(MAG, MAG_Rank, .desc = FALSE),
             genome_labels = fct_reorder(genome_labels, MAG_Rank, .desc = FALSE)) %>%
      mutate(Versatility = dplyr::case_when(
        subpathway %in% c('reductive_glycine-serine_pathway', 'reductive_glycine-acetyl_pathway', "CAZy-Chitin") ~ 'Both C & N',
        metapath == 'nitrogen' ~ 'Nitrogen',
        TRUE ~ 'Carbon'),
        metapath = ifelse(grepl("Both", Versatility), "both C & N", metapath)
      ) %>%
      mutate(Versatility = factor(Versatility, levels = c("Carbon", "Both C & N", "Nitrogen")),
             metapath = factor(metapath, levels = c("cazymes", "C degradation", "C fermentation", "methanogenesis", "C fixation", "both C & N","nitrogen")))
    
    y_title_list <- list(
      d15N_peat = paste("MAGs (Ranked by correlation to \u03B4\U00B9\u2075", "N)"),
      d13C_CO2 = paste("MAGs (Ranked by correlation to \u03B4\U00B9\u00B3", "C CO\u2082)"),
      d13C_CH4 =  paste("MAGs (Ranked by correlation to \u03B4\U00B9\u00B3","C CH\u2084)"),
      d13C_peat =  paste("MAGs (Ranked by correlation to \u03B4\U00B9\u00B3","C peat)"),
      alphaC =  paste("MAGs (Ranked by correlation to \u0251C)"),
      N.percent =  paste("MAGs (Ranked by correlation to total %N)"),
      C.percent =  paste("MAGs (Ranked by correlation to total %C)"),
      CO2.mM =  paste("MAGs (Ranked by correlation to [CO\u2082](mM))"),
      TN.mM =  paste("MAGs (Ranked by correlation to total N (mM))"),
      DOC.mM =  paste("MAGs (Ranked by correlation to Dissolved Organic Carbon (mM))")
    )
    yaxis <- y_title_list[[EnvVar]]
    
    p2 <- ggplot(data = mag_CN_plot_df, aes(y = genome_labels, x = subpathway_label)) +
      geom_point(aes(fill = specialisation), shape = 21, size = rel(5)) +
      facet_grid(.~ Versatility + metapath, scales = "free_x", space = "free_x") +
      scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                        guide = guide_legend(override.aes = list(size = rel(5), shape = 21))) +
      ylab(yaxis) + xlab("Functions") +
      theme_bw() +
      theme(text = element_text(size = rel(5)),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y.right =  element_text(angle = 0, hjust = 0),
            strip.background = element_blank(),
            panel.grid = element_blank())
    
    
    # Network stats
    if(DatOrig == "All") {habitat <- c("Palsa", "Bog", "Fen")} else {habitat <- DatOrig}
    mag_deg_cent <- network_stats %>%
      left_join(genome_sizes, by = c("genome_id" = "representative")) %>%
      filter(genome_id %in% mag_CN_plot_df$MAG) %>%
      filter(Habitat %in% habitat) %>%
      rename(MAG = genome_id) %>%
      left_join(input_counts$taxonomy, by = c("MAG" = "genome")) %>%
      mutate(across(Domain:Species, ~gsub("[dpcofgs]__", "", .))) %>%
      mutate(genome_labels = paste0(Phylum, ": ", Genus, " spp.")) %>%
      group_by(MAG) %>%
      mutate(tally = sum(!is.na(degree))) %>%
      # Change names of previously named genomes
      mutate(genome_labels = case_when(MAG == "20170700_S25_26" ~ paste0(Phylum, ": ", "Ca. A. stordalenmirensis"),
                                       MAG == "PLGY01" ~ paste0(Phylum, ": ", "Ca. M. stordalenmirensis"),
                                       MAG == "PMEE01" ~ paste0(Phylum, ": ", "Ca. M. crilli"),
                                       MAG == "20120700_S1X_2" ~ paste0(Phylum, ": ", "Ca. Acididegradans cibantis"),
                                       MAG == "20170700_S25_14" ~ paste0(Phylum, ": ", "Ca. Acidimolecula latum"),
                                       MAG == "3300037175_6" ~ paste0(Phylum, ": ", "Ca. Acidiacetogenum stordalenmirensis"),
                                       MAG == "20170700_S25_9" ~ paste0(Phylum, ": ", "Ca. Acidimolecula amplum"), 
                                       MAG == "20110800_E3D_9" ~ paste0(Phylum, ": ", "Ca. Acididegradans virginiarichiae"),
                                       MAG == "20110800_S1D_14" ~ paste0(Phylum, ": ", "Ca. Abiskobacterium ruthvarnerae"),
                                       MAG == "20160700_S15_2" ~ paste0(Phylum, ": ", "Ca. Acidiflorens viviennsis"),
                                       MAG == "20100900_E1D_10" ~ paste0(Phylum, ": ", "Ca. Acidiflorens anabellensis"),
                                       MAG == "3300036929_5" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum latum"),
                                       MAG == "3300037085_4" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum fermentum"),
                                       MAG %in% c("20110800_E2D_15", "20111000_P1M_1","20111000_P2D_6", "20120500_P15_4",
                                                  "20120500_P23_1", "20120500_P24_7", "20120500_P26_8", "20120500_P36_7",
                                                  "20120700_P2D_10", "20120700_P3D_89", "20130700_C2X_13","20140700_E22_12",
                                                  "20150700_S23_3", "20170700_S25_25", "PKWB01", "PLNE01", "PLYD01", "PLYZ01",
                                                  "PLZA01", "PMOL01", "PNAG01", "PNBE01") ~ paste0(Phylum, ": ", "Ca. Changshengia", " spp."),
                                       TRUE ~ genome_labels)) %>%
      #mutate(genome_labels = paste0(genome_labels, ": ", MAG)) %>% %>%
      #filter(tally > 1) %>% # don't keep MAGs that only appear in network in 1 year
      pivot_longer(matches("centrality|degree|mean_adj_size"), names_to = "NetworkStat", values_to = "value") %>%
      group_by(Habitat, MAG, genome_labels, NetworkStat) %>%
      summarize(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), n = n(),
                se = sd / sqrt(n)) %>%
      ungroup() %>%
      left_join(mag_CN_plot_df %>% select(MAG, MAG_Rank), by = "MAG") %>%
      mutate(MAG = fct_reorder(MAG, MAG_Rank))
    mag_deg_cent.plot <- mag_deg_cent %>%
      select(MAG, mean, NetworkStat, se, sd, n, Habitat) %>%
      distinct() %>%
      ggplot(aes(x = MAG, y = mean, group = NetworkStat)) +
      geom_bar(stat = "identity", aes(fill = Habitat), position = position_dodge2(preserve = "single")) +
      geom_errorbar(data = mag_deg_cent %>% filter(NetworkStat!= "mean_adj_size") %>%
                      select(MAG, mean, NetworkStat, se, sd, n, Habitat) %>%
                      distinct(),
                    aes(ymin = mean - se, ymax = mean + se,
                        group = Habitat), position=position_dodge2(width = 0.1, padding = 0.3, preserve = "single")) +
      scale_fill_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
      #geom_line(aes(color = MAG, group = MAG)) +
      #coord_flip() +
      facet_wrap(~ NetworkStat, scale = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(hjust = 1, angle = 45))
    
    
    return(list(cn_plot = p2, networkplot = mag_deg_cent.plot))
    
  } else {
    mag_CN_plot_df <- network_stats %>%
      filter(genome_id %in% MAGs) %>%
      filter(Habitat == DatOrig) %>%
      rename(MAG = genome_id) %>%
      pivot_longer(matches("centrality|degree"), names_to = "NetworkStat", values_to = "value") %>%
      group_by(MAG, Habitat, NetworkStat) %>%
      summarize(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), n = n(),
                se = sd / sqrt(n)) %>%
      ungroup() %>%
      left_join(input_counts$taxonomy, by = c("MAG" = "genome")) %>%
      mutate(across(Domain:Species, ~gsub("[dpcofgs]__", "", .))) %>%
      mutate(Rank = as.numeric(fct_reorder(MAG, desc(mean)))) %>%
      select(-mean, -sd, -se, -n) %>%
      mutate(genome_labels = paste0(Phylum, ": ", Genus, " spp.", Rank)) %>%
      # Change names of previously named genomes
      mutate(genome_labels = case_when(MAG == "20170700_S25_26" ~ paste0(Phylum, ": ", "Ca. A. stordalenmirensis"),
                                       MAG == "PLGY01" ~ paste0(Phylum, ": ", "Ca. M. stordalenmirensis"),
                                       MAG == "PMEE01" ~ paste0(Phylum, ": ", "Ca. M. crilli"),
                                       MAG == "20120700_S1X_2" ~ paste0(Phylum, ": ", "Ca. Acididegradans cibantis"),
                                       MAG == "20170700_S25_14" ~ paste0(Phylum, ": ", "Ca. Acidimolecula latum"),
                                       MAG == "3300037175_6" ~ paste0(Phylum, ": ", "Ca. Acidiacetogenum stordalenmirensis"),
                                       MAG == "20170700_S25_9" ~ paste0(Phylum, ": ", "Ca. Acidimolecula amplum"), 
                                       MAG == "20110800_E3D_9" ~ paste0(Phylum, ": ", "Ca. Acididegradans virginiarichiae"),
                                       MAG == "20110800_S1D_14" ~ paste0(Phylum, ": ", "Ca. Abiskobacterium ruthvarnerae"),
                                       MAG == "20160700_S15_2" ~ paste0(Phylum, ": ", "Ca. Acidiflorens viviennsis"),
                                       MAG == "20100900_E1D_10" ~ paste0(Phylum, ": ", "Ca. Acidiflorens anabellensis"),
                                       MAG == "3300036929_5" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum latum"),
                                       MAG == "3300037085_4" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum fermentum"),
                                       MAG %in% c("20110800_E2D_15", "20111000_P1M_1","20111000_P2D_6", "20120500_P15_4",
                                                  "20120500_P23_1", "20120500_P24_7", "20120500_P26_8", "20120500_P36_7",
                                                  "20120700_P2D_10", "20120700_P3D_89", "20130700_C2X_13","20140700_E22_12",
                                                  "20150700_S23_3", "20170700_S25_25", "PKWB01", "PLNE01", "PLYD01", "PLYZ01",
                                                  "PLZA01", "PMOL01", "PNAG01", "PNBE01") ~ paste0(Phylum, ": ", "Ca. Changshengia", " spp."),
                                       TRUE ~ genome_labels)) %>%
      #mutate(genome_labels = paste0(genome_labels, ": ", MAG)) %>% %>%
      left_join(pathways_CN$cn_paths_long, by = c("MAG" = "genome")) %>%
      left_join(specialisations, by = c("MAG" = "genome")) %>%
      mutate(specialisations = factor(specialisation, levels = special_levels),
             MAG = fct_reorder(MAG, Rank, .desc = FALSE),
             genome_labels = fct_reorder(genome_labels, Rank, .desc = FALSE)) %>%
      mutate(Versatility = dplyr::case_when(
        subpathway %in% c('reductive_glycine-serine_pathway', 'reductive_glycine-acetyl_pathway', "CAZy-Chitin") ~ 'Both C & N',
        metapath == 'nitrogen' ~ 'Nitrogen',
        TRUE ~ 'Carbon'),
        metapath = ifelse(grepl("Both", Versatility), "both C & N", metapath)
      ) %>%
      mutate(Versatility = factor(Versatility, levels = c("Carbon", "Both C & N", "Nitrogen")),
             metapath = factor(metapath, levels = c("cazymes", "C degradation", "C fermentation", "methanogenesis", "C fixation", "both C & N","nitrogen")))  
    
    xaxis <- "MAGs (Ranked by Increasing Average Degree)"
    
    p2 <- ggplot(data = mag_CN_plot_df, aes(x = genome_labels, y = subpathway_label)) +
      geom_point(aes(fill = specialisation), shape = 21, size = 3) +
      facet_grid(Versatility + metapath ~., scales = "free_y", space = "free_y") +
      scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill,
                        guide = guide_legend(override.aes = list(size = 3, shape = 21))) +
      xlab(xaxis) + ylab("Functions") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y.right =  element_text(angle = 0, hjust = 0),
            strip.background = element_blank(),
            panel.grid = element_blank())
    
    
    
    # Show network info
    # extract order of MAGs
    mag_order <- mag_CN_plot_df %>%
      select(MAG, Rank) %>% distinct()
    mag_deg_cent <- network_stats %>%
      left_join(genome_sizes, by = c("genome_id" = "representative")) %>%
      filter(genome_id %in% MAGs) %>%
      filter(Habitat == DatOrig) %>%
      rename(MAG = genome_id) %>%
      left_join(input_counts$taxonomy, by = c("MAG" = "genome")) %>%
      mutate(across(Domain:Species, ~gsub("[dpcofgs]__", "", .))) %>%
      mutate(genome_labels = paste0(Phylum, ": ", Genus, " spp.")) %>%
      group_by(MAG) %>%
      mutate(tally = sum(!is.na(degree))) %>%
      # Change names of  named genomes
      mutate(genome_labels = case_when(MAG == "20170700_S25_26" ~ paste0(Phylum, ": ", "Ca. A. stordalenmirensis"),
                                       MAG == "PLGY01" ~ paste0(Phylum, ": ", "Ca. M. stordalenmirensis"),
                                       MAG == "PMEE01" ~ paste0(Phylum, ": ", "Ca. M. crilli"),
                                       MAG == "20120700_S1X_2" ~ paste0(Phylum, ": ", "Ca. Acididegradans cibantis"),
                                       MAG == "20170700_S25_14" ~ paste0(Phylum, ": ", "Ca. Acidimolecula latum"),
                                       MAG == "3300037175_6" ~ paste0(Phylum, ": ", "Ca. Acidiacetogenum stordalenmirensis"),
                                       MAG == "20170700_S25_9" ~ paste0(Phylum, ": ", "Ca. Acidimolecula amplum"), 
                                       MAG == "20110800_E3D_9" ~ paste0(Phylum, ": ", "Ca. Acididegradans virginiarichiae"),
                                       MAG == "20110800_S1D_14" ~ paste0(Phylum, ": ", "Ca. Abiskobacterium ruthvarnerae"),
                                       MAG == "20160700_S15_2" ~ paste0(Phylum, ": ", "Ca. Acidiflorens viviennsis"),
                                       MAG == "20100900_E1D_10" ~ paste0(Phylum, ": ", "Ca. Acidiflorens anabellensis"),
                                       MAG == "3300036929_5" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum latum"),
                                       MAG == "3300037085_4" ~ paste0(Phylum, ": ", "Ca. Acidiaspiraculum fermentum"),
                                       MAG %in% c("20110800_E2D_15", "20111000_P1M_1","20111000_P2D_6", "20120500_P15_4",
                                                  "20120500_P23_1", "20120500_P24_7", "20120500_P26_8", "20120500_P36_7",
                                                  "20120700_P2D_10", "20120700_P3D_89", "20130700_C2X_13","20140700_E22_12",
                                                  "20150700_S23_3", "20170700_S25_25", "PKWB01", "PLNE01", "PLYD01", "PLYZ01",
                                                  "PLZA01", "PMOL01", "PNAG01", "PNBE01") ~ paste0(Phylum, ": ", "Ca. Changshengia", " spp."),
                                       TRUE ~ genome_labels)) %>%
      #mutate(genome_labels = paste0(genome_labels, ": ", MAG)) %>%
      #filter(tally > 1) %>% # don't keep MAGs that only appear in network in 1 year
      pivot_longer(matches("centrality|degree|mean_adj_size"), names_to = "NetworkStat", values_to = "value") %>%
      group_by(MAG, Habitat, genome_labels, NetworkStat) %>%
      summarize(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), n = n(),
                se = sd / sqrt(n)) %>%
      ungroup() %>%
      left_join(mag_order, by = "MAG") %>%
      mutate(MAG = fct_reorder(MAG, Rank))
    mag_deg_cent.plot <- mag_deg_cent %>%
      ggplot(aes(x = MAG, y = mean, group = NetworkStat)) +
      geom_bar(stat = "identity") +
      geom_errorbar(data = mag_deg_cent %>% filter(NetworkStat!= "mean_adj_size"),
                    aes(ymin = mean - se, ymax = mean + se), width = 0.5) +
      #geom_point(aes(color = MAG)) +
      #geom_line(aes(color = MAG, group = MAG)) +
      #coord_flip() +
      facet_wrap(~ NetworkStat, scale = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(hjust = 1, angle = 45))
    
    
    
    return(list(cn_plot = p2, networkplot = mag_deg_cent.plot))
    
  }

}


#test <- plot_CN_met(EnvVar = "d13C_CH4", DatOrig = "All")
# Save CN plots
# Save all CN plot combos:
save_all_CN <- function(combos){
  DatOrig <- combos[[1]]
  if(DatOrig == "All") {hab <- "all habitats"} else {hab <- paste(str_to_lower(DatOrig), "samples")}
  Mod_col <- combos[[2]]
  EnvVar <- combos[[3]]
  
  CN_combo_plot <- plot_CN_met(EnvVar = EnvVar, DatOrig = DatOrig, 
                               MAGs = NULL, Significant_only = T)
  ggsave(plot = CN_combo_plot$cn_plot, str_c(DatOrig, Mod_col, EnvVar, "_VIP_CN_plot.png", 
                                            sep = "_"), 
         height = 11, width = 14,
         path = paste0(output_dir, "/longNames"), dpi = 300)
  ggsave(plot = CN_combo_plot$networkplot, str_c(DatOrig, Mod_col, EnvVar, "_VIP_network_plot.png", 
                                             sep = "_"), 
         height = 10, width = 10,
         path = paste0(output_dir, "/longNames"), dpi = 300)
}



wgcna_combos %>%
  split(1:nrow(.)) %>%
#  slice_head() %>%
  purrr::map(., ~save_all_CN(.))


CN_combo_plot <- plot_CN_met(EnvVar = "alphaC", DatOrig = "Bog", 
                             MAGs = NULL, Significant_only = T)
ggsave(plot = CN_combo_plot$cn_plot + theme(axis.text.y = element_text(face = "bold.italic")), 
       str_c("Bog", "blue", "alphaC", "_VIP_CN_plot.png", 
                                           sep = "_"), 
       height = 13, width = 18,
       path = paste0(output_dir, "/longNames"), dpi = 300)

################################################################################
# Combining CN metabolism and VIP plot for Fen and bog
################################################################################
# Fen Yellow C% plot significant
fen.cpercent_hive <- plot_hive(EnvVar = "C.percent", DatOrig = "Fen", check_cor = TRUE,
          y_title = "Significant correlations\n to %C across fen samples", Significant_only = T)

sig_mags_fen <- fen.cpercent_hive$data %>%
  filter(p.adj < 0.05)
fen.cpercent_cn <- plot_CN_met(EnvVar = "C.percent", DatOrig = "Fen",
                               Significant_only = T, sigMags = unique(sig_mags_fen$genome))


p1 <- fen.cpercent_hive$plot
p2 <- fen.cpercent_cn$cn_plot
fen.cpercent_grid <- plot_grid(NULL, p1, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2,7))

ggsave(plot = fen.cpercent_grid, "fen.cpercent_derek_tears.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)


# Bog Blue alphaC plot
bog.alphac_hive <- plot_hive(EnvVar = "alphaC", DatOrig = "Bog", check_cor = TRUE,
                               y_title = expression(paste("Correlation to ", alpha[c], " across bog samples")),
                               Significant_only = T, mag_text = F)
sig_mags_bog <- bog.alphac_hive$data %>%
  filter(p.value < 0.05)


bog.alphaC_cn <- plot_CN_met(EnvVar = "alphaC", DatOrig = "Bog", 
                             Significant_only = T, sigMags = unique(sig_mags_bog$genome))


p1 <- bog.alphac_hive$plot
p2 <- bog.alphaC_cn$cn_plot
bog.alphaC_grid <- plot_grid(NULL, p1, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2.1,7))

ggsave(plot = bog.alphaC_grid, "bog.alphaC_derek_tears.png", 
       height = 16, width = 15,
       path = output_dir, dpi = 300)


# All including non-significant
# Fen Yellow C% plot
fen.cpercent_hive <- plot_hive(EnvVar = "C.percent", DatOrig = "Fen", check_cor = TRUE,
                               y_title = "Correlations\n to %C across fen samples", Significant_only = F)
fen.cpercent_cn <- plot_CN_met(EnvVar = "C.percent", DatOrig = "Fen")


p1 <- fen.cpercent_hive$plot
p2 <- fen.cpercent_cn$cn_plot
fen.cpercent_grid_ns <- plot_grid(NULL, p1, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2,7))

ggsave(plot = fen.cpercent_grid_ns, "fen.cpercent_derek_tears_nonsig.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)


# Bog Blue alphaC plot
bog.alphac_hive <- plot_hive(EnvVar = "alphaC", DatOrig = "Bog", check_cor = TRUE,
                             y_title = expression(paste("Correlations\n to ", alpha[c], "  across bog samples")))

bog.alphaC_cn <- plot_CN_met(EnvVar = "alphaC", DatOrig = "Bog")


p1 <- bog.alphac_hive$plot
p2 <- bog.alphaC_cn$cn_plot
bog.alphaC_grid_ns <- plot_grid(NULL, p1, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2.1,7))

ggsave(plot = bog.alphaC_grid_ns, "bog.alphaC_derek_tears_nonsig.png", 
       height = 16, width = 15,
       path = output_dir, dpi = 300)

################################################################################
# Tears for all modules
################################################################################
so_many_tears_modules <- function(combos) {
  
  DatOrig <- combos[[1]]
  if(DatOrig == "All") {hab <- "all habitats"} else {hab <- paste(str_to_lower(DatOrig), "samples")}
  Mod_col <- combos[[2]]
  EnvVar <- combos[[3]]
  significant_only <-combos[[4]]
  # setup appropriate y-title
  y_title_list <- list(
    d15N_peat = paste("Correlation to \u03B4\U00B9\u2075", "N across", hab),
    d13C_CO2 = paste("Correlation to \u03B4\U00B9\u00B3", "C CO\u2082 across", hab),
    d13C_CH4 =  paste("Correlation to \u03B4\U00B9\u00B3","C CH\u2084", " across ", hab),
    d13C_peat =  paste("Correlation to \u03B4\U00B9\u00B3","C peat", " across ", hab),
    alphaC =  paste("Correlation to \u0251C", " across ", hab),
    N.percent =  paste("Correlation to total %N across ", hab),
    C.percent =  paste("Correlation to total %C across ", hab),
    CO2.mM =  paste("Correlation to [CO\u2082](mM) across ", hab),
    TN.mM =  paste("Correlation to total N (mM) across ", hab),
    DOC.mM =  paste("Correlation to Dissolved Organic Carbon (mM) across ", hab)
  )
  y_title <- y_title_list[[EnvVar]]
  
  
  
  module_hive <- plot_hive(EnvVar = EnvVar, DatOrig = DatOrig, check_cor = TRUE,
                           y_title = y_title, Significant_only = significant_only)
  
  # setup plotting of only significant mags
  if(significant_only) {
    sig_title <- "sigonly"
    sig_mags_data <- module_hive$data %>%
      filter(p.adj < 0.05)
    
    # exit if no significant mags
    if(nrow(sig_mags_data)<1) {
      print(paste0("No significant correlations in ", DatOrig, " ", Mod_col, " ", EnvVar))
      return(combos)
    }
    
    module_cn <- plot_CN_met(EnvVar = EnvVar, DatOrig = DatOrig, 
                                 Significant_only = T, sigMags = unique(sig_mags_data$genome))
  } else {
    sig_title <- ""
    module_cn <- plot_CN_met(EnvVar = EnvVar, DatOrig = DatOrig, Significant_only = F)
  }
  
  # arrange plots
  p1 <- module_hive$plot
  p2 <- module_cn$cn_plot
  module_grid <- plot_grid(NULL, p1, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2.1,7))
  
  ggsave(plot = module_grid, paste0(DatOrig,"_", EnvVar,"_",Mod_col, "_derek_tears_",sig_title,".png"), 
         height = 16, width = 15,
         path = output_dir, dpi = 300)
}

# All tears significant only
wgcna_combos %>%
  mutate(Significant_only = T) %>%
  split(1:nrow(.)) %>%
  #  slice_head() %>%
  purrr::map(., ~so_many_tears_modules(.))

# All tears, including non-significant
wgcna_combos %>%
  mutate(Significant_only = F) %>%
  split(1:nrow(.)) %>%
  #  slice_head() %>%
  purrr::map(., ~so_many_tears_modules(.))

################################################################################
# CN and network plots for methanogens
################################################################################
# Methanogens only

fen_methanogens <- c("20110700_S2D_21", "20110800_E1S_5", "20110800_E3D_10", "20110800_E3D_36", "20120600_E2D_39", "3300037104_19", 
                     "PMGM01", "20120800_E2X_5", "3300037365_14", "PMKW01", "20140700_E14_1", "PLGY01", "PMMV01", "3300037062_11", 
                     "PMEE01", "PMNG01")

bog_methanogens <- c("20110700_S2D_21", "20120700_S1D_59", "PLGY01", "PLTK01", "PMEE01", "PMGM01")

p2 <- plot_CN_met(DatOrig = "Fen", MAGs = fen_methanogens)

ggsave(plot = p2$cn_plot, "fen_methanogens_cn_plot.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)

ggsave(plot = p2$networkplot, "fen_methanogens_network_connection.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)

p3 <- plot_CN_met(DatOrig = "Bog", MAGs = bog_methanogens)

ggsave(plot = p3$cn_plot, "bog_methanogens_cn_plot.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)

ggsave(plot = p3$networkplot, "Bog_methanogens_network_connection.png", 
       height = 10, width = 10,
       path = output_dir, dpi = 300)



q <- plot_grid(NULL, p, p2, ncol = 1, align = "v", axis = "lrb", rel_heights = c(0.5,2,7))
q
ggsave(plot = q, str_c("comb", environmental_var, habitat, "VIP_specials.png", 
                       sep = "_"), 
       height = 10, width = 10,
       path = output_dir, dpi = 300)

################################################################################
# Genome Size and Functions Plots
################################################################################

methanogen_origin_list <- bind_rows(
  data.frame(methanogen_origin = "fen_methanogens", genome = fen_methanogens),
  data.frame(methanogen_origin = "bog_methanogens", genome = bog_methanogens))

CN_genome_size.plot <- genome_sizes %>%
  rename(genome = representative) %>%
  left_join(pathways_CN$cn_paths_long, by = c("genome")) %>%
  mutate(Versatility = ifelse(subpathway_label %in% c("chitin", "rGlycine serine", "rGlycine acetyl"), "both C & N", Versatility)) %>%
  left_join(specialisations, by = "genome") %>%
  select(genome, mean_adj_size, specialisation, Versatility) %>%
  distinct() %>%
  ggplot(aes(y = mean_adj_size, x = specialisation)) +
  geom_point(aes(color = Versatility), shape = 16, position = position_jitterdodge(jitter.width = 0.3), alpha = 0.1) +
  geom_boxplot(aes(color = Versatility), fill = NA, position = position_dodge2()) +
  ylab("genome size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = genome_size_CN_capability, str_c("methanogen_cn_capability_genome_size.png", 
                                               sep = "_"), 
       height = 10, width = 12,
       path = output_dir, dpi = 300)



metapath_genome_size_boxplot <- genome_sizes %>%
  rename(genome = representative) %>%
  left_join(pathways_CN$cn_paths_long, by = c("genome")) %>%
  mutate(Versatility = ifelse(subpathway_label %in% c("chitin", "rGlycine serine", "rGlycine acetyl"), "both C & N", Versatility)) %>%
  group_by(genome, mean_adj_size, Versatility, metapath) %>%
  tally() %>%
  rename(n_metapath = n) %>%
  left_join(specialisations, by = "genome") %>%
  ungroup() %>%
  select(genome, mean_adj_size, specialisation, metapath, n_metapath) %>%
  mutate(metapath = factor(metapath, levels = metapathway_levels)) %>%
  distinct() %>%
  ggplot(aes(y = mean_adj_size, x = specialisation)) +
  geom_point(aes(color = metapath), shape = 16, position = position_jitterdodge(jitter.width = 0.3), alpha = 0.1) +
  geom_boxplot(aes(fill = metapath), alpha = 0.3, position = position_dodge2(width = 0.3)) +
  guides(color = "none") +
  scale_fill_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  scale_color_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  ylab("genome size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
metapath_genome_size_boxplot

ggsave(plot = metapath_genome_size_boxplot, str_c("methanogen_genome_size_allgenomes_boxplot.png", 
                                               sep = "_"), 
       height = 10, width = 12,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)


# Test for differences between genome sizes

metapath_genome_size.df <- genome_sizes %>%
  rename(genome = representative) %>%
  left_join(pathways_CN$cn_paths_long, by = c("genome")) %>%
  mutate(Versatility = ifelse(subpathway_label %in% c("chitin", "rGlycine serine", "rGlycine acetyl"), "both C & N", Versatility)) %>%
  group_by(genome, mean_adj_size, metapath) %>%
  count(.drop = F) %>% ungroup() %>%
  complete(nesting(genome, mean_adj_size), metapath) %>% # get a count of metapaths that are 0 counts
  filter(!is.na(metapath)) %>%
#  mutate(n = ifelse(is.na(n), 0, n)) %>%
  rename(n_metapath = n) %>%
  left_join(specialisations, by = "genome") %>%
  ungroup() %>%
  select(genome, mean_adj_size, specialisation, metapath, n_metapath) %>%
  mutate(metapath = factor(metapath, levels = metapathway_levels)) %>%
  distinct() %>%
  mutate(Methanogen_origin = ifelse(genome %in% intersect(bog_methanogens, fen_methanogens), 
                                    "Fen and Bog", 
                             ifelse(genome %in% bog_methanogens, "Bog", 
                             ifelse(genome %in% fen_methanogens, "Fen", NA))))
metapath_genome_size.plot <- metapath_genome_size.df %>%
  ggplot(aes(y = mean_adj_size, x = n_metapath)) +
  geom_point(aes(color = metapath, fill = metapath), shape = 16, alpha = 0.1) +
  geom_point(data = metapath_genome_size.df %>% filter(!is.na(Methanogen_origin)), 
             aes(fill = metapath, shape = Methanogen_origin), alpha = 0.3) +
  geom_smooth(aes(color = metapath), method = "lm", se = F) +
  facet_wrap(~specialisation, scale = "free") + 
  scale_color_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  scale_fill_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  scale_shape_manual(name = "Methanogen Origin", breaks = c("Fen and Bog", "Fen", "Bog"), values = c(21, 22,24)) +
  ylab("genome size") + xlab("number of metapathways") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
metapath_genome_size.plot

ggsave(plot = metapath_genome_size.plot, str_c("metapath_genome_size_allgenomes.png", 
                                               sep = "_"), 
       height = 10, width = 12,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)


### Testing for differences in size by specialisation only (no metapaths)
genome_size_specialization.df <- genome_sizes %>%
  rename(genome = representative) %>%
  left_join(specialisations, by = "genome") %>%
  filter(!is.na(specialisation)) %>%
  select(genome, mean_adj_size, specialisation) %>%
  mutate(Methanogen_origin = ifelse(genome %in% intersect(bog_methanogens, fen_methanogens), 
                                    "Fen and Bog", 
                                    ifelse(genome %in% bog_methanogens, "Bog", 
                                           ifelse(genome %in% fen_methanogens, "Fen", NA))))


genome_size_test.df <- genome_size_specialization.df %>%
  select(specialisation, mean_adj_size, genome) %>%
  distinct() %>%
  group_by(specialisation) %>%
  mutate(mean_genome_grp = mean(mean_adj_size), sd_grp = sd(mean_adj_size),
         se = sd_grp/ sqrt(n()))

write_tsv(genome_size_test.df %>%
            select(specialisation, mean_genome_grp, sd_grp, se) %>%
            distinct() %>%
            mutate(across(where(is.numeric), ~signif(.x, digits = 3))),
          file = paste0(output_dir, "/genome_size_plots", "/mean_std_dev_se_genome_size_specialisation.txt"))


# anova:
genome_special.aov <- aov(mean_adj_size ~ specialisation, data = genome_size_test.df)
summary(genome_special.aov)
#plot(genome_special.aov)

genome_special_tukey <- TukeyHSD(genome_special.aov)

write_tsv(genome_special_tukey$specialisation %>%
            data.frame()%>%
            rownames_to_column(var = "comparison") %>%
            rename(genome_size_difference = `diff`,
                   lower_95_conf_int = `lwr`,
                   upper_95_conf_int = `upr`) %>%
            mutate(across(where(is.numeric), ~signif(.x, digits = 3))),
          file = paste0(output_dir, "/genome_size_plots", "/tukey_test_results_genome_size_specialisation.txt"))

# Tukey labels
generate_label_df <- function(TUKEY, variable){
  library(multcompView)
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

labels <- generate_label_df(genome_special_tukey, "specialisation")#generate labels using function

names(labels)<-c('Letters','specialisation')#rename columns for merging

tukey_letters <- genome_size_test.df %>%
  group_by(specialisation, mean_genome_grp) %>%
  summarize(max = max(mean_adj_size)) %>%
  mutate(mean_adj_size = max + max*0.2) %>% # make label position 20% higher than maxiumum value
  left_join(labels, by = c("specialisation"))

specialisation_genome_size.plot <- genome_size_test.df %>%
  ggplot(aes(y = mean_adj_size, x = fct_reorder(specialisation, mean_genome_grp, .desc = TRUE))) +
  geom_boxplot() +
  geom_point(aes(fill = specialisation), shape = 21, color = "black", alpha = 0.1,
             position = "jitter") +
  geom_text(data = tukey_letters, aes(label = Letters)) +
  scale_color_manual(name = "Specialisation", breaks = special_levels, values = special_colour) +
  scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_colour) +
  ylab("Genome size") + xlab("Specialisations") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
specialisation_genome_size.plot

ggsave(plot = specialisation_genome_size.plot, str_c("specialisation_genome_size_allgenomes.png", 
                                               sep = "_"), 
       height = 6, width = 12,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)








genome_size_CN_capability <- methanogen_origin_list %>%
  left_join(pathways_CN$cn_paths_long, by = "genome") %>%
  mutate(Versatility = ifelse(subpathway_label %in% c("chitin", "rGlycine serine", "rGlycine acetyl"), "both C & N", Versatility)) %>%
  group_by(methanogen_origin, genome, Versatility) %>% add_tally() %>%
  rename(`NoOfCNPathways` = n) %>% ungroup() %>%
  distinct() %>%
  complete(nesting(genome, methanogen_origin), Versatility) %>% # get a count of metapaths that are 0 counts
  filter(!is.na(Versatility)) %>%
  mutate(NoOfCNPathways = ifelse(is.na(NoOfCNPathways), 0, NoOfCNPathways)) %>%
  ungroup() %>%
  left_join(genome_sizes, by = c("genome" = "representative")) %>%
  select(NoOfCNPathways, mean_adj_size, genome, Versatility, methanogen_origin) %>%
  distinct() %>%
  ggplot(aes(x = NoOfCNPathways, y = mean_adj_size)) +
  geom_point() +
  geom_smooth(alpha = 0.4, method = "lm") +
  geom_text_repel(aes(label = genome)) +
  facet_grid(methanogen_origin~Versatility) +
  ylab("genome size") + xlab("Number of pathways in carbon, nitrogen, or joint categories") +
  theme_bw()

ggsave(plot = genome_size_CN_capability, str_c("methanogen_cn_capability_genome_size.png", 
                       sep = "_"), 
       height = 10, width = 12,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)


genome_size_CN_capability_metapath <- methanogen_origin_list %>%
  left_join(pathways_CN$cn_paths_long, by = "genome") %>%
  mutate(metapath = ifelse(subpathway_label %in% c("chitin", "rGlycine serine", "rGlycine acetyl"), "both C & N", metapath)) %>%
  select(-Versatility) %>% distinct() %>%
  group_by(methanogen_origin, genome, metapath) %>% add_tally() %>%
  rename(`NoOfMetaPathways` = n) %>%
  ungroup() %>%
  complete(nesting(genome, methanogen_origin), metapath) %>% # get a count of metapaths that are 0 counts
  filter(!is.na(metapath)) %>%
  mutate(NoOfMetaPathways = ifelse(is.na(NoOfMetaPathways), 0, NoOfMetaPathways)) %>%
  ungroup() %>%
  select(genome, methanogen_origin, metapath, NoOfMetaPathways) %>%
  distinct() %>%
  left_join(specialisations, by = c("genome")) %>%
  left_join(genome_sizes, by = c("genome" = "representative")) %>%
  select(NoOfMetaPathways, mean_adj_size, genome, metapath, methanogen_origin) %>%
  ggplot(aes(x = NoOfMetaPathways, y = mean_adj_size)) +
  geom_point() +
  geom_smooth(alpha = 0.4, method = "lm") +
  geom_text_repel(aes(label = genome), size = rel(3), alpha = 0.5) +
  facet_grid(methanogen_origin~metapath) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  ylab("genome size") + xlab("Number of pathways in each metapathway") +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank()) 
ggsave(plot = genome_size_CN_capability_metapath, str_c("methanogen_cn_capability_genome_size_metapath.png", 
                                               sep = "_"), 
       height = 6, width = 20,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)

################################################################################
# Genome Size and network metrics
################################################################################
network_stats_summary <- network_stats %>%
  group_by(genome_id, Habitat) %>%
  summarize(across(matches("centrality|degree"), list(mean = ~mean(., na.rm = T), sd = ~sd(.x, na.rm = T),
                                                      se = ~sd(.x, na.rm = T)/ sqrt(n()))))

metapath_genome_size.df <- genome_sizes %>%
  rename(genome = representative) %>%
  left_join(pathways_CN$cn_paths_long, by = c("genome")) %>%
  group_by(genome, mean_adj_size, metapath) %>%
  count(.drop = F) %>% ungroup() %>%
  complete(nesting(genome, mean_adj_size), metapath) %>% # get a count of metapaths that are 0 counts
  filter(!is.na(metapath)) %>%
  mutate(n = ifelse(is.na(n), 0, n), # replace NAs with 0 counts for metapathways
         n = ifelse(n>0, 1, n)) %>% # make any count >0 = 1 for presence/absence of metapathways
  rename(metapath_pa = n) %>%
  distinct() %>%
  left_join(specialisations, by = "genome") %>%
  ungroup() %>%
  select(genome, mean_adj_size, specialisation, metapath, metapath_pa) %>%
  mutate(metapath = factor(metapath, levels = metapathway_levels)) %>%
  distinct() %>%
  left_join(network_stats_summary, by = c("genome" = "genome_id")) %>%
  # filter organisms not found in any habitat's network
  filter(!is.na(Habitat)) %>%
  filter(!is.na(specialisation))

metapath_genome_size.plot <- metapath_genome_size.df %>%
  ggplot(aes(y = mean_adj_size, x = degree_mean)) +
  geom_point(aes(color = metapath, fill = metapath, shape = as.factor(metapath_pa)), alpha = 0.2) +
  geom_smooth(aes(color = metapath, linetype = as.factor(metapath_pa)), method = "lm", se = F) +
  facet_wrap(Habitat~specialisation, scale = "free", drop = T) + 
  scale_color_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  scale_fill_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  ylab("genome size") + xlab("Degree (average across years)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
metapath_genome_size.plot


# No metapath colors
degree_genome_size.plot <- metapath_genome_size.df %>%
  select(-metapath, -metapath_pa) %>% distinct() %>%
  ggplot(aes(y = mean_adj_size, x = degree_mean)) +
  geom_point( alpha = 0.2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(Habitat~specialisation, scale = "free", drop = T) + 
  #scale_color_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  #scale_fill_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  ylab("genome size") + xlab("Degree (average across years)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
degree_genome_size.plot

ggsave(plot = degree_genome_size.plot, str_c("degree_genome_size.png", 
                                                        sep = "_"), 
       height = 6, width = 20,
       path = paste0(output_dir, "/genome_size_plots"), dpi = 300)


metapath_genome_size.df %>%
  select(-metapath, -metapath_pa) %>% distinct() %>%
  ggplot(aes(y = mean_adj_size, x = betweeness_centrality_mean)) +
  geom_point( alpha = 0.2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(Habitat~specialisation, scale = "free", drop = T) + 
  #scale_color_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  #scale_fill_manual(name = "Metapathways", breaks = metapathway_levels, values = colour_metapathway) +
  ylab("genome size") + xlab("Betweenness Centrality (average across years)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Running linear models with Metapathways, degree, and centrality

# Check that degree and centrality are not too correlated
ggplot(data = metapath_genome_size.df, 
       aes(x = degree_mean, y = betweeness_centrality_mean)) +
  geom_point(aes(color = Habitat)) + 
  geom_smooth(method = "lm")
cor.test(data = metapath_genome_size.df, ~ degree_mean + betweeness_centrality_mean)

# refromat dataframe to pivot wide metapath presence/absence
genome_size_metric_data <- metapath_genome_size.df %>%
  pivot_wider(names_from = metapath, values_from = metapath_pa) %>%
  mutate(log_mean_adj_size = log(mean_adj_size))
genome_size_metric_data$specialisation
names(genome_size_metric_data)
genome_size_metric_data %>%
  ggplot(aes(x = log(mean_adj_size))) +
  geom_histogram() +
  facet_wrap(~Habitat)

# seems a bit skewed, probably could take the log, or do a poisson distribution (we'll do log)

library(lme4)

methanogen.genome_size_metric_data <- metapath_genome_size.df %>% 
  filter(metapath_pa == 1) %>%
  mutate(log_mean_adj_size = log(mean_adj_size),
         metapath = factor(metapath, levels = metapathway_levels)) %>%
  filter(specialisation == "methanogen") %>%
  filter(Habitat == "Bog") %>%
  mutate(across(all_of(c("degree_mean", "betweeness_centrality_mean")), ~as.vector(scale(.)), .names = "{.col}_scale"))

genome_size_mod_all <- lmer(log_mean_adj_size ~ degree_mean_scale + betweeness_centrality_mean_scale  + (1 | metapath),
                            data = methanogen.genome_size_metric_data)
summary(genome_size_mod_all)

tdat <- data.frame(predicted=predict(genome_size_mod_all), 
                   residual = residuals(genome_size_mod_all), 
                   methanogen.genome_size_metric_data)
# Resid/predict; looks great
tdat %>% 
ggplot(aes(x = predicted, y = residual, color = metapath)) + geom_point() + geom_hline(yintercept = 0, lty = 3) +
  scale_color_manual(values = colour_metapathway, breaks = metapathway_levels)

ggplot(tdat,aes(x=residual)) + geom_histogram(bins=20, color="black")

# qq plot, not great at the ends... could be better
ggplot(tdat,aes(sample=residual)) + stat_qq() + stat_qq_line()


# prediction plot
ggplot(tdat,aes(x=degree_mean_scale,y=predicted,colour=metapath, group=metapath)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = colour_metapathway, breaks = metapathway_levels) +
  ylab("predicted genome size") +
  theme(legend.position="bottom", legend.direction = "horizontal")




# Generalists
generalist.genome_size_metric_data <- metapath_genome_size.df %>% 
  filter(metapath_pa == 1) %>%
  mutate(log_mean_adj_size = log(mean_adj_size),
         metapath = factor(metapath, levels = metapathway_levels)) %>%
  filter(specialisation == "generalist") %>%
  mutate(across(all_of(c("degree_mean", "betweeness_centrality_mean")), ~as.vector(scale(.)), .names = "{.col}_scale"))

genome_size_mod_all <- lmer(log_mean_adj_size ~ degree_mean_scale + betweeness_centrality_mean_scale + 
                              (1| Habitat) + (1 | metapath),
                            data = generalist.genome_size_metric_data)
summary(genome_size_mod_all)

tdat <- data.frame(predicted=predict(genome_size_mod_all), 
                   residual = residuals(genome_size_mod_all), 
                   generalist.genome_size_metric_data)
# Resid/predict; looks great
tdat %>% 
  ggplot(aes(x = predicted, y = residual, color = metapath, shape = Habitat)) + geom_point() + geom_hline(yintercept = 0, lty = 3) +
  scale_color_manual(values = colour_metapathway, breaks = metapathway_levels)

ggplot(tdat,aes(x=residual)) + geom_histogram(bins=20, color="black")

# qq plot, not great at the ends... could be better
ggplot(tdat,aes(sample=residual)) + stat_qq() + stat_qq_line()


# prediction plot
ggplot(tdat,aes(x=degree_mean_scale,y=predicted,colour=metapath, group=metapath)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = colour_metapathway, breaks = metapathway_levels) +
  ylab("predicted genome size") +
  theme(legend.position="bottom", legend.direction = "horizontal")




genome_size_mod_simple <- lmer(log_mean_adj_size ~  `C fixation_scale` + `C fermentation_scale` +
                                 `cazymes_scale` + `nitrogen_scale` + degree_mean_scale + betweeness_centrality_mean_scale + 
                                 (1| Habitat),
                               data = methanogen.genome_size_metric_data)
summary(genome_size_mod_all)


genome_size_mod_simple <- lmer(log_mean_adj_size ~  `C fixation_scale` + `C fermentation_scale` +
                              `cazymes_scale` + `nitrogen_scale` + degree_mean_scale + betweeness_centrality_mean_scale + 
                              (1| Habitat),
                            data = methanogen.genome_size_metric_data)
summary(genome_size_mod_all)


#################################################################################
# Functional redundancy - How many macromolecule degradiers are encoded by multiple generalisits 
#################################################################################
genome_redundnancy_in_c_processing <- specialisations %>%
  filter(specialisation == "generalist") %>%
  select(genome) %>% distinct() %>%
  mutate(total_generalists = n()) %>%
  left_join(product_refined %>% filter(call), by = "genome") %>%
  right_join(metapathway_groups, by = "subpathway") %>%
  select(metapath, subpathway_label, total_generalists, genome) %>%
  distinct() %>%
  group_by(metapath, subpathway_label, total_generalists) %>%
  tally() %>%
  mutate(proportion_generalists_containing_pathway = n/total_generalists) %>%
  filter(metapath %in% c("cazymes", "C degradation", "C fermentation"))


genome_redundnancy_in_c_processing %>% 
  filter(proportion_generalists_containing_pathway > 0.5) # 


genome_redundnancy_in_c_processing %>%
  ungroup() %>%
  select(metapath, n) %>%
  rename(count = n) %>%
  group_by(metapath) %>%  
  summarize(mean_num_generalists = mean(count),
            median_num_generalists = median(count),
            total_subpathways_in_metapath = n(),
            sd = sd(count))



# Collector's curve 
func_red_matrix <- specialisations %>%
  filter(specialisation == "generalist") %>%
  select(genome) %>% distinct() %>%
  mutate(total_generalists = n()) %>%
  left_join(product_refined %>% filter(call), by = "genome") %>%
  right_join(metapathway_groups, by = "subpathway") %>%
  select(metapath, subpathway_label, genome) %>%
  distinct() %>%
  filter(metapath %in% c("cazymes", "C degradation")) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "subpathway_label", id_cols = "genome", values_from = "value", values_fill = 0) %>%
  filter(!is.na(genome)) %>%
  column_to_rownames(var = "genome") 


func_red_matrix <- specialisations %>%
  select(specialisation, genome) %>% distinct() %>%
  group_by(specialisation) %>%
  mutate(total_in_specialisation = n()) %>%
  left_join(product_refined %>% filter(call), by = "genome") %>%
  right_join(metapathway_groups, by = "subpathway") %>%
  select(specialisation, metapath, subpathway_label, genome) %>%
  distinct() %>%
  mutate(value = 1) %>%
  group_by(specialisation) %>% nest() %>%
  mutate(data = purrr::map(data, ~ .x %>%
      pivot_wider(names_from = "subpathway_label", id_cols = "genome", values_from = "value", values_fill = 0) %>% 
        filter(!is.na(genome)) %>%
        column_to_rownames(var = "genome")
      )) %>% 
  mutate(specaccum = purrr::map(data, ~ specaccum(.x)),
         number_of_subpathways = purrr::map(specaccum, ~ .x[[4]]),
         organisms = purrr::map(specaccum, ~ .x[[3]]),
         sd = purrr::map(specaccum, ~ .x[[5]])) %>%
  unnest(c(number_of_subpathways, organisms, sd))
  
ggplot(func_red_matrix,
       aes(x = organisms, y = number_of_subpathways)) +
  geom_point(aes(color = specialisation)) +
  scale_color_manual(values = special_colour, breaks = special_levels) +
  theme_bw()


rare_curve.plot <- vegan::rarecurve(func_red_matrix, label = F, se = T, tidy = T)

spec_accu <- specaccum(func_red_matrix)
ggplot(data.frame(No_subpathway = spec_accu$richness, organisms = spec_accu$sites, spec_accu$sd),
       aes(x = organisms, y = No_subpathway)) +
  geom_point() +
  theme_bw()

