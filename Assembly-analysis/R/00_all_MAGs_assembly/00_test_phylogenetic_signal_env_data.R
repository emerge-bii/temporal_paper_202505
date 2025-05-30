#' ## Test phylogenetic signal
#' This tests the phylogenetic signal in our data to see if assembly analysis is
#' appropriate; Note: this script contains parts that likely will need to be run
#' on super computer. You may kill your R session if you attempt to run it locally.

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(ape)
library(vegan)
library(analogue)
library(parallel)
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

# Remove any outputs that were already generated
lapply( list.files(here(figures.fp, "env_var_histograms"), pattern = ".png", 
                   full.names = T), function(x) {
  file.remove(x)
})
lapply( list.files(here(figures.fp, "optima_histograms"), pattern = ".png", 
                   full.names = T), function(x) {
                     file.remove(x)
                   })
lapply( list.files(here(figures.fp, "optima_plots"), pattern = "_opt_filttree.png", 
                   full.names = T), function(x) {
                     file.remove(x)
                   })

lapply( list.files(here(figures.fp, "optima_plots"), pattern = "_opt_tree.png", 
                   full.names = T), function(x) {
                     file.remove(x)
                   })

lapply( list.files(outputs.fp, pattern = "mcorr_", full.names = T), function(x) {
  file.remove(x)
})

file.remove(here(outputs.fp, "/mcorr_list.RDS"))
file.remove(paste0(figures.fp, "/phylogentic_dist_mantel_correlogram.png"))
file.remove(paste0(figures.fp, "/phylogentic_dist_mantel_correlogram_size.png"))

#### ====================================================================== ####
#### Define some helper functions
#### ====================================================================== ####
# Calculate niche optima

calc_niche_optima <- function(otu_table = input_counts$otu_table[,-1],
                              env_data_table = input_counts$sample_metadata,
                              column_to_calc_opt = "DepthAvg__",
                              sample_id_column = "temporal_sample_id",
                              OTU_ubiq_cutoff = 0.15,
                              tol = FALSE, verbose = FALSE) {
  # arguments
  # otu_table ---------- the species by sample table (with species as rownames, and 
  #                      samples as columns)
  # env_data_table ----- the environmental data table with a sample id column matching
  #                      the sample id column names that are in the otu table
  # column_to_calc_opt - the environmental data that will be used to caluclate niche
  #                      optima. The column should be continuous and numeric. 
  #                      Samples that are NA will be filtered out prior to optima 
  #                      calculation; should be specified as a string.
  # sample_id_column --- the sample id column name. 
  # OTU_ubiq_cutoff ---- the percentage of samples that an OTU must be present 
  #                      in to be considered in the analysis; default is 5% (0.05)
  # tol ---------------- tolerance: should environmental tolerance also be 
  #                      calculated? Default: False
  # verbose ------------ if true the user will get verbose messages. default: False
  
  # First: Tell the user what they've asked for:
  writeLines(paste0("Calculating environmental optima for ", column_to_calc_opt))
  
  # Second: Filter out samples that are NA for this environmental gradient
  # For sanity because testing tibble column types is weird, convert env_data_table
  # to data frame
  env_data_table <- env_data_table %>% as.data.frame()
  
  # Do some sanity checks on the column to make sure it's numeric
  if(!is.numeric(env_data_table[,column_to_calc_opt])) {
    stop("The column you specified is not numeric")
  }
  
  # Count the number of samples that are NA that you will remove
  num_nas <- sum(is.na(env_data_table[,column_to_calc_opt]))
  
  # Tell the user what you've found
  if(num_nas > 0) {
    if(verbose) {
      writeLines(paste0(num_nas, " samples out of ", 
                        length(env_data_table[,column_to_calc_opt]),
                        " are NA and will be filtered out"))
    }
    
    prop_samples_retained <- 1 - (num_nas/length(env_data_table[,column_to_calc_opt]))
    nsamples_remaining <- length(env_data_table[,column_to_calc_opt]) - num_nas
    
    # If there are either less than 5% of samples remaining or fewer than 10 
    # samples, we will not bother calculating the optima for this variable
    if( prop_samples_retained < 0.05 | nsamples_remaining < 10 ) {
      writeLines("Too few samples to calculate optima; returning null.")
      return(NULL)
    }
  }
  
  # Select sample id column and variable to calculate optima for 
  # and filter it to remove NAs
  env_data_table_filt <- env_data_table %>%
    select(all_of(c(sample_id_column, column_to_calc_opt))) %>%
    filter_at(vars(all_of(column_to_calc_opt)), any_vars(!is.na(.))) %>%
    column_to_rownames(var = sample_id_column) %>% as.data.frame()
  
  # Tell user how many sampels were filtered out
  perc_samples_remaining <- round((100-(100*num_nas/length(env_data_table[,column_to_calc_opt]))), 
                              digits = 2)
  if(verbose) {
    writeLines(paste0(perc_samples_remaining, 
                      "% of samples were used to calculate ASV optima"))
  }
  
  # Now that you've removed samples, make sure that otu table and env data 
  # samples are the same
  otu_table <- otu_table[,rownames(env_data_table_filt)]
  
  # Third: Filter to keep otus that are in more than OTU_ubiq_cutoff percent of 
  # samples
  # 
  otu_table_pa <- otu_table # convert to presence absence
  otu_table_pa[otu_table > 0] <- 1 # convert to presence absence
  
  # Get number of samples an OTU must be present in
  nsample_cutoff <- OTU_ubiq_cutoff * ncol(otu_table_pa)
  
  # Filter OTU table
  otu_table_filt <- otu_table[rowSums(otu_table_pa) > nsample_cutoff,]
  
  if(verbose) {
    writeLines(paste0((nrow(otu_table) - nrow(otu_table_filt)), " OTUs were filtered",
                      " out due to not being present in more than ", OTU_ubiq_cutoff*100,
                      "% of samples"))
  }
  
  # transform OTU table to sample-by-species table
  otu_table_t <- t(otu_table_filt)
  
  # Tell the user the number of OTUs and samples being used
  writeLines(paste0("The environmental optima for ", column_to_calc_opt, " will",
                    " be calculated for ", ncol(otu_table_t), " organisms using ",
                    nrow(otu_table_t), " samples."))

  # Calculate the environmental optima  
  
  # Convert environmental variable to named numeric vector
  env <- env_data_table_filt[, column_to_calc_opt]
  names(env) <- rownames(env_data_table_filt)

  env_optima <- analogue::optima(otu_table_t, 
                                 env)
  
  if(tol) {
    env_tol <- tolerance(otu_table_t, env, 
                         useN2 = TRUE)
    env_opt_tol.df <- data.frame(env_optima, tolerance = as.numeric(env_tol))
    return(env_optima = env_opt_tol.df)
  }
  
  env_optima.df <- as.data.frame(env_optima)
  names(env_optima.df) <- paste0(column_to_calc_opt, "_opt")
  return(env_optima = env_optima.df)
}

# Generate mantel correlogram for an optima distance matrix and phylogenetic 
# distance matrix
calculate_mantel_correlogram <- function(env_opt_dm = env_opt_distmat$all, 
                                            pat_dist_norm = pat_dist_normalized,
                                            n_dist_class = 50) {
  # env_opt_dm ---- The distance matrix for environmental optima.  Expects 'dist' 
  #                 object of pairwise comparisons between ASVs. default: 
  #                 env_opt_distmat[[1]] 
  # pat_dist_norm - The distance matrix between tips of the phylogenetic tree. 
  #                 This will be subset to include only tips that exist in 
  #                 env_opt_dm; It should be normalized to scale from 0-1; 
  #                 default: pat_dist_normalized
  # n_dist_class -- the number of distance classes to feed to mantel.correlog
  
  # Filter patristic distance matrix
  rownms <- rownames(as.matrix(env_opt_dm))
  pat_dist_norm_filt <- pat_dist_norm[rownms, rownms]
  
  mcorr <- mantel.correlog(env_opt_dm, pat_dist_norm_filt, 
                  n.class = n_dist_class, r.type="pearson", cutoff = FALSE)
  return(mcorr)
}

# Plot tree
plot_optima <- function(tree = input_ra$tree,
                        taxonomy = input_ra$taxonomy,
                        env_optima = env_optima_list[[1]],
                        opt_column = "DepthAvg___opt",
                        filter_tree = TRUE) {
  # Requires:
  library(tidytree)
  library(ape)
  library(ggtree)
  library(ggnewscale)
  library(ggstance)
  # First drop any tips that are not in the module list
  env_optima <- env_optima %>% 
    rownames_to_column(var = "genome") %>%
    left_join(taxonomy, by = "genome")
  
  if(filter_tree) {
    tree.filt <- keep.tip(tree, env_optima$genome)
    
    Ntip(tree.filt) == nrow(env_optima) # check filtering successful; Should be true
  } else {
    tree.filt <- tree
  }
  
  
  # Next, add the module membership to the tree
  tib_tree <- as_tibble(tree.filt)
  
  tib_tree_data.circ <- left_join(tib_tree, 
                                  env_optima,
                                  by = c("label" = "genome")) %>%
    as.treedata()
  
  # Heatmap for circular tree
  heatmap <- env_optima %>% 
    select(genome, all_of(opt_column)) %>% 
    column_to_rownames(var = "genome")
  
  # Plot vertical tree
  p <- ggtree(tree.filt)
  pdat <- p
  
  # module_palette_circ <- 1:length(module_palette)
  # names(module_palette_circ) <- module_palette
  
  # Plot circular tree
  #### =============== ####
  # Generate clade_labels for phyla
  phylumlist <- env_optima$Phylum %>%
    unique()
  classlist <- env_optima %>%
    filter(Phylum == "p__Proteobacteria") %>%
    filter(Class != "NA") %>%
    select(Class) %>%
    unique()
  cladelist <- c(phylumlist[phylumlist!="p__Proteobacteria"], classlist$Class)
  
  # remove NA
  cladelist <- cladelist[!is.na(cladelist)]
  clade_labels <- rep(NA, length = length(cladelist))
  
  # IMPORTANT: This only labels phylum with more than 1 genomes in them. There are
  #  phyla that have fewer ASVs that are unlabelled (this helps with the spacing of the phylum labels on the tree).
  # For my own knowledge: "M.stor:PLGY01", "A.stor:20170700_S25_26"
  for (i in 1:length(cladelist)) {
    # if proteobacteria
    if( grepl("proteobacteria", cladelist[i]) ) {
      tip_list <- env_optima %>%
        filter(Class == cladelist[i])
    } else { # not proteobacteria
      tip_list <- env_optima %>%
        filter(Phylum == cladelist[i])
    }
    # Debugging
    #print(paste0("Group ", cladelist[i]))
    #print(length(tip_list$genome))
    # Get the most recent common ancestor of clades
    clade_labels[i] <- ifelse(length(tip_list$genome) > 1, getMRCA(tree.filt,
                                                                   tip = tip_list$genome),NA)
  }
  names(clade_labels) <- cladelist
  
  clade_labels <- na.omit(clade_labels)
  clade_labels <- sort(clade_labels)
  
  # make colors for phyla
  clade_labels_color <- rep(c("grey40", "grey20", "grey80"), 
                            length.out = length(clade_labels))
  #clade_labels_color <- c("red", "orange", "yellow", "green", "blue", "purple", "pink", "black", "grey40", "grey20")
  
  # Plot base tree
  p <- ggtree(tib_tree_data.circ, # layout = "rectangular") +
              layout = "fan", open.angle = 5) +
    geom_point2(aes(subset=(node==1209)), # this is the archaea
                shape=21, size=2, fill='red') +
    geom_tippoint(aes(subset=(label %in% c("PLGY01", "20170700_S25_26"))),
                  shape = 21, size = 2, fill = "yellow")

  # p
  
  # Annotate clades
  p1 <- p
  
  for (i in 1:length(clade_labels)) {
    p1 <- p1 + geom_cladelabel(node = clade_labels[i],
                               label = names(clade_labels)[i],
                               align = TRUE, barsize = 2,
                               offset = 1,
                               offset.text = 0.3,
                               fontsize = 3,
                               hjust = 0,
                               #alpha = 0.5,
                               angle = "auto",
                               color = "grey20")
    #color = clade_labels_color[i])
  }
  # Troubleshooting
  # p1 +
  #   geom_nodelab(aes(label = node)) +
  #   geom_hilight(node = 80 )#+
  #geom_tiplab(aes(subset = Class == "c__Alphaproteobacteria"),
  #            align = TRUE, offset = rel(0.2))
  
  # Plot Heatmap
  p2 <- p1 + new_scale_fill()
  p3 <- gheatmap(p2, heatmap, offset = .2,
                 width = .2, colnames = FALSE) +
    scale_fill_viridis_c(name = opt_column, 
                      na.value = "white") +
    xlim_tree(1.5)
  # Debugging
  p3 
  
  return(p3)
}

#### ====================================================================== ####
# Step 1: Subset environmental data to only continuous values (that can be tested)
# and include only those that make sense to test (for example month and year should
# be excluded)
#### ====================================================================== ####
env_data_table_num <- input_counts$sample_metadata %>%
  column_to_rownames(var = "temporal_sample_id") %>%
  select_if(is.numeric) %>%
  select(!matches(c("Core__", "Month__", "Year__", "DepthMin__", "DepthMax__",
                   # WTD and ALD have artificial low values where value was below detection limit
                    "WTD", "ALD",
                   # Remove CO2.percent and CH4.percent; they have very few values
                   # and for CH4.percent all values are 0
                   "CO2.percent", "CH4.percent",
                   # Remove CH4/CO2 ratio/percent, and fluxes, they are derived values with many NAs
                   # and the nested nature of the fluxes is not great either
                   "wk_before", "CH4_CO2_porewater",
                   # Remove non-relevant AirTemperature columns
                   "samplingdate", "^n_AirTemperature", "sd_AirTemperature",
                   "median_AirTemperature", "mean_AirTemperature",
                  "UpdateDate__.DepthInfo", "UpdateDate__.Biogeochemistry"))) %>%
  # The number of NA values must be less than 85% of total observations
  select_if(~sum(is.na(.)) < 179)

# Print the names of the columns
names(env_data_table_num)
# Checkout variable distributions
# Checkout optima distributions
lapply(2:ncol(env_data_table_num), 
       function(x) {
         env_hist.df <- env_data_table_num[,c(1,x)]
         env.hist.plot <- env_hist.df  %>%
           ggplot(aes(x = .data[[names(env_hist.df)[2]]])) + 
           geom_histogram()
  
  # Create directory for plots if it doesn't exist
  if(!dir.exists(here(figures.fp, "env_var_histograms"))) {
    dir.create(here(figures.fp, "env_var_histograms"))}
  
  # Save plots
  plot_filename <- here(figures.fp, "env_var_histograms", 
                        paste0(names(env_data_table_num[x]), "_hist.png"))
  # fix special characters
  plot_filename <- gsub("%", "Perc", plot_filename)
  
  ggsave(env.hist.plot, file = plot_filename)
})

#### ====================================================================== ####

# Step 2: Approximate the niche optima for each otu across each environmental 
# data. An otu is considered to show an optimum at 
# the environmental variable where it most often occurs
#### ====================================================================== ####
numeric_var_names <- names(env_data_table_num)
numeric_var_names

env_data_table_num <- env_data_table_num %>%
  rownames_to_column(var = "temporal_sample_id")

env_optima_list <- lapply(numeric_var_names, function(x) {
  calc_niche_optima(otu_table = input_counts$otu_table[,-1],
                    env_data_table = env_data_table_num,
                    sample_id_column = "temporal_sample_id",
                    column_to_calc_opt = x,
                    OTU_ubiq_cutoff = 0.05,
                    tol = FALSE)})
# add names 
names(env_optima_list) <- numeric_var_names # add names to list

# Drop null values from list:
env_optima_list <- env_optima_list[ !sapply(env_optima_list,is.null) ]

# Create a dataframe with all optima
all_opt <- lapply(env_optima_list, function(x) {x$genome <- rownames(x); return(x)}) %>%
  reduce(full_join, by = "genome") %>%
  select(genome, everything()) %>%
  filter_all(all_vars(!is.na(.))) %>%
  column_to_rownames(var = "genome")

# Add all optima df to the list
env_optima_list$all <- all_opt

# Checkout optima distributions
lapply(seq_along(env_optima_list), function(x) {
  opt_hist.df <- env_optima_list[[x]] %>% 
    rownames_to_column(var = "genome")
  opt.hist.plot <- opt_hist.df  %>%
    ggplot(aes(x = .data[[names(opt_hist.df)[2]]])) + 
      geom_histogram()
  
  # Create directory for plots if it doesn't exist
  if(!dir.exists(here(figures.fp, "optima_histograms"))) {
    dir.create(here(figures.fp, "optima_histograms"))}
  
  # Save plots
  plot_filename <- here(figures.fp, "optima_histograms", 
                        paste0(names(env_optima_list[x]), "_opt", "_hist.png"))
  # fix special characters
  plot_filename <- gsub("%", "Perc", plot_filename)
  
  ggsave(opt.hist.plot, file = plot_filename)
})

#### ====================================================================== ####

# Step 3: Calculate distance matrix between the env. optima of OTUs
#### ====================================================================== ####

env_opt_distmat <- lapply(seq_along(env_optima_list), function(x) {
  #dist(x, method = "manhattan")
  print(names(env_optima_list[x]))
  opt_df <- env_optima_list[[x]]
  vegdist(opt_df, method = "gower")
})
names(env_opt_distmat) <- names(env_optima_list)
sapply(env_opt_distmat, function(x) {
  list(min(x), max(x))
})

writeLines("Number of OTUS being tested for each variable: ")
env_opt_notu <- sapply(env_opt_distmat, function(x) {
  length(attr(x, "Labels")) # report notus 
  }) 
print(env_opt_notu)
#### ====================================================================== ####

# Step 4: Calculate patristic distance of tree and normalize
#### ====================================================================== ####
pat_dist_all <- cophenetic(input_counts$tree)

# Normalize distances from 0-1; 
# Question for later: Should this be done on the master tree of 1662 tips or 
# on the sub trees which have subsets of the 1662 tips?; Currently we do this on 
# the whole tree which would make the phylogenetic distances comparable even when
# only a subset of OTUs have been used to calculate optima
max_pat_dist = max(pat_dist_all)
min_pat_dist = min(pat_dist_all)
pat_dist_normalized <- apply(pat_dist_all, MARGIN = c(1,2), 
                             FUN = function(x) {((x-min_pat_dist)/(max_pat_dist-min_pat_dist))})

#### ====================================================================== ####

# Step 5: Calculate Mantel Correlograms
#### ====================================================================== ####
# use simple forked cluster mclapply; won't work on windows systems, and can't 
# use multiple nodes
ncore <- detectCores()
writeLines(paste0("using ", ncore, " cores, starting loop."))
mcorr_list <- mclapply(seq_along(env_opt_distmat), function(x) {
  var_name <- names(env_opt_distmat)[[x]]
  # Debugging
  print(var_name)
  
  # reorder to match
  dim(pat_dist_normalized)
  dist_mat <- as.matrix(env_opt_distmat[[x]])
  dim(dist_mat)
  pat_dist_filt <- pat_dist_normalized[rownames(dist_mat), colnames(dist_mat)]
  dim(pat_dist_filt)
  # print(rownames(dist_mat) == rownames(pat_dist_filt))
  
  mcorr <- calculate_mantel_correlogram(env_opt_dm = as.dist(dist_mat),
                               pat_dist_norm = pat_dist_filt,
                               n_dist_class = 50)
  saveRDS(mcorr, paste0(outputs.fp, "/mcorr_", var_name,".RDS"))
  }, mc.cores = ncore)
# 
# data.frame(Pat_dist = as.vector(pat_dist_filt), 
#            +            dist_mat = as.vector(dist_mat)) %>% View()
#   ggplot(aes(x = Pat_dist, y = dist_mat)) + geom_point()

names(mcorr_list) <- names(env_opt_distmat)
#### ====================================================================== ####

# Step 6: Plot correlograms
# load from server output
#### ====================================================================== ####
writeLines("Plotting Correlelograms")
mcorr_files <- list.files(outputs.fp, pattern = "mcorr_", full.names = T)
mcorr_list <- lapply(mcorr_files, 
       function(x) {readRDS(x)})

mcorr_names <- sapply(mcorr_files, function(x) {
  y <- gsub(paste0(outputs.fp, "/mcorr_"), "", x)
  z <- gsub(".RDS", "", y)
  return(z)
})

names(mcorr_list) <- mcorr_names

correlogram.plot.df <- lapply(seq_along(mcorr_list), function(x) {
  mantel.res <- mcorr_list[[x]]$mantel.res %>% 
    as.data.frame() %>%
    mutate(EnvVar = names(mcorr_list[x]))
  return(mantel.res)}) %>%
  reduce(bind_rows) %>%
  mutate(signif = ifelse(`Pr(corrected)` < 0.05, TRUE, FALSE)) %>%
  left_join(data.frame(NOTU = env_opt_notu) %>%
              rownames_to_column(var = "EnvVar"), by = "EnvVar") %>%
  mutate(EnvVarFacet = paste0(EnvVar, " (MAGs = ", NOTU, ")"))

mc.plot <- ggplot(data = correlogram.plot.df, 
                  aes(x = class.index, y = Mantel.cor)) +
  geom_point(aes(shape = signif), color = "black") +
  scale_shape_manual(values = c(0,15)) +
  #scale_color_discrete(name = "Environmental Gradient") +
  geom_line(color = "black") +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  facet_wrap(~EnvVarFacet) + 
  theme_bw() +
  xlab("Phylogenetic Distance (normalized)") + 
  ylab("Mantel correlation") +
  guides(color=guide_legend(title="Environmental Gradient"),
         shape = "none")
mc.plot  

mc.plot_size <- ggplot(data = correlogram.plot.df, 
                  aes(x = class.index, y = Mantel.cor)) +
  geom_point(aes(color = EnvVar, shape = signif, size = n.dist)) +
  scale_shape_manual(values = c(0,15)) +
  #scale_color_discrete(name = "Environmental Gradient") +
  geom_line(aes(color = EnvVar)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  facet_wrap(~EnvVarFacet) + 
  theme_bw() +
  xlab("Phylogenetic Distance (normalized)") + 
  ylab("Mantel correlation") +
  guides(color=guide_legend(title="Environmental Gradient"),
         shape = "none")
mc.plot_size

#### ====================================================================== ####

# Attempting plot signifiant correlations on phylogenetic tree
#### ====================================================================== ####
# Plot trees filtered by OTUs
ncore <- detectCores()
writeLines(paste0("detected ", ncore, " cores, starting loop."))

env_opt_tree_plots <- mclapply(seq_along(env_optima_list), function(x) {
  writeLines(names(env_optima_list[x]))
  plot <- plot_optima(tree = input_ra$tree,
              taxonomy = input_ra$taxonomy,
              env_optima = env_optima_list[[x]],
              opt_column = paste0(names(env_optima_list[x]), "_opt"),
              filter_tree = TRUE)

  # Create directory for plots if it doesn't exist
  if(!dir.exists(here(figures.fp, "optima_plots"))) {
    dir.create(here(figures.fp, "optima_plots"))}
  
  # Save plots
  plot_filename <- here(figures.fp, "optima_plots", 
                        paste0(names(env_optima_list[x]), "_opt", "_filttree.png"))
  # fix special characters
  plot_filename <- gsub("%", "Perc", plot_filename)
  
  ggsave(plot, file = plot_filename)
  return(plot)
}, mc.cores = ncore)

# Plot trees filtered by OTUs
env_opt_tree_plots_allotu <- mclapply(seq_along(env_optima_list), function(x) {
  writeLines(names(env_optima_list[x]))
  plot <- plot_optima(tree = input_ra$tree,
                      taxonomy = input_ra$taxonomy,
                      env_optima = env_optima_list[[x]],
                      opt_column = paste0(names(env_optima_list[x]), "_opt"),
                      filter_tree = F)
  
  # Create directory for plots if it doesn't exist
  if(!dir.exists(here(figures.fp, "optima_plots"))) {
    dir.create(here(figures.fp, "optima_plots"))}
  
  # Save plots
  plot_filename <- here(figures.fp, "optima_plots", 
                        paste0(names(env_optima_list[x]), "_opt", "_tree.png"))
  # fix special characters
  plot_filename <- gsub("%", "Perc", plot_filename)
  
  ggsave(plot, file = plot_filename)
  return(plot)
}, mc.cores = ncore)


#### ====================================================================== ####

#### Save Data and Figures
#### ====================================================================== ####

# Data
# Mantel Correlogram List 
saveRDS(mcorr_list,
          paste0(outputs.fp, "/mcorr_list.RDS"))

# Figures
ggsave(mc.plot, device = "png", width = 14, height = 7,
       filename = paste0(figures.fp, "/phylogentic_dist_mantel_correlogram.png"))

ggsave(mc.plot_size, device = "png", width = 14, height = 7,
       filename = paste0(figures.fp, "/phylogentic_dist_mantel_correlogram_size.png"))