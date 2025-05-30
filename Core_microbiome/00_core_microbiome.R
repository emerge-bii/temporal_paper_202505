#' ## Defining the core microbiome as recommended in Shade and Stopnisek 2019: https://www.sciencedirect.com/science/article/pii/S1369527419300426 
#' 

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
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("Core_microbiome", "outputs")
figures.fp <- here("Core_microbiome", "figures")

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
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), 
                                 c("#737373", "#FFFFFF", "#000000")), 
                          c("blue", "darkblue", "green", "darkgreen", "red", 
                            "darkred", "orange", "darkorange", "purple", 
                            "darkpurple", "yellow", "brown", "grey", "white", 
                            "black"))
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

# Define some helper functions
#### ====================================================================== ####
calculate_bc_core_contribution_time <- function(otu = fen_otus[,-1], 
                                     map = input_counts$sample_metadata, 
                                     n_top_otus = 500, display_plots = T, grp_fct = "DepthLumping") {
  
  # Calculate occurance, presence absence, relative abundance and combine occupancy and abundance
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # average relative abundance across all samples of each OTU  
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('OTU_ID')
  
  if(display_plots) {
    # Occupancy abundance plot:
    ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
      geom_point(pch=21, fill='white') +
      labs(x="log10(mean relative abundance)", y="Occupancy")
  }
  
  
  PresenceSum <- tibble(OTU_ID = as.factor(row.names(otu)), otu) %>% 
    pivot_longer(-OTU_ID, names_to = "temporal_sample_id", values_to = "abun") %>%
    left_join(map, by = 'temporal_sample_id') %>%
    group_by(OTU_ID, !!as.name(grp_fct)) %>% # we want the core that persists over time therefore group by depth, NOT by time
    summarise(grp_freq=sum(abun>0)/length(abun),        # frequency of detection between time points (# of times otu detected/ # of observations for each time point)
              coreGrp=ifelse(grp_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific time (all samples of a timepoint have OTU in them), 0 if not
    group_by(OTU_ID) %>% 
    summarise(sumF=sum(grp_freq), # sums frequency of detection, frequency is at most the number of time points (7)
              sumG=sum(coreGrp), # sums if an OTU is found in all points at a given time (binary); max is 7 
              nS=length(!!as.name(grp_fct)), # number of years
              Index=(sumF+sumG)/(nS*2)) # calculating weighting Index based on number of time points detected and frequency 
                                    # If OTU present all the time in all samples = 1 otherwise <1; 2 is there to scale between 0 and 1
                              
  if(display_plots) {
    # Plot weighing index
    ggplot(PresenceSum, aes(y = sumG, x = sumF)) +
      geom_point(pch=21, aes(fill= Index)) +
      scale_fill_viridis(name = "Weighing index") +
      labs(x=paste0("Frequency within ", grp_fct), y=paste0("Occupancy within ", grp_fct))
  }
  
  # Create the ranking list of OTUS
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by='OTU_ID') %>%
    transmute(OTU_ID=OTU_ID,                           
              rank=Index) %>%
    arrange(desc(rank))
  
  
  # Calculate contribution of each OTU to BC dissimilarity; start with top OTU, then calculate until n_top_otus reached
  
  # Calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL
  
  # calculating BC dissimilarity based on the 1st ranked OTU
  ReadDepth <- colSums(otu) # get read depths for each sample (for caluclating bray curtis)
  otu_start <- otu_ranked$OTU_ID[1]                   
  start_matrix <- as.matrix(as.matrix(otu)[otu_start,]) # OTU table with only highest-ranked otu
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(ReadDepth[x[1]] + ReadDepth[x[2]])) # bray-curtis based on single OTU for each sample comparison
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  # Now that the top OTU is finished, do it for the 2nd through n_top_otus otus
  # calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th. 
  # Can be set to the entire length of OTUs in the dataset, however it might take 
  # some time if more than 5000 OTUs are included.
  for(i in 2:n_top_otus){                              
    otu_add=otu_ranked$OTU_ID[i]                       
    add_matrix <- as.matrix(as.matrix(otu)[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix) # Use both previous matrix and add new otus for bray-curtis calculation
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(ReadDepth[x[1]] + ReadDepth[x[2]]))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  
  # calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
  if(n_top_otus < nrow(otu_ranked)) {
    x <-  apply(combn(ncol(as.matrix(otu)), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(ReadDepth[x[1]] + ReadDepth[x[2]]))   
    x_names <- apply(combn(ncol(as.matrix(otu)), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
    df_full <- data.frame(x_names,x)
    names(df_full)[2] <- length(rownames(otu))
    BCfull <- left_join(BCaddition,df_full, by='x_names')
  } else {
    BCfull <- BCaddition
  }
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  temp_BC_matrix_t <- as_tibble(t(temp_BC_matrix))
  
  BC_ranked <- tibble(rank = as.factor(row.names(t(temp_BC_matrix))),temp_BC_matrix_t) %>% 
    pivot_longer(-rank, names_to = "comparison", values_to = "BC") %>%
    group_by(rank) %>%
    summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
    arrange(desc(-MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
  Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  
  if(n_top_otus < nrow(otu_ranked)) {
    BC_ranked <- BC_ranked[-nrow(BC_ranked),] # removes the final value if not all OTUS are included in n_top_otus 
  }
  
  BC_contribution <- bind_cols(BC_ranked, otu_ranked[1:n_top_otus,] %>% rename(raw_rank = rank)) %>%
    select(OTU_ID, rank, raw_rank, MeanBC, proportionBC, IncreaseBC)
  
  return(list(BC_ranks = BC_contribution, occ_abun = occ_abun))
}

# Creates thresholds for core inclusion
create_core_thresholds <- function(BC_ranks = BC_contribution, 
                                   method = "Both", display_plots = T) {
  
  # Method can be "A", "B", or "Both"
  if(method %in% c("A", "Both")) {
    #A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
    fo_difference <- function(pos){
      left <- (BC_ranks[pos, "MeanBC"] - BC_ranks[1, "MeanBC"]) / pos
      right <- (BC_ranks[nrow(BC_ranks), "MeanBC"] - BC_ranks[pos,"MeanBC"]) / (nrow(BC_ranks) - pos)
      return(left - right)
    }
    BC_ranks$fo_diffs <- unlist(sapply(1:nrow(BC_ranks), fo_difference))
    
    elbow <- which.max(BC_ranks$fo_diffs)
  }
  
  if(method %in% c("B", "Both")) {
    #B) Final increase in BC similarity of equal or greater then 2% 
    lastCall <- last(as.numeric(as.character(BC_ranks$rank[(BC_ranks$IncreaseBC>=1.02)])))
  }
  
  if(display_plots) {
    #Creating plot of Bray-Curtis similarity
    p <- BC_ranks[1:nrow(BC_ranks),] %>%
      mutate(rank = factor(rank, levels=c(1:n()))) %>%
      ggplot(aes(x=rank)) +
      geom_point(aes(y=proportionBC)) +
      theme_classic() + 
      theme(strip.background = element_blank(),
            axis.text.x = element_text(size=7, angle=45)) +
      labs(x='ranked OTUs',y='Bray-Curtis similarity')
    
    if(method %in% c("A", "Both")) {
      p <- p +
        geom_vline(xintercept=elbow, linewidth=3, col='red', cex=.5) +
        annotate(geom="text", x=elbow+14, y=.1, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")
    }
    
    if(method %in% c("B", "Both")) {
      p <- p + 
        geom_vline(xintercept=lastCall, linewidth=3, col='blue', cex=.5) +
        annotate(geom="text", x=lastCall+3, y=.5, label=paste("Last 2% increase (",lastCall,")",sep=''), color="blue")
    }
      
    p  
  }
  
  # Return list construction
  if(method == "A") {
    return_list <- list(elbow = elbow, lastCall = NULL)
  }
  if(method == "B") {
    return_list <- list(elbow = NULL, lastCall = lastCall)
  }
  if(method == "Both") {
    return_list <- list(elbow = elbow, lastCall = lastCall)
  }
  
  # Add plots if there
  if(display_plots) {
    return_list <- append(return_list, list(plot = p))
  }
  
  return(return_list)
}

apply_sloan_neutral_model <- function(otu, map, BC_ranks = BC_contribution, 
                                      occ_abun, cutoffs, grp_fct = "DepthLumping",
                                      fill_option = "None", group_colors = colour_depth,
                                      taxonomy = NULL) {
  # Source function list
  source(here("Core_microbiome", "sncm.fit_metagenome.R")) # Source the sncm.fit functions
  
  # Get core cutoffs
  elbow <- cutoffs$elbow
  lastCall <- cutoffs$lastCall
  
  # Use Sloan neutral model to prioritize OTUs
  # Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
  spp <- t(otu)
  taxon=as.vector(rownames(otu))
  
  if(!is.null(taxonomy)) {
    taxon <- data.frame(genome = taxon) %>%
      left_join(taxonomy, by = "genome") %>%
      column_to_rownames("genome")
  }
  
  # HANNAH - figure out how to deal with non-integer numbers and sequencing depth differences
  #Models for the whole community
  obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL, whichN = "mean")
  sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL, whichN = "mean")
  
  above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs above prediction 
  below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs below prediction
  
  #Create a column defining "core" OTUs
  occ_abun$fill <- 'no'
  if(!is.na(lastCall)) {
    occ_abun$fill[occ_abun$OTU_ID %in% BC_ranks$OTU_ID[1:lastCall]] <- 'core'
  } else {
    writeLines("2 % increase not found, using elbow method instead")
    occ_abun$fill[occ_abun$OTU_ID %in% BC_ranks$OTU_ID[1:elbow]] <- 'core'
  }
  
  
  if(fill_option == "None") {
    point_fill <- "black"
  } else {
    point_fill <- colour_habitat[which(habitat_levels == fill_option)]
  }
  
  model_plot <- ggplot() +
    geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
    geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill=point_fill, size=1.8) +
    geom_line(color='black', data=obs.np, linewidth=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
    annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label = paste0("   ","Sloan Neutral Model\n",
                            "   ","R2 = ", round(sta.np$Rsqr, 2), "\n",
                            "   ","AIC = ", round(sta.np$AIC, 2), "\n",
                            "   ","m = ", round(sta.np$m, 2), " +/- ", round(sta.np$m.ci, 2))) +
    labs(x="log10(mean relative abundance)", y="Occupancy")
  model_plot
  
  binom_model_plot <- ggplot() +
    geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
    geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill=point_fill, size=1.8) +
    geom_line(color='black', data=obs.np, linewidth=1, aes(y=obs.np$bino.pred, x=log10(obs.np$p)), alpha=.25) +
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$bino.upr, x=log10(obs.np$p)), alpha=.25)+
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$bino.lwr, x=log10(obs.np$p)), alpha=.25)+ 
    annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label = paste0("   ","Binomial Model\n",
                            "   ","R2 = ", round(sta.np$Rsqr.bino, 2), "\n",
                            "   ","AIC = ", round(sta.np$AIC.bino, 2))) +
    labs(x="log10(mean relative abundance)", y="Occupancy")
  binom_model_plot
  
  poisson_model_plot <- ggplot() +
    geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
    geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill=point_fill, size=1.8) +
    geom_line(color='black', data=obs.np, linewidth=1, aes(y=obs.np$pois.pred, x=log10(obs.np$p)), alpha=.25) +
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$pois.upr, x=log10(obs.np$p)), alpha=.25)+
    geom_line(color='black', lty='twodash', linewidth=1, data=obs.np, aes(y=obs.np$pois.lwr, x=log10(obs.np$p)), alpha=.25)+
    annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label = paste0("   ","Poisson Model\n",
                            "   ","R2 = ", round(sta.np$Rsqr.pois, 2), "\n",
                            "   ","AIC = ", round(sta.np$AIC.pois, 2))) +
    labs(x="log10(mean relative abundance)", y="Occupancy")
  poisson_model_plot
  
  
  #' Exercise 2:
  #' Highlight on the same occupancy-abundance plot OTUs that are above and below  
  #' the natural model prediction.
  #' Extra: add to the plot the above.pred and below.pred values
  
  
  #Creating a plot of core taxa occupancy by time point
  core <- occ_abun$OTU_ID[occ_abun$fill == 'core']
  
  otu_relabun <- decostand(otu, method="total", MARGIN=2)
  
  plotDF <- tibble(OTU_ID = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
    pivot_longer(names_to = "temporal_sample_id", values_to = "relabun", -OTU_ID) %>%
    left_join(map, by = 'temporal_sample_id') %>%
    left_join(BC_ranks, by='OTU_ID') %>%
    filter(OTU_ID %in% core) %>% 
    group_by(OTU_ID, !!as.name(grp_fct)) %>%
    summarise(grp_freq=sum(relabun>0)/length(relabun),        
              coreGrp=ifelse(grp_freq == 1, 1, 0),      
              detect=ifelse(grp_freq > 0, 1, 0))
  
  plotDF$OTU_ID <- factor(plotDF$OTU_ID, levels=BC_ranks$OTU_ID[1:length(core)])
  
  
  # Plot occupancy by group
  occup_core_plot <- ggplot(plotDF,aes(x=OTU_ID, grp_freq)) +    
    geom_bar(aes(fill=factor(!!as.name(grp_fct))),
    stat = 'identity', position = position_dodge(-.9)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(plotDF$OTU_ID))) +
    theme(axis.text = element_text(size=6)) +
    labs(x='Ranked OTUs', y=paste0('Occupancy by ', grp_fct), fill=grp_fct) +
    #scale_fill_manual(values = colour_year, breaks = year_levels) +
    theme_bw()
  
  # Add colors to plot if provided
  if(!is.null(group_colors)) {
    occup_core_plot <- occup_core_plot + 
      scale_fill_manual(values = group_colors)
  }
  
  occup_core_plot
  
  return(list(model_plot = model_plot, binom_model_plot = binom_model_plot, 
              poisson_model_plot = poisson_model_plot,
              occup_core_plot = occup_core_plot, 
              model_comp_stats = sta.np,
              model_pred_output = obs.np,
              occ_abun_core = occ_abun))
}

#### ====================================================================== ####


# For Palsa
#### ====================================================================== ####
# Prepare Palsa Data
habitat <- "Palsa"
palsa_samples <- input_counts$sample_metadata %>% filter(Habitat__ == "Palsa") %>% 
  select(temporal_sample_id)
palsa_otus <- input_counts$otu_table %>%
  select(genome, any_of(palsa_samples$temporal_sample_id))


palsa_core_bc <- calculate_bc_core_contribution_time(otu = palsa_otus[,-1],
                                                     map = input_counts$sample_metadata,
                                                     n_top_otus = 1662, display_plots = T)

palsa_core_th <- create_core_thresholds(BC_ranks = palsa_core_bc$BC_ranks)
palsa_core_th$plot

ggsave(file = paste0(figures.fp, "/core_th_", habitat, ".png"),
       plot = palsa_core_th$plot, device = "png")


palsa_neut_model <- apply_sloan_neutral_model(otu = palsa_otus[,-1],
                                              map = input_counts$sample_metadata,
                                              BC_ranks = palsa_core_bc$BC_ranks,
                                              occ_abun = palsa_core_bc$occ_abun,
                                              cutoffs = palsa_core_th, 
                                              fill_option = "Palsa")
palsa_neut_model$occup_core_plot
palsa_neut_model$model_plot
palsa_neut_model$binom_model_plot
palsa_neut_model$poisson_model_plot
palsa_neut_model$model_pred_output %>% View()


# Join occ_abund and BC ranks dataframe
palsa_core_otu_table <- left_join(palsa_core_bc$BC_ranks, palsa_neut_model$occ_abun_core, by = "OTU_ID") %>%
  rename(AverageRelativeAbundance = otu_rel,
         OccupancyProp = otu_occ,
         CoreStatus = fill)

write.table(palsa_core_otu_table, file = paste0(outputs.fp, "/core_otus_table_", habitat, ".txt"),
            sep = "\t", row.names = F)

# Save the occupancy core plot
ggsave(plot = palsa_neut_model$occup_core_plot,
       file = paste0(figures.fp, "/occup_core_", habitat, ".png"),
       device = "png")
ggsave(plot = palsa_neut_model$model_plot,
       file = paste0(figures.fp, "/neutral_model_", habitat, ".png"),
       device = "png")
ggsave(plot = palsa_neut_model$binom_model_plot,
       file = paste0(figures.fp, "/bionom_model_", habitat, ".png"),
       device = "png")
ggsave(plot = palsa_neut_model$poisson_model_plot,
       file = paste0(figures.fp, "/poisson_model_", habitat, ".png"),
       device = "png")


# # Artificially rarefying appears to give much the same results as not
# palsa_otus_rarefied <- input_counts$otu_table %>%
#   select(genome, any_of(palsa_samples$temporal_sample_id)) %>%
#   mutate(across(where(is.numeric), ~floor(10*.))) %>% # use floor rather than round so we simply make the detection limit more corse rather than messing with the ; can adjust the detection limit by multiplying column by factor of 10
#   select_if(negate(function(col) is.numeric(col) && sum(col) < 3084)) %>% #dim()
#   select(-genome) %>% #colSums() %>% sort() %>% hist()
#   t() %>%
#   rrarefy(., sample = 3084) %>% t() %>%
#   as_tibble(rownames = "genome") %>%
#   #mutate(across(where(is.numeric), ~./10)) %>%
#   mutate(genome2 = genome) %>%
#   column_to_rownames("genome2")
# 
# 
# 
# 
# palsa_core_bc_rare <- calculate_bc_core_contribution_time(otu = palsa_otus_rarefied[,-1],
#                                                      map = input_counts$sample_metadata,
#                                                      n_top_otus = 1662, display_plots = T)
# 
# palsa_core_th_rare <- create_core_thresholds(BC_ranks = palsa_core_bc_rare$BC_ranks)
# palsa_core_th_rare$plot
# 
# palsa_neut_model <- apply_sloan_neutral_model(otu = palsa_otus_rarefied[,-1],
#                                               map = input_counts$sample_metadata,
#                                               BC_ranks = palsa_core_bc_rare$BC_ranks,
#                                               occ_abun = palsa_core_bc_rare$occ_abun,
#                                               cutoffs = palsa_core_th_rare, 
#                                               fill_option = "Palsa")
# palsa_neut_model_rare$occup_core_plot
# palsa_neut_model_rare$model_plot
# palsa_neut_model_rare$binom_model_plot
# palsa_neut_model_rare$poisson_model_plot

# Dispersal limitation is important to model

# Palsa by depth:
habitat <- "Palsa"
Palsa_list <- lapply(depth_levels, function(x) {
  which_samples <- input_counts$sample_metadata %>% filter(Habitat__ == habitat) %>% 
    filter(DepthLumping == x) %>%
    select(temporal_sample_id)
  which_otus <- input_counts$otu_table %>%
    select(genome, any_of(which_samples$temporal_sample_id))
  
  # Step 1: calculate the core BC of each of the top 500 taxa
  core_bc <- calculate_bc_core_contribution_time(otu = which_otus[,-1],
                                                       map = input_counts$sample_metadata,
                                                       n_top_otus = 1662, display_plots = T)
  
  # Step 2: calculate the core BC of each of the top 500 taxa
  core_th <- create_core_thresholds(BC_ranks = core_bc$BC_ranks)
  
  ggsave(file = paste0(figures.fp, "/core_th_", habitat, "_", x, ".png"),
         plot = core_th$plot, device = "png")
  
  # Step 3: Apply the sloan neutral model
  neut_model <- apply_sloan_neutral_model(otu = which_otus[,-1],
                                                map = input_counts$sample_metadata,
                                                BC_ranks = core_bc$BC_ranks,
                                                occ_abun = core_bc$occ_abun,
                                                cutoffs = core_th, 
                                                fill_option = habitat)
  
  # Save the occupancy core plot
  ggsave(plot = neut_model$occup_core_plot,
         file = paste0(figures.fp, "/occup_core_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$model_plot,
         file = paste0(figures.fp, "/neutral_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$binom_model_plot,
         file = paste0(figures.fp, "/bionom_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$poisson_model_plot,
         file = paste0(figures.fp, "/poisson_model_", habitat, "_", x, ".png"),
         device = "png")
  
  
  core_otu_table <- left_join(core_bc$BC_ranks, neut_model$occ_abun_core, by = "OTU_ID") %>%
    rename(AverageRelativeAbundance = otu_rel,
           OccupancyProp = otu_occ,
           CoreStatus = fill)
  
  write.table(core_otu_table, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_", x, ".txt"),
              sep = "\t", row.names = F)
  
  
  return(list(BC_increase = core_bc$BC_ranks,
              model_pred_output = neut_model$model_pred_output,
              model_comp_stats = neut_model$model_comp_stats,
              core_otus = core_otu_table))
})

names(Palsa_list) <- paste0(habitat, "_", depth_levels)


PalsaBCCore <- bind_rows(Palsa_list[[1]]$core_otus %>% mutate(Habitat = "Palsa", Depth = "0-9"),
          Palsa_list[[2]]$core_otus %>% mutate(Habitat = "Palsa", Depth = "10-19"),
          Palsa_list[[3]]$core_otus %>% mutate(Habitat = "Palsa", Depth = "20-29"),
          Palsa_list[[4]]$core_otus %>% mutate(Habitat = "Palsa", Depth = "30-39"))

write.table(PalsaBCCore, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_by_depth.txt"),
            sep = "\t", row.names = F)
#### ====================================================================== ####

# For Bog
#### ====================================================================== ####
habitat <- "Bog"
bog_samples <- input_counts$sample_metadata %>% filter(Habitat__ == "Bog") %>% 
  select(temporal_sample_id)
bog_otus <- input_counts$otu_table %>%
  select(genome, any_of(bog_samples$temporal_sample_id))

bog_core_bc <- calculate_bc_core_contribution_time(otu = bog_otus[,-1],
                                                     map = input_counts$sample_metadata,
                                                     n_top_otus = 500, display_plots = T)

bog_core_th <- create_core_thresholds(BC_ranks = bog_core_bc$BC_ranks)
bog_core_th$plot

ggsave(file = paste0(figures.fp, "/core_th_", habitat, ".png"),
       plot = bog_core_th$plot, device = "png")

bog_neut_model <- apply_sloan_neutral_model(otu = bog_otus[,-1],
                                              map = input_counts$sample_metadata,
                                              BC_ranks = bog_core_bc$BC_ranks,
                                              occ_abun = bog_core_bc$occ_abun,
                                              cutoffs = bog_core_th,
                                              fill_option = "Bog")
bog_neut_model$occup_core_plot
bog_neut_model$model_plot
bog_neut_model$binom_model_plot
bog_neut_model$poisson_model_plot

write.table(bog_core_otu_table, file = paste0(outputs.fp, "/core_otus_table_", "bog", ".txt"),
            sep = "\t", row.names = F)

# Save the occupancy core plot
ggsave(plot = bog_neut_model$occup_core_plot,
       file = paste0(figures.fp, "/occup_core_", habitat, ".png"),
       device = "png")
ggsave(plot = bog_neut_model$model_plot,
       file = paste0(figures.fp, "/neutral_model_", habitat, ".png"),
       device = "png")
ggsave(plot = bog_neut_model$binom_model_plot,
       file = paste0(figures.fp, "/bionom_model_", habitat, ".png"),
       device = "png")
ggsave(plot = bog_neut_model$poisson_model_plot,
       file = paste0(figures.fp, "/poisson_model_", habitat, ".png"),
       device = "png")


# Bog by depth:
habitat <- "Bog"
Bog_list <- lapply(depth_levels, function(x) {
  which_samples <- input_counts$sample_metadata %>% filter(Habitat__ == habitat) %>% 
    filter(DepthLumping == x) %>%
    select(temporal_sample_id)
  which_otus <- input_counts$otu_table %>%
    select(genome, any_of(which_samples$temporal_sample_id))
  
  # Step 1: calculate the core BC of each of the top 500 taxa
  core_bc <- calculate_bc_core_contribution_time(otu = which_otus[,-1],
                                                 map = input_counts$sample_metadata,
                                                 n_top_otus = 1662, display_plots = T)
  
  # Step 2: calculate the core BC of each of the top 500 taxa
  core_th <- create_core_thresholds(BC_ranks = core_bc$BC_ranks)
  
  ggsave(file = paste0(figures.fp, "/core_th_", habitat, "_", x, ".png"),
         plot = core_th$plot, device = "png")
  
  # Step 3: Apply the sloan neutral model
  neut_model <- apply_sloan_neutral_model(otu = which_otus[,-1],
                                          map = input_counts$sample_metadata,
                                          BC_ranks = core_bc$BC_ranks,
                                          occ_abun = core_bc$occ_abun,
                                          cutoffs = core_th, 
                                          fill_option = habitat)
  
  # Save the occupancy core plot
  ggsave(plot = neut_model$occup_core_plot,
         file = paste0(figures.fp, "/occup_core_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$model_plot,
         file = paste0(figures.fp, "/neutral_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$binom_model_plot,
         file = paste0(figures.fp, "/bionom_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$poisson_model_plot,
         file = paste0(figures.fp, "/poisson_model_", habitat, "_", x, ".png"),
         device = "png")
  
  
  core_otu_table <- left_join(core_bc$BC_ranks, neut_model$occ_abun_core, by = "OTU_ID") %>%
    rename(AverageRelativeAbundance = otu_rel,
           OccupancyProp = otu_occ,
           CoreStatus = fill)
  
  write.table(core_otu_table, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_", x, ".txt"),
              sep = "\t", row.names = F)
  
  
  return(list(BC_increase = core_bc$BC_ranks,
              model_pred_output = neut_model$model_pred_output,
              model_comp_stats = neut_model$model_comp_stats,
              core_otus = core_otu_table))
})

names(Bog_list) <- paste0(habitat, "_", depth_levels)


BogBCCore <- bind_rows(Bog_list[[1]]$core_otus %>% mutate(Habitat = "Bog", Depth = "0-9"),
                         Bog_list[[2]]$core_otus %>% mutate(Habitat = "Bog", Depth = "10-19"),
                         Bog_list[[3]]$core_otus %>% mutate(Habitat = "Bog", Depth = "20-29"),
                         Bog_list[[4]]$core_otus %>% mutate(Habitat = "Bog", Depth = "30-39"))
write.table(BogBCCore, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_by_depth.txt"),
            sep = "\t", row.names = F)
#### ====================================================================== ####


# For Fen
#### ====================================================================== ####
habitat <- "Fen"
fen_samples <- input_counts$sample_metadata %>% filter(Habitat__ == "Fen") %>% 
  select(temporal_sample_id)
fen_otus <- input_counts$otu_table %>%
  select(genome, any_of(fen_samples$temporal_sample_id))
fen_otus_rarefied <- input_counts$otu_table %>%
  select(genome, any_of(fen_samples$temporal_sample_id)) %>%
  mutate(across(where(is.numeric), ~floor(10*.))) %>% # use floor rather than round so we simply make the detection limit more corse rather than messing with the ; can adjust the detection limit by multiplying column by factor of 10
  select_if(negate(function(col) is.numeric(col) && sum(col) < 2000)) %>% #dim()
  select(-genome) %>% #colSums() %>% sort() %>% hist()
  t() %>%
  rrarefy(., sample = 2000) %>% t() %>%
  as_tibble(rownames = "genome") %>%
  mutate(across(where(is.numeric), ~./10)) %>%
  mutate(genome2 = genome) %>%
  column_to_rownames("genome2")


fen_core_bc <- calculate_bc_core_contribution_time(otu = fen_otus[,-1],
                                                   map = input_counts$sample_metadata,
                                                   n_top_otus = 500, display_plots = T)

fen_core_th <- create_core_thresholds(BC_ranks = fen_core_bc$BC_ranks)
fen_core_th$plot

ggsave(file = paste0(figures.fp, "/core_th_", habitat, ".png"),
       plot = fen_core_th$plot, device = "png")


fen_neut_model <- apply_sloan_neutral_model(otu = fen_otus[,-1],
                                            map = input_counts$sample_metadata,
                                            BC_ranks = fen_core_bc$BC_ranks,
                                            occ_abun = fen_core_bc$occ_abun,
                                            cutoffs = fen_core_th)
fen_neut_model$occup_core_plot
fen_neut_model$model_plot
fen_neut_model$binom_model_plot
fen_neut_model$poisson_model_plot

write.table(fen_core_otu_table, file = paste0(outputs.fp, "/core_otus_table_", "fen", ".txt"),
            sep = "\t", row.names = F)

# Save the occupancy core plot
ggsave(plot = fen_neut_model$occup_core_plot,
       file = paste0(figures.fp, "/occup_core_", habitat, ".png"),
       device = "png")
ggsave(plot = fen_neut_model$model_plot,
       file = paste0(figures.fp, "/neutral_model_", habitat, ".png"),
       device = "png")
ggsave(plot = fen_neut_model$binom_model_plot,
       file = paste0(figures.fp, "/bionom_model_", habitat, ".png"),
       device = "png")
ggsave(plot = fen_neut_model$poisson_model_plot,
       file = paste0(figures.fp, "/poisson_model_", habitat, ".png"),
       device = "png")



# Fen by depth:
habitat <- "Fen"
Fen_list <- lapply(depth_levels, function(x) {
  which_samples <- input_counts$sample_metadata %>% filter(Habitat__ == habitat) %>% 
    filter(DepthLumping == x) %>%
    select(temporal_sample_id)
  which_otus <- input_counts$otu_table %>%
    select(genome, any_of(which_samples$temporal_sample_id))
  
  # Step 1: calculate the core BC of each of the top 500 taxa
  core_bc <- calculate_bc_core_contribution_time(otu = which_otus[,-1],
                                                 map = input_counts$sample_metadata,
                                                 n_top_otus = 1662, display_plots = T)
  
  # Step 2: calculate the core BC of each of the top 500 taxa
  core_th <- create_core_thresholds(BC_ranks = core_bc$BC_ranks)
  
  ggsave(file = paste0(figures.fp, "/core_th_", habitat, "_", x, ".png"),
         plot = core_th$plot, device = "png")
  
  # Step 3: Apply the sloan neutral model
  neut_model <- apply_sloan_neutral_model(otu = which_otus[,-1],
                                          map = input_counts$sample_metadata,
                                          BC_ranks = core_bc$BC_ranks,
                                          occ_abun = core_bc$occ_abun,
                                          cutoffs = core_th, 
                                          fill_option = habitat)
  
  # Save the occupancy core plot
  ggsave(plot = neut_model$occup_core_plot,
         file = paste0(figures.fp, "/occup_core_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$model_plot,
         file = paste0(figures.fp, "/neutral_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$binom_model_plot,
         file = paste0(figures.fp, "/bionom_model_", habitat, "_", x, ".png"),
         device = "png")
  ggsave(plot = neut_model$poisson_model_plot,
         file = paste0(figures.fp, "/poisson_model_", habitat, "_", x, ".png"),
         device = "png")
  
  
  core_otu_table <- left_join(core_bc$BC_ranks, neut_model$occ_abun_core, by = "OTU_ID") %>%
    rename(AverageRelativeAbundance = otu_rel,
           OccupancyProp = otu_occ,
           CoreStatus = fill)
  
  write.table(core_otu_table, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_", x, ".txt"),
              sep = "\t", row.names = F)
  
  
  return(list(BC_increase = core_bc$BC_ranks,
              model_pred_output = neut_model$model_pred_output,
              model_comp_stats = neut_model$model_comp_stats,
              core_otus = core_otu_table))
})

names(Fen_list) <- paste0(habitat, "_", depth_levels)


FenBCCore <- bind_rows(Fen_list[[1]]$core_otus %>% mutate(Habitat = "Fen", Depth = "0-9"),
                       Fen_list[[2]]$core_otus %>% mutate(Habitat = "Fen", Depth = "10-19"),
                       Fen_list[[3]]$core_otus %>% mutate(Habitat = "Fen", Depth = "20-29"),
                       Fen_list[[4]]$core_otus %>% mutate(Habitat = "Fen", Depth = "30-39"))

write.table(FenBCCore, file = paste0(outputs.fp, "/core_otu_table_", habitat, "_by_depth.txt"),
            sep = "\t", row.names = F)
#### ====================================================================== ####