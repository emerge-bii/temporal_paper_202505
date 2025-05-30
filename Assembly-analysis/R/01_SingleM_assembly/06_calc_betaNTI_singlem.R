#' ## Calculate betaNTI
#' This step calcluates the betaNTI or phylogenetic turnover among samples
#' It requires a phylogenetic tree and a species-by-sample table
#' This step takes a long time and would be much better run on a super computer
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
library(picante); packageVersion("picante") # handles tree manipulations
library(vegan); packageVersion("vegan") # for ecological applications
library(ape)
library(here)
library(doMPI)

# Load required data
source(here("Assembly-analysis", "R", "05_prepare_singlem_input.R")) # this will also source the setup script

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("Assembly-analysis", "outputs")
figures.fp <- here("Assembly-analysis", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### =================================================================================================== ####
#' First we'll set up user settings and command line arguments: 
#' Users should modify these settings for their own environment.  
#'  - Set the number of replicate randomizations for the BetaMNTD calculations `beta.reps`
#'  - Set the relative abundance flag if your OTU table is not integers but real numbers instead
#'  - If you're working in an environment with parallel computing set `paral <- TRUE`
# User settings and command line arguments
#### =================================================================================================== ####
beta.reps <- 10000; # number of randomizations for beta MNTD

rel_abund <- TRUE # is the otu table in relative abundances such as read mappings?

paral <- TRUE

# Tell the user their settings:
writeLines(paste("User settings are as follows: \n",
                 "Parallel:", paral, "\n",
                 "Replicates:", beta.reps, "\n",
                 "Relative Abundances:", rel_abund))

# If parallel is true, detect available cores
if(paral==T) {
  library(Rmpi)  # R implementation of MPI interface
  library(doMPI)
  n.cores <- max(1, mpi.universe.size()-1) # set n.cores to 1 or mpi allocation size, whichever is bigger
  #print(n.cores)
  print(paste(n.cores, "worker nodes are available."))
  
  if (n.cores > 24) {
    # we are using more than 1 node, so really run in parallel mode
    cl <- startMPIcluster(count = n.cores)
    registerDoMPI(cl) # tell foreach about the cluster
    print(paste("Running in parallel mode on",n.cores,"worker nodes."))
    
  } else if (n.cores > 1) {
    library(doParallel)
    cl <- parallel::makeForkCluster(n.cores)
    registerDoParallel(cl)
    print(paste("Running in parallel mode without mpi on", n.cores, "worker nodes."))
  } else {
    registerDoSEQ() # tells foreach to use sequential mode
    paral <- FALSE
    print(paste("Only", n.cores,
                "available. Not enough to run in parallel, running in sequential mode."))
  }
} else {
  registerDoSEQ()
  print("Running in sequential mode.")
} 

# Command Line arguments - Specify the singlem Gene - either S3.11 or S3.5 currently:
args = commandArgs(trailingOnly=TRUE)
writeLines("Command Args are...")
print(args)

# test if there is at least one argument: if not, return an error
if (length(args)==1) { # args will always be >= 1 b/c of the --show-progress flag in mpiexec
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  # Check if the specified singlem input variable exists
  input_var <- paste0(args[1], "_input_ra")
  exists(input_var)
  writeLines(paste0("Calculating BetaNTI for ", args[1]))
  
  # Assign input to be the user-specified variable
  input <- eval(as.symbol(input_var)) 
  # Read in tree
  magtree <- input$tree
}

lowmem_cophenetic <- function(tree, chnksize = 1000) {
  # Re-writing cophenetic function based on
  # Liam Revell's solution: http://blog.phytools.org/2015/10/new-reasonably-fast-method-to-compute.html
  # Make combinations of tip labels
  tip_lab_combn <- combn(tree$tip.label, m = 2)
  tip_lab_combn <- t(tip_lab_combn)
  
  # Define the iterator interpretor function
  # this one will return slices of a list no larger than chunkSize
  idivix <- function(n, chunkSize, comb_list) {
    i <- 1
    it <- idiv(n, chunkSize=chunkSize)
    nextEl <- function() {
      m <- nextElem(it)  # may throw 'StopIterator'
      
      value <- list(i=i, m=m, comb = comb_list[seq(i, length.out=m),])
      i <<- i + m
      value
    }
    obj <- list(nextElem=nextEl)
    class(obj) <- c('abstractiter', 'iter')
    obj
  }
  dist_calc <- foreach(a=idivix(nrow(tip_lab_combn), 
                                          chunkSize=chnksize,
                                          comb_list=tip_lab_combn),
                                 .packages = c("phytools"),
                                 .init = NULL,
                                 .combine = rbind,
                                 .multicombine = TRUE,
                                 .maxcombine = 100,
                                 .inorder = FALSE,
                                 .errorhandling = "remove",
                                 .export = c("tip_lab_combn", "tree"),
                                 .verbose = FALSE) %dopar% {
                do.call('rbind', lapply(seq(1, length.out=a$m), function(i) {
                   x <- a$comb[i,]
                   phytools::fastDist(tree, x[1], x[2])}))
                                 }
    
   #dist_calc <- apply(tip_lab_combn, MARGIN = 1, tr = tree,
    #                 FUN = function(x, tr) {
    #tree$tip.label[1]
    #dist <- phytools::fastDist(tr, x[1], x[2])
                       
  #})
  #print(head(dist_calc))
  dist_calc <- data.frame(tip_lab_combn, dist = dist_calc)
  self_comp <- data.frame(X1 = tree$tip.label,
                          X2 = tree$tip.label, dist = 0)
  #print(head(self_comp))
  dist_calc <- dist_calc %>%
    bind_rows(self_comp) %>%
    pivot_wider(names_from = "X1", values_from = "dist") %>% 
    column_to_rownames(var = "X2")
  dist_calc <- dist_calc[tree$tip.label, tree$tip.label]
  dist_calc <- as.matrix(as.dist(dist_calc))
  return(dist_calc)
}

# testtree <- rtree(n = 40)
# testtree <- rtree(n = 100)
# testtree <- rtree(n = 1000)
# testtree <- rtree(n = 10000)
# fc <- lowmem_cophenetic(tree = testtree)
# fc <- lowmem_cophenetic_2(tree = testtree)

#### =================================================================================================== ####


#' ## Read in data 
#### =================================================================================================== ####
# Read in OTU table and env data
writeLines("reading in input data")

# Verify the everything is matching
match.phylo.otu <- match.phylo.data(magtree, input$otu_table[-1])

#### =================================================================================================== ####

#' ## Step 1: First we calculate the empirical betaMNTD
#' This is the observed betaMNTD which we will compare to our null distribution that we will generate below
#### =================================================================================================== ####

# Get the patristic distance for the tree:
#phy_pat_dist <- lowmem_cophenetic(match.phylo.otu$phy)

# phy_pat_dist <- castor::get_all_pairwise_distances(match.phylo.otu$phy,
#                                                    only_clades = c(1:Ntip(magtree)))
# 
# 
# 
# data("example.data")
# 
# phy_pat_dist <- iCAMP::pdist.big(example.data$tree, 
#                                  wd = "~/Downloads", 
#                                  nworker = 2, memory.G = 8, output = TRUE)

tmp_pat <- here(outputs.fp, "tmp_pat")
if(!dir.exists(tmp_pat)) {
  dir.create(tmp_pat)
} else { # if file exits remove contents before running iCAMP
    file.remove(dir(tmp_pat, full.names = TRUE))
  }

phy_pat_dist <- iCAMP::pdist.big(match.phylo.otu$phy, 
                                 wd = tmp_pat, 
                                 nworker = n.cores, memory.G = 128, # amount of memory per node on premise
                                 output = TRUE)

head(phy_pat_dist)

# Calculate non-null expectation

beta.mntd.weighted <- as.matrix(comdistnt(t(match.phylo.otu$data),
                                          phy_pat_dist,
                                          abundance.weighted=T))
dim(beta.mntd.weighted)

#### =================================================================================================== ####

#' ## Step 2: Calculate null expectation (randomized) betaMNTD
#' This is our null distribution. Calculating this will take a while - best to run on server
#### =================================================================================================== ####
null_bMNTD <- NULL # make a null object to hold list of results
null_bMNTD <- foreach (i=1:beta.reps, .packages = c("picante")) %dopar% {
  as.matrix(comdistnt(t(match.phylo.otu$data),
                      taxaShuffle(phy_pat_dist),
                      abundance.weighted=T,
                      exclude.conspecifics = F))
  } # end of parallel computation

# reformat data output from foreach
# Unlisting the list of null beta.reps; This creates an 3d array
# of size: # samples x # samples x # beta.reps
rand.weighted.bMNTD.comp <- array(as.numeric(unlist(null_bMNTD)), 
                                  dim=c(ncol(match.phylo.otu$data),
                                        ncol(match.phylo.otu$data), 
                                        beta.reps))

dim(rand.weighted.bMNTD.comp) # should be # samples x # samples x # beta.reps

#### =================================================================================================== ####

#' ## Step 3: Now calculate betaNTI by normalizing observed against the null expectation
#' Reminder: 
#' -2 < bNTI < 2 indicates stochastic process which can be further refined using values from RCbray in next step.
#' bNTI < -2 indicates homogeneous selection
#' bNTI > 2 indicates heterogeneous selection
#### =================================================================================================== ####
# Save time by setting up matrix first
weighted.bNTI <- matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));

dim(weighted.bNTI); # should be # samples x # samples

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) { # columns + 1 so samples are not compared to themselves
    
    # Pull out the randomized betaMNTD values
    rand.vals <- rand.weighted.bMNTD.comp[rows,columns,]; # length of this object = beta.reps
    # Normalize the observed against the randomized values
    weighted.bNTI[rows,columns] <- (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals); # normalize
    rm("rand.vals");
    
  };
};

# Rename columns and rows to get a samplexsample table
rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
# weighted.bNTI;
#### =================================================================================================== ####

#' Step 4: Verification - Plot Null and empirical distributions
#### =================================================================================================== ####
dist_plot.df <- map_df(null_bMNTD, ~{as.data.frame(.x)}, .id = "null_id") %>% 
  mutate(null_id = paste0("null_", null_id)) %>% 
  bind_rows(.,beta.mntd.weighted %>%
              as_tibble() %>%
              mutate(null_id = "empirical")) %>%
  #group_by(null_id) %>%
  #nest()
  mutate(null_id = as.factor(null_id)) %>%
  pivot_longer(cols = !matches("null_id"),
               names_to = "rep", values_to = "betaMNTD")

dist_plot <- dist_plot.df %>%
  filter(null_id != "empirical") %>%
  ggplot(aes(x = betaMNTD, fill = null_id)) +
  geom_density(alpha = 0.3) +
  scale_fill_grey() +
  geom_density(data = dist_plot.df %>%
                 filter(null_id == "empirical"),
               aes(x = betaMNTD),
               fill = "red", alpha = 0.3) +
  theme(legend.position = "none")

#### =================================================================================================== ####

#### Save Data and Figures
#### =================================================================================================== ####

# Data
# BetaNTI Matrix
write.csv(weighted.bNTI,
          paste0(outputs.fp, "/", args[1], "_weighted_bNTI.csv"),
          quote=FALSE)

# Figures
ggsave(dist_plot, device = "png",
       filename = paste0(figures.fp,  "/", args[1], "_weighted_bNTI_dist.png"))
