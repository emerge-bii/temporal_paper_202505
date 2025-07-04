#' ## Calculate RCBC
#' This step calculates the raup-crick bray-curtis turnover across samples. 
#' This step takes a long time and would be much better run on a super computer
#' Code adapted from Danczak et al. 2020; https://journals.asm.org/doi/full/10.1128/mSystems.00098-20


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
library(data.table); # for loading in large files
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

#### =================================================================================================== ####

#' First we'll set up user settings and command line arguments: 
#' Users should modify these settings for their own environment.
#'  - Set the number of replicate randomizations for rcbc
#'  - set the relative abundance flag
#'  - if you're working in an environment with parallel computing set `paral <- TRUE`
#'  - if speedup <- TRUE only non-significant betaNTI calculations will be run in this step
#'  - use the max_chunk_size to control how much memory is allocated in total
#'  - if you have run this program once and already have the rc-bray-curtis null, you can specify it here
#### =================================================================================================== ####
beta.reps <- 10000; # number of randomizations for rcbc

rel_abund <- TRUE # is the otu table in relative abundances such as read mappings?

speedup <- FALSE # speed up the process by filtering the community comparisions for only those that aren't considered under selection from betaNTI calculations; currently this option is not implemented in the code below

paral <- TRUE

# Set this if you're having memory blow-out problems
max_chunk_size <- 1000000

# Should the program skip the generation of null bray curtis, and instead read it in from a file?
skip_null <- FALSE

# Tell the user their settings:
writeLines(paste("User settings are as follows: \n",
                 "Parallel:", paral, "\n",
                 "Speedup:", speedup, "\n",
                 "Replicates:", beta.reps, "\n",
                 "Relative Abundances:", rel_abund, "\n",
                 "Maximum chunk size:", max_chunk_size, "\n",
                 "Skip Null Bray Curtis:", skip_null))

if(paral==T) {
  library(Rmpi)  # R implementation of MPI interface
  library(doMPI)
  library(iterators)
  n.cores <- max(1, mpi.universe.size()-1) # set n.cores to 1 or mpi allocation size, whichever is bigger
  print(n.cores)
  print(paste(n.cores, "worker nodes are available."))
  
  if (n.cores > 24) {
    # we are using more than 1 CPU, so really run in parallel mode
    cl <- startMPIcluster(count = n.cores)
    registerDoMPI(cl) # tell foreach about the cluster
    print(paste("Running in parallel mode on",n.cores,"worker nodes."))
    
  } else if (n.cores > 1) {
    library(doParallel)
    library(iterators)
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

# Set the number of cores available for data.table to use
setDTthreads(threads = n.cores)

#### =================================================================================================== ####

#### Read in data 
#### =================================================================================================== ####
# Read in OTU table and env data
input <- input_counts

# Read in tree
magtree <- input_counts$tree

# Read in betaNTI results, if we want to speedup based on this step
if(speedup == TRUE) {
  betaNTI <- read_csv(file = paste0(outputs.fp, "/weighted_bNTI.csv")) %>%
    column_to_rownames(var = "X1")
}


#### =================================================================================================== ####

#' ## Step 1: Prepare data to calculate null RCBC communities
#### =================================================================================================== ####
# Setting up the data to mesh with the Stegen et al. code
#spXsite = t(input$otu_table[,2:10]) # test with smaller sample set
#dim(spXsite)
spXsite = t(input$otu_table[-1]) 
dim(spXsite)

# If speedup is true, we will only do comparisons for the comparisons that were
# identified as stochastic in part one with BetaNTI
if(speedup == T) {
  # Prepare list of stochastic sample comparisons to filter sample list by
  betaNTI_stoc <- betaNTI #filter to match spXsite
  betaNTI_stoc[betaNTI > 2 | betaNTI < -2 ] <- NA
  # Convert row/col names to indices
  betaNTI_stoc <- betaNTI_stoc[rownames(spXsite), 
                               rownames(spXsite)] #filter to match spXsite
  colnames(betaNTI_stoc) <- 1:ncol(betaNTI_stoc) 
  rownames(betaNTI_stoc) <- 1:nrow(betaNTI_stoc)
}

# if relative abundance is true add 1 to all values greater than zero otherwise 
# raup crick gets cranky
if(rel_abund == TRUE){
  spXsite[spXsite>0] = spXsite[spXsite>0]+1
} # Adding one to my RPKM value to generate an "RPKM+1" stat. This is necessary to work around Raup-Crick

# Count number of sites and total species richness across all plots (gamma diversity)
n_sites = nrow(spXsite)
gamma = ncol(spXsite)

# Build a site-by-site matrix to hold the results, with the names of the sites in the row and col names:
results = matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

# Make the spXsite matrix into a new, pres/abs. matrix to get species occurences:
spXsite.inc = ceiling(spXsite/max(spXsite))

# Create an occurrence vector - used to give more weight to widely distributed species in the null model
occur = apply(spXsite.inc, MARGIN=2, FUN=sum)

# Create an abundance vector - used to give more weight to abundant species in the second step of the null model
abundance = apply(spXsite, MARGIN=2, FUN=sum)
#### =================================================================================================== ####

#' ## Step 2: Create null pairwise comparisons (beta.reps number of times) for each sample
#### =================================================================================================== ####
# Tell the user the time/date and what we are doing
print(date())
if(!skip_null){
   writeLines("Calculating null comparisons")
}

# Function definitions for step 2
# First set up functions to create independent site_comb combinations so each can be run
# in parallel on a super computer; the calculation for each site-site combination is independent so this is fine
################################################################################
# create_spec_comb creates a list of pairwise comparisons for every sample and
# repeats that beta.reps number of times
create_spec_comb <- function(nsites = nrow(spXsite), nreps = beta.reps) {
  # This function creates a matrix of site combinations that
  # describe the lower-half of a sitexsite matrix, while excluding the 
  # self-comparison diagonal. In addition, it repeats each sitexsite comparison
  # nreps number of times. The list can then be used in the create_null function
  # below to parallelize (via apply) each site-by-site null comparison.
  z <- sequence(nsites)

  mat_comb <- cbind(
    null.two = unlist(lapply(2:nsites, function(x) x:nsites), use.names = FALSE),
    null.one = rep(z[-length(z)], times = rev(tail(z, -1))-1)) # a two-site matrix for our null comparisons
  
  if(speedup == TRUE) {
    stoc_samples <- unlist(betaNTI_stoc[lower.tri(betaNTI_stoc)]) # vector of samples where non-stochastic are NA and stochastic are betaNTI score
    mat_comb <- mat_comb[!is.na(stoc_samples),] # only keep samples with betaNTI between -2 and 2
  }
  
  mat_rep <- matrix(unlist(lapply(mat_comb, rep, nreps)), ncol = 2)
  nrow_mat_comb <- nrow(mat_comb)
  # Replicate the site combination nreps number of times so that sites combinations
  # are in order
  mat_rep_iter <- cbind(mat_rep[,c(1,2)], iter=rep(1:nreps, times = nrow_mat_comb))
  
  return(mat_rep_iter)
}

# create null creates two null communities and calculates bray-curtis distance between them for a site1-site2-iteration combination
create_null <- function(site_comb, spXsite = spXsite) {
  # This function creates the null comparisons between sites
  # by first comparing the observed number of species in each of the two sites and
  # and weighing them by their occurrence frequencies, but otherwise randomizing;
  # the species that are in the community. Then the two sites are combined into a
  # matrix, and bray-curtis distance is measured for the two-site comparision. The
  # results are saved along with the sites being compared and the beta.rep iteration.
  #create_null <- function(site_comb, spXsite = spXsite) {
  # This function creates the null comparisons between sites
  # by first comparing the observed number of species in each of the two sites and
  # and weighing them by their occurrence frequencies, but otherwise randomizing;
  # the species that are in the community. Then the two sites are combined into a
  # matrix, and bray-curtis distance is measured for the two-site comparision. The
  # results are saved along with the sites being compared and the beta.rep iteration.
  #
  # ### Function Arguments: ###
  # site_comb: is a list of site combinations and iterations, created by
  #            create_spec_comb() function, it takes the format of a list of 
  #            vectors each vector has a length of three and the first item is 
  #            "null.two" (ranges from 2:number_of_sites), the second item is 
  #            "null.one" (ranges from 1:number_of_sites-1), the third is the 
  #            replicate for the null distribution (ranges from 
  #            1:number_of_beta_reps). Together null.one and null.two describe 
  #            bottom triangle of a site-by-site matrix (excluding 
  #            self-comparisons)
  # Here's an example of what site_comb could look like:
  # 
  # List of 2340
  # $ X1   : num [1:3] 2 1 1
  # $ X2   : num [1:3] 2 1 2
  # $ X3   : num [1:3] 2 1 3
  # $ X4   : num [1:3] 2 1 4
  # etc....
  # 
  # spXsite: is the species by site matrix
  null.two <- unlist(site_comb)[1]
  null.one <- unlist(site_comb)[2]
  iter <- unlist(site_comb)[3]
  #print(paste(iter, null.one, null.two))
  #writeLines("Data read in")
  
  # Generates two empty communities of size gamma (# species in regional pool)
  com1<-rep(0,gamma)
  com2<-rep(0,gamma)
  #writeLines("Empty communities generated")
  # Add observed number of species to com1, weighting by species occurrence frequencies
  # In first sample, choose a random subset of species to form the community: x = choose from any of the gamma # of species in the region, size = choose a vector of species with a size equal to the number of species in original community; Weigh the probability of choosing a particular species by its regional occurrence frequency ("occur"); Then assign this random subset of regional species a presence value of "1"
  com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
  
  # using the pool of randomly chosen species in com1, randomly sample x species, where x is difference in # of total individuals (or abundance counts) in the real community one and the total number of observations (pres/abs) in the real community (which is equivalent to sum(com1), and sum(spXite.inc[null.one,])).
  # The chance of choosing a species is propotionate to their abundances in the real regional community,however species may be chosen with replacement, so it is possible to pick the same species twice or more
  com1.samp.sp <- sample(which(com1>0), (sum(spXsite[null.one,])-sum(com1)), 
                         replace=TRUE, prob=abundance[which(com1>0)]); # random draws of species from the com1 community; chosen propotionate to their regional abundances. The total size of the community is bounded by the size of the real community membership
  com1.samp.sp <- cbind(com1.samp.sp,1); # hist(com1.samp.sp[,1]) = randomly generated species distribution for com1
  com1.samp.sp
  # count the number of occurences for each species
  com1.sp.counts <- as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
  com1.sp.counts$sp <- as.numeric(rownames(com1.sp.counts));
  com1[com1.sp.counts$sp] <- com1[com1.sp.counts$sp] + com1.sp.counts$counts;
  #sum(com1) - sum(spXsite[null.one,]); # This should be zero if everything worked properly
  rm('com1.samp.sp','com1.sp.counts');			
  #writeLines("Community 1 calculated")
  
  # Again for com2
  com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
  com2.samp.sp <- sample(which(com2>0), (sum(spXsite[null.two,])-sum(com2)), replace=TRUE, prob=abundance[which(com2>0)]);
  com2.samp.sp <- cbind(com2.samp.sp,1);
  com2.sp.counts <- as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
  com2.sp.counts$sp <- as.numeric(rownames(com2.sp.counts));
  com2[com2.sp.counts$sp] <- com2[com2.sp.counts$sp] + com2.sp.counts$counts;
  # sum(com2) - sum(spXsite[null.two,]); # This should be zero if everything worked properly
  rm('com2.samp.sp','com2.sp.counts');
  #writeLines("Community 2 calculated")
  
  # Bind the two null communities together to create 2-samplexspecies table for the null
  # communities
  null.spXsite <- rbind(com1,com2); # Null.spXsite
  #writeLines("null df created") # check to report to the user - won't work when doMPI is running
  
  # Calculates the null Bray-Curtis
  temp_null_bray_curtis <- as.numeric(vegdist(null.spXsite, method='bray'));
  #writeLines("ran vegan") # check to report to the user - won't work when doMPI is running
  
  # Final report out:
  null_results <- list("null_bc" = temp_null_bray_curtis, 
                       "site1" = null.one,
                       "site2" = null.two,
                       "null_iter" = iter)
  return(null_results)
}

################################################################################
# Now run those functions: 
if(!skip_null) {
  # Create list of site-by-site comparisons (quick doesn't need to be run in parallel)
  site_comb_list <- create_spec_comb(nsites = nrow(spXsite), nreps = beta.reps)

  # Tell the user what's going down
  writeLines(paste("We will be running a total of", nrow(site_comb_list), "comparisons."))
  writeLines(paste("These are the product of", beta.reps, "replications of", nrow(spXsite),   "pairwise sample comparisons with self-comparisons removed. If speedup is set to TRUE, non   -stochastic comparisions have also been removed reducing the total number further."))
}

# Create an iterator for the foreach loop, that speeds up parallelization
# This iterator chops up the comb_list into pieces to hand off to each worker
# it also saves the iterator i and m for use within the foreach function
# Currently the non-parallel version of this is not tested, but probably will
# not work since I haven't updated the iterator.
if(paral == TRUE) {
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
  } else {
  iter_seq <- 1:length(site_comb_list)
  }

# Set the chunk size
chnk_size <- min(floor(nrow(site_comb_list)/n.cores), max_chunk_size)
writeLines(paste("To speed up processing, chunk size is set to", chnk_size, "which means there will be", nrow(site_comb_list)/chnk_size, "tasks to complete."))

# print(nrow(site_comb_list)) # for troubleshooting
# str(site_comb_list) # for troubleshooting

if(!skip_null) {
  # Set up the tmp directory to report out - this is important otherwise
  # the final combination object might be too big to hold in memory. 
  # Writing the outputs to temporary files and then reading them in solves that problem
  tmp_null <- here(outputs.fp, "tmp")
  if(!dir.exists(tmp_null)) {
     dir.create(tmp_null)
  } else {
     # Remove any files already in the tmp directory
     f <- list.files(tmp_null, include.dirs = F, full.names = T, recursive = T)
     file.remove(f)
  }
  
  null_bray_curtis.ls <- foreach(a=idivix(nrow(site_comb_list), 
                                          chunkSize=chnk_size,
                                          comb_list=site_comb_list),
                                 .packages = c("vegan"),
                                 .init = NULL,
                                 .combine = rbind,
                                 .multicombine = TRUE,
                                 .maxcombine = 100,
                                 .inorder = FALSE,
                                 .errorhandling = "remove",
                                 .verbose = TRUE) %dopar% {
      null.tmp <- do.call('rbind', lapply(seq(1, length.out=a$m), function(i) {
        site_comb_temp <- as.list(a$comb[i,])
        create_null(site_comb = site_comb_temp, spXsite = spXsite)}))
      # Because the files can get large, write the null comparisons to disk
      result_file_name <- paste0(a$comb[1,], collapse = "_")
      tmp <- matrix(unlist(null.tmp), ncol = 4)
      write.csv(tmp, file = paste0(tmp_null,"/", result_file_name, ".csv"),
                row.names = FALSE)                           
  }
  writeLines("Finished calculating null replicates at:")
  print(date())
  
  # Read in temporary files
  null_bray_curtis <- list.files(path = tmp_null, pattern = "*.csv", full.name = T) %>% 
                        map_df(~fread(.))
  writeLines("Null_bray_curtis")
  # head(null_bray_curtis) # Check
  # tail(null_bray_curtis) # Check
  # str(null_bray_curtis) # Check
  
  ## Prepare outputs to pass along to next step
  colnames(null_bray_curtis) <- c("null_bc", "site1", "site2", "null_iter")
  null_bray_curtis <- as.data.table(null_bray_curtis)
  # Sort null bray curtis, in case it's not already sorted
  null_bray_curtis <- null_bray_curtis[order(site1, site2, null_iter, decreasing = TRUE)]
  
  # head(null_bray_curtis) # debugging check, should be a dataframe with 4 columns
  # str(null_bray_curtis) # debugging check
  # null_bray_curtis[ which(null_bray_curtis$site1==1 & null_bray_curtis$site2==2), "null_bc"] # debugging check; should be length of beta.reps and show the null comparisions for sites 1 and 2

  gc() # memory control - unclear if necessary but gives date and time, so it's nice

  ## Save objects mid-process incase we have to restart.

  saveRDS(null_bray_curtis,  paste0(outputs.fp, "/null_bray_curtis.RDS"))
  write.csv(null_bray_curtis,
            paste0(outputs.fp, "/null_bray_curtis_lf.csv"), row.names = FALSE,
            quote=FALSE)
}
gc()
#### =================================================================================================== ####


#' ## Step 3: Compare null comparisons to the observed comparisons
#### =================================================================================================== ####
# Tell the user the date and what we are doing
writeLines("comparing null to obs")
print(date())

# If skipping null generation:
if(skip_null) {
  writeLines("Reading in null bray curtis from:")
  print(paste0(outputs.fp, "/null_bray_curtis_lf.csv"))
  writeLines(paste0("null_bray_curtis_lf.csv was created: ", 
             file.info(paste0(outputs.fp, "/null_bray_curtis_lf.csv"))$ctime))
  
  null_bray_curtis <- fread(paste0(outputs.fp, "/null_bray_curtis_lf.csv"))
  null_bray_curtis <- null_bray_curtis[, !c("V1")] # get rid of rownames column
  # Debugging
  # str(null_bray_curtis)
  # head(null_bray_curtis)
  null_bray_curtis <- null_bray_curtis[order(site1, site2, null_iter, decreasing = TRUE)]
}


# Function definition
################################################################################
compare_null <- function(site_comb, null_bray_curtis, beta.reps = beta.reps, spXsite = spXsite) {
  # This function compares the null comparisons with the observations
  #
  # ### Function Arguments: ###
  # site_comb: is a list of site combinations and iterations, created by
  #            create_spec_comb() function, it takes the format of a list of 
  #            vectors each vector has a length of three and the first item is 
  #            "null.two" (ranges from 2:number_of_sites), the second item is 
  #            "null.one" (ranges from 1:number_of_sites-1), the third is the 
  #            replicate for the null distribution (ranges from 
  #            1:number_of_beta_reps). Together null.one and null.two describe 
  #            bottom triangle of a site-by-site matrix (excluding 
  #            self-comparisons)
  # Here's an example of what site_comb could look like:
  # 
  # List of 2340
  # $ X1   : num [1:3] 2 1 1
  # $ X2   : num [1:3] 2 1 2
  # $ X3   : num [1:3] 2 1 3
  # $ X4   : num [1:3] 2 1 4
  # etc....
  #
  # null_bray_curtis = the null_bray_curtis data table
  # beta.reps = number of replicates for each combination.
  #
  # spXsite: is the species by site matrix
  null.one <- unlist(site_comb)[1]
  null.two <- unlist(site_comb)[2]
  
  # Calculates the observed Bray-Curtis
  obs.bray = vegdist(spXsite[c(null.one,null.two),], method='bray');
  
  # Extract the null bray curtis vector for the site under calculation
  null.bc.site <- null_bray_curtis[site1==null.one & site2==null.two, "null_bc"]

  # How many null observations is the observed value tied with?
  num_exact_matching_in_null = sum(null.bc.site==obs.bray);
  
  # How many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null = sum(null.bc.site<obs.bray);
  
  rc = ((num_less_than_in_null + (num_exact_matching_in_null)/2)/beta.reps) # This variation of rc splits ties
  
  rc = (rc-.5)*2 # Adjusts the range of the  Raup-Crick caclulation to -1 to 1
  
  # Final report out:
  results.lf <- list("RCBC" = rc, 
                     "site1" = null.one,
                     "site2" = null.two)
  return(results.lf)
}
################################################################################
# Create list of site-by-site comparisons
# nreps should = 1 since we only need to do a null-obs comparision for each pair of sites one time
rm(site_comb_list)
site_comb_list <- create_spec_comb(nsites = nrow(spXsite), nreps = 1)

# Use foreach to compare null to obs in parallel
chnk_size <- min(floor(nrow(site_comb_list)/n.cores), max_chunk_size)
chnk_size <- max(n.cores, chnk_size) # n.cores represents the low end of tasks we want to run
writeLines(paste("To speed up processing, chunk size is set to", chnk_size, "which means there will be", nrow(site_comb_list)/chnk_size, "tasks to complete."))

if(paral == TRUE) {
# Define the iterator interpretor function
# this one will return slices of a list no larger than chunkSize
# Importantly, it needs to subset null_bray_curtis to a small enough
# size to operate on, BUT all iterations for each null site combination
# must be included so the length.out can't be smaller than beta.reps.
idivix <- function(n, chunkSize, null_bray_curtis, beta.reps) {
  i <- 1
  it <- idiv(n, chunkSize=chunkSize)
  
  nextEl <- function() {
    m <- nextElem(it)  # may throw 'StopIterator'
    # make sure that each iteration starts on beta.reps + 1
    beta_start <- ifelse(i != 1, ((i-1)*beta.reps) + 1, i)
    nullbc <- null_bray_curtis[seq(beta_start, length.out=m*beta.reps),]
    comb <- unique(nullbc[,c("site1", "site2")])
    # Create a list of iterator outputs
    # i = start of the iterator site combination
    # m = length of the number of site combinations in this iterator
    # comb = a data.table with the unique combinations of sites (no replicates)
    # nullbc = the sliced null_bray_curtis data.table for this subset of sites
    # the number of rows in comb (aka the number of site combinations handled at one time)
    value <- list(i=i, m=m, comb = comb, nullbc = nullbc, ncomb = nrow(comb))
    i <<- i + m
    value
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}
}

# Set up the tmp directory to report out - this is important otherwise
# the final combination object might be too big to hold in memory. 
# Writing the outputs to temporary files and then reading them in solves that problem
tmp_comp <- here(outputs.fp, "tmp_2")
  if(!dir.exists(tmp_comp)) {
     dir.create(tmp_comp)
  } else {
     # Remove any files already in the tmp directory
     f <- list.files(tmp_comp, include.dirs = F, full.names = T, recursive = T)
     file.remove(f)
  }

writeLines("Starting results loop")
results.lf.ls <- foreach(a=idivix(nrow(site_comb_list), 
                                  chunkSize=chnk_size,
                                  null_bray_curtis=null_bray_curtis,
                                  beta.reps = beta.reps),
                         .packages = c("vegan", "data.table"),
 	       	       	       .init = NULL,
                         .combine = rbind,
                         .multicombine = TRUE,
                         .maxcombine = 100,
                         .inorder = FALSE,
                         .errorhandling = "remove",
                         .verbose = TRUE) %dopar% {
  comp.tmp <- do.call('rbind', lapply(seq(1, length.out=a$ncomb), function(i) {
    site_comb_temp <- as.list(a$comb[i,])
    compare_null(site_comb = site_comb_temp, beta.reps = beta.reps, null_bray_curtis = a$nullbc, spXsite = spXsite)
    })) 
  # Because memory can be large write to disk
  result_file_name <- paste0(a$comb[1,], collapse = "_")
  tmp <- matrix(unlist(comp.tmp), ncol = 3)
  write.csv(tmp, file = paste0(tmp_comp,"/", result_file_name, ".csv"),
                row.names = FALSE)
}

# Read in temporary files
writeLines("reading in tmp files")
results.lf <- list.files(path = tmp_comp, pattern = "*.csv", full.name = T) %>%
                        map_df(~fread(.))
colnames(results.lf) <- c("RCBC", "site1", "site2")

# Convert results data.table (in long format) to distance matrix (in wide format)
# for convenience. We will save and write out both of these products
self_comp <- data.frame(site1 = rep(1:nrow(spXsite)),
                        site2 = rep(1:nrow(spXsite)),
                        RCBC = NA) # self-comparisons

writeLines("Reformatting results")
# combine self comparisons to results and spread to distance matrix format
results.wf <- results.lf %>%
  select(RCBC, site1, site2) %>%
  bind_rows(self_comp) %>%
  arrange(site1, site2) %>%
  mutate(site1 = as.character(site1),
         site2 = as.character(site2)) %>%
  pivot_wider(names_from = c(site1),
              values_from = RCBC) %>%
  column_to_rownames("site2")

#### =================================================================================================== ####

#### Save Data and Figures
#### =================================================================================================== ####
# 
# Data
# RCBC matrix
saveRDS(results.wf, paste0(outputs.fp, "/rcbc_matrix.RDS"))
write.csv(results.wf,
          paste0(outputs.fp, "/rcbc_matrix.csv"),
          quote=FALSE)

# Long formats outputs from foreach loops
write.csv(results.lf,
          paste0(outputs.fp, "/rcbc_long_form.csv"),
          quote=FALSE, row.names = FALSE)


#### Cleanup parallel
#### =================================================================================================== ####

if(paral == TRUE) {
  if (n.cores > 24) {
     writeLines("closing mpi cluster")
     # Note: due to this problem: https://stackoverflow.com/questions/41007564/stopcluster-in-r-snow-freeze
     # using closeCluster will cause the job to hang. To bruteforce it, I will use mpi.exit/mpi.quit instead
     #doMPI::closeCluster(cl) # due to this problem: 

     writeLines("Stopping mpi")
     Rmpi::mpi.exit()
     Rmpi::mpi.quit()  # or mpi.quit(), which quits R as well

  } else if (n.cores > 1) {
    writeLines("closing parallel fork cluster")
    parallel::stopCluster(cl)
    Rmpi::mpi.quit()
  }
}
