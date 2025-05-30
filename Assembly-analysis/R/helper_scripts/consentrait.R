library("data.table")
library("adephylo")
library("ape")

args <- commandArgs(TRUE)


calc_consentrait <- function(tree, traittable, out_dir = get_wd()) {
  # Modified from the Martiny et al. script
  # Modifications include: 
  # - allowing user to set output directory
  # - sensing if a tree is rooted or not and only rooting if not already rooted
  # - sensing if a tree is a multi-tree or not
  
  # Options are:
  # tree ------- a multitree or ape tree read into r
  # traittable - a trait table where rows are tips in the phylogeny and columns 
  #              are trait values. The first column must be a column of rownames
  # outdir ----- an output director to use for the files that will be generated
  
  # Required libraries
  library("data.table")
  library("adephylo")
  library("ape")
  
  
  tree_all <- tree # can be multitree
  table <- traittable
  
  # Sense if multiphylo or not
  if(class(tree_all) == "multiPhylo") {
    loop_length <- length(tree_all)
    # Initialize output matrix
    Mean_all<-matrix(nrow=ncol(table)-1,ncol=loop_length)
  } else {
    loop_length = 1
    Mean_all<-matrix(nrow=ncol(table)-1,ncol=loop_length)
  }
  
  for (m in 1:loop_length) {#loop through all trees
    if(class(tree_all) == "multiPhylo") {
      tree<-tree_all[[m]]
    } else {
      tree <- tree_all
    }
    # testing if table and tree contain the same entries - else drop tips
    z<-subset(tree$tip.label,!(tree$tip.label %in% table[,1]))
    if (length(z) > 0) {
      drop.tip(tree,z)
    }
    
    #rooting tree with first taxon - change if different root
    if(!is.rooted(tree)) {
      root_tree<-root(tree,1,resolve.root=T)
    } else {
      root_tree <- tree
    }
    
    #replacing negative branch lengths - e.g., from PHYLIP
    root_tree$edge.length[root_tree$edge.length <= 0] =  0.00001
    subtree<-subtrees(root_tree, wait=FALSE)
    
    
    cluster_mean<-numeric(length=0)
    
    # loop through all traits
    for (j in 2:ncol(table)) {
      print(c("Analyzing",names(table[j]),"!!!!!!!!!!!!!!!!!!!!!!!"))
      #Loading trait table
      table_tmp<-table[,c(1,j)]
      colnames(table_tmp)[1]<-"ID";
      colnames(table_tmp)[2]<-"Trait";
      
      # removing all entries not in tree 
      table_tmp2<-data.table(table_tmp)
      setkey(table_tmp2,ID)
      table2<-table_tmp2[intersect(table_tmp2$ID,root_tree$tip.label)]
      setkey(table2,ID)
      
      #initializing result vectors and file names
      positives<-vector(mode="list",length=0)
      cluster_size<-numeric(length=0)
      cluster_size_file<-paste(out_dir, "/R_cluster_size_",names(table[j]),".txt",sep="")
      
      cluster_dist<-numeric(length=0)
      cluster_dist_file<-paste(out_dir, "/R_cluster_dist_",names(table[j]),".txt",sep="")
      
      #initalizing files
      if (m == 1) {
        cat(c("trait","tree","distance","cluster_size"), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
        cat("\n", file = cluster_size_file, fill = FALSE, labels = NULL,append = TRUE)
        
        cat(c("trait","tree","distance"), file = cluster_dist_file, sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
        cat("\n", file = cluster_dist_file, fill = FALSE, labels = NULL,append = TRUE)
      }
      
      
      #loop through all subtrees and determining if any subtrees have >90% positives
      for (i in 1:length(subtree)){
        tip_names<-subtree[[i]]$tip.label
        if (mean(table2[tip_names][,Trait]) > 0.9 ) {#change this value if you want a new threshold
          match_test<-match(tip_names,positives)
          if (all(is.na(match_test))) {
            positives<-c(positives,tip_names)
            cluster_dist<-distRoot(subtree[[i]],tip_names, method=c("p"))
            cluster_size<-c(cluster_size,mean(cluster_dist))
            
            # printing to files###
            cat(c(names(table[j]),m,mean(cluster_dist),length(cluster_dist)), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
            cat("\n", file = cluster_size_file, fill = FALSE, labels = NULL,append = TRUE)
            
            cat(names(table[j]),m,cluster_dist, file = cluster_dist_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
            cat("\n", file = cluster_dist_file, fill = FALSE, labels = NULL,append = TRUE)
            
            
            #print(cluster_dist)
          }
          else if (any(is.na(match_test))) {
            print("some NAs - something is weird")
          }
          else {
            #print(tip_names)
            #print("found cluster before")
          }
        }
      }
      
      
      ##### find singletons ######
      a<-table2[table2$Trait == 1,][,ID]
      g<-as.character(a)
      
      singletons_names = setdiff(g,positives)
      if (length(singletons_names) > 0) {
        for (h in 1:length(singletons_names)){
          # weigh singletons with half
          singleton_edges = 0.5*root_tree$edge.length[which.edge(root_tree,singletons_names[h])] #here we use half the distance for singletons
          cluster_size<-c(cluster_size,singleton_edges)
          
          cat(c(names(table[j]),m,singleton_edges,1), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
          cat("\n", file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
        }
        
      }
      Mean_all[j-1,m] = mean(cluster_size)
    }
    
    
  }
  
  #output file
  rownames(Mean_all) <- colnames(table)[2:ncol(table)]
  write.table(Mean_all,paste0(out_dir, "/Mean_all_bootstrap2.txt"), sep = "\t",
              col.names = FALSE, row.names = TRUE, quote = FALSE)
}

