#' ## Calculate Random Forest Assembly
#' This step runs random forest models on assembly process

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
library(ranger) # fast implementation of randomforest
library(fastshap) # for shapley values
library(caret)
library(viridis)
library(janitor)
library(here)

# Load required data
#source(here("setup.R"))
source(here("Assembly-analysis", "R", "00_all_MAGs_assembly","03_prepare_pairwise_data.R")) # note this also sources setup.R

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
inputs.fp <- here("Assembly-analysis", "outputs")
outputs.fp <- here("Assembly-analysis", "outputs", "funct_rf_model")
figures.fp <- here("Assembly-analysis", "figures")

# If makes sense, why things might happen, how compares to microbial data

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

assembly_levels <- c("Homogenous.selection", "Heterogenous.selection",
                     "Homogenizing.dispersal", "Dispersal.limitation.and.drift",
                     "Drift")

#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Read in OTU table and env data
input <- input_ra
input$otu_table <- input$otu_table[-1]

# Read in differential data
betanull.lf.diff <- read_csv(file = paste0(inputs.fp, "/betanull.lf.diff.csv"))


#### ====================================================================== ####
# Setup paralell processing
#### ====================================================================== ####
server <- TRUE
retune_models <- FALSE
# because tuning models uses caret which hates mpi, and shapley takes forever, 
# we only run shapley if models have been pre-tuned and we can take advantage of 
# mpi 
if(retune_models) {
  run_shapley = FALSE 
} else {
  run_shapley = TRUE
}
set.seed(50) # set a seed so training and testing will always be the same between runs

# Tell the user their settings:
writeLines(paste("User settings are as follows: \n",
                 "Server:", server,
                 "Retune models? ", retune_models,
                 "Shapley values: ", run_shapley))

#+ include = FALSE
# Set up the paralelization;
# Install packages
#packReq <- c("doMPI")

#Install and load all required packages
#lapply(packReq, function(x) {
#  print(x)
#  if (require(x, character.only = TRUE) == FALSE) {
#    install.packages(x)
#    library(x, character.only = TRUE)
#  }})

# If parallel is true, detect available cores
if(server==T) {
  library(Rmpi)  # R implementation of MPI interface
  library(doMPI)
  n.cores <- max(1, mpi.universe.size()-1) # set n.cores to 1 or mpi allocation size, whichever is bigger
  if(n.cores == 1) {
    n.cores <- parallel::detectCores() -1
  }
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

#### ====================================================================== ####

# Useful functions
#### ====================================================================== ####

# Random forest tunning and prediction function for categorical Y
run_categ_random_forest <- function(data = data.rf,
                                    y_column = "Assembly_Process",
                                    use_weights = TRUE, 
                                    tune_model = TRUE,
                                    tuned_model = NULL,
                                    remove_hi_cor = TRUE,
                                    remove_nzv = TRUE,
                                    shapley = run_shapley,
                                    seed.val = NULL,
                                    rang_n_trees = 10) {
  
  #attach(environment(caret::train)) # necessary to work with  doMPI see: https://github.com/topepo/caret/issues/1017
  
  # Prepare data
  # Remove na values and assign training and testing data (with proportional groups)
  if(!is.null(seed.val)) {
    set.seed(seed.val)
  }
  
  data <- data %>%
    # Remove NAs
    na.omit()
  
  if(remove_nzv) {
    # identifying near zero variance predictors. see: https://rdrr.io/cran/caret/man/nearZeroVar.html 
    nzv <- data %>%
      select(-!!as.name(y_column)) %>%
      select_if(is.numeric) %>% 
      caret::nearZeroVar(freqCut = 15,
                         uniqueCut = 10,
                         saveMetrics = TRUE)
    
    
    col_to_remove <- nzv %>% filter(nzv) %>%
      rownames()
    
    writeLines(paste("Removing columns with near zero variance. Columns removed: "))
    writeLines(col_to_remove)
    
    data <- data %>% 
      select(!all_of(col_to_remove))
    
  }
    
  if(remove_hi_cor) {
    # Idenfying variables with correlation greater than 0.7
    cor_vals <- data %>% 
      select(-!!as.name(y_column)) %>%
      select_if(is.numeric) %>%
      cor()
    
    hicor <- cor_vals %>%
      caret::findCorrelation(cutoff = 0.9, names = TRUE, exact = TRUE)
    
    cor_vals[lower.tri(cor_vals,diag=TRUE)] <- NA
    
    hicor_report <- cor_vals %>% data.frame() %>% 
      rownames_to_column(var = "Var1") %>% 
      pivot_longer(-Var1, names_to = "Var2", values_to = "cor") %>% 
      filter(!is.na(cor)) %>%
      mutate(removed1 = ifelse(Var1 %in% hicor, TRUE, FALSE),
             removed2 = ifelse(Var2 %in% hicor, TRUE, FALSE),
             removed = ifelse(removed1|removed2, TRUE, FALSE))
    
    writeLines(paste("Removing columns with high correlations. Columns removed: "))
    writeLines(hicor)
    
    data <- data %>% 
      select(!all_of(hicor))
  }
  
  # Partition Data into training and testing sets
  data <- data %>%
    # add a training and testing column
    group_by(!!as.name(y_column)) %>%
    mutate(traintest = sample(c(0,1), size = n(), replace = T, prob = c(0.8,0.2))) %>%
    ungroup() %>%
    mutate(traintest = ifelse(traintest == 0, "train", "test"))
  
  
  # Separate into training and testing
  train.rf <- data %>%
    filter(traintest == "train") %>%
    select(-traintest, -Site1, -Site2)
  
  test.rf <- data %>%
    filter(traintest == "test") %>%
    select(-traintest, -Site1, -Site2)
  
  test.site <- data %>% 
    filter(traintest == "test") %>%
    select(all_of(y_column), Site1, Site2, traintest)
  
  # Weights
  if(use_weights) {
    writeLines(paste0("Using weights for ", y_column))
    # Compute weights
    weights.rf <- train.rf %>%
      group_by(!!as.name(y_column)) %>% add_count() %>% ungroup() %>%
      mutate(recip = 1/n) %>%
      select(!!as.name(y_column), recip, n) %>% 
      distinct() %>%
      mutate(weights = recip/sum(recip)) %>%
      select(!!as.name(y_column), weights) %>%
      right_join(train.rf, by = y_column) %>%
      select(!!as.name(y_column), weights)
    
    print(weights.rf %>% unique())
    
    # # Check the number in each group --> we should use weights
    # data %>%
    #   group_by(traintest, !!as.name(y_column)) %>%
    #   tally() %>%
    #   mutate(Percent = 100* n/sum(n)) %>%
    #   rename(Count = n) %>%
    #   ungroup() %>% group_by(traintest) %>% 
    #   mutate(Sum = sum(Count)) %>%
    #   mutate(traintest = paste0(traintest, " (n = ", Sum, ")")) %>%
    #   ggplot(aes_string(y = "Percent", fill = y_column, x = 1,
    #                     group = y_column)) +
    #   theme_bw() #+ 
    #   #barchart_plotting_layers + facet_wrap(~traintest)
  }
  
  if(tune_model) {
    # Set up Model tuning controls
    # cross-validation
    cv <- trainControl(method = "cv",
                       number = 5,
                       verboseIter = TRUE,
                       classProbs = TRUE) # for probability tree
    
    
    # create hyperparameter grid
    n_features <- ncol(train.rf) - 1
    hyper_grid <- expand.grid(
      mtry = c(sqrt(n_features), floor(n_features * c(.05, .15, .25, .333, .4))), # 5 evenly spaced across range 2 to # predictors, including the default sqrt(p) value
      min.node.size = c(1, 3, 5, 10),
      splitrule = c("gini", "extratrees")
    )
    
    if(use_weights) {
      weight_setting <- weights.rf$weights
    } else {
      weight_setting <- NULL
    }
    rang_form <- as.formula(paste0(y_column, "~ ."))
    
    writeLines("Training RF model, columns used: ")
    writeLines(names(train.rf))
    
    # Run the random forest model using tuning to find the best parameters
    rf_tuned <-  caret::train(rang_form,
                                   data=train.rf,
                                   method="ranger",
                                   metric = "Accuracy",
                                   num.trees = rang_n_trees*n_features,
                                   #tuneLength = 4,
                                   tuneGrid = hyper_grid,
                                   trControl = cv,
                                   weights = weight_setting,
                                   importance = 'impurity') 
    
    print(rf_tuned)
    rf_tuned_plot <- plot(rf_tuned) 
    
    tuned_output <- list(tuned_model = rf_tuned,
         model_tunning = rf_tuned$bestTune,
         rf_tuned_plot = rf_tuned_plot,
         test.rf = test.rf,
         train.rf = train.rf,
         weights.rf = weights.rf,
         test.site = test.site)
    
    if(remove_hi_cor) {
      hicor_list <- list(hicor = hicor_report)
      tuned_output <- append(tuned_output, hicor_list)
    }
    
    if(remove_nzv) {
      nzv_list <- list(nzv = nzv)
      tuned_output <- append(tuned_output, nzv_list)
    }
    writeLines("i'm in tuned output")
    
    return(tuned_output)
  } else {
    # run tuned model
    # requires that tuned_model != NULL
    if(is.null(tuned_model)) {
      stop("tuned_model is not provided")
    }
    writeLines("I'm in hard-code model")
    
    # Clean names for ranger, which is more fussy about that sort of thing
    data_rang <- janitor::clean_names(tuned_model$train.rf, case = "none")
    
    n_features <- ncol(data_rang) - 1
    
    rang_form <- as.formula(paste0(y_column, "~ ."))
    
    writeLines("Running tuned RF model [impurity], columns used: ")
    writeLines(names(train.rf))
    
    rf_model_impurity <-  ranger(rang_form,
                                  data=data_rang,
                                  num.trees = rang_n_trees*n_features,
                                  mtry = tuned_model$bestTune$mtry,
                                  splitrule = tuned_model$bestTune$splitrule,
                                  min.node.size = tuned_model$bestTune$min.node.size,
                                  case.weights = tuned_model$weights.rf$weights,
                                  importance = 'impurity',
                                  probability = TRUE,
                                  seed = seed.val)

    writeLines("Running tuned RF model [permutation], columns used: ")
    writeLines(names(train.rf))
    
    rf_model_permutation <-  ranger(rang_form,
                                 data=data_rang,
                                 num.trees = rang_n_trees*n_features,
                                 mtry = tuned_model$bestTune$mtry,
                                 splitrule = tuned_model$bestTune$splitrule,
                                 min.node.size = tuned_model$bestTune$min.node.size,
                                 case.weights = tuned_model$weights.rf$weights,
                                 importance = 'permutation',
                                 probability = TRUE,
                                 seed = seed.val)
    
    # Measure Goodness of Fit:
    # Predict and evaluate model performance
    pred_impure <- predict(rf_model_impurity, tuned_model$test.rf)
    pred_permute <- predict(rf_model_permutation, tuned_model$test.rf)
    # check <- bind_cols(tuned_model$test.rf, pred) %>% 
    #   rename(pred = last_col())
    
    
    writeLines("Computing variable importance p-values")
    vip_permute_pval <- importance_pvalues(rf_model_permutation, method = "altmann", formula = rang_form,
                                           data=data_rang)
    
    vip_imp_pval <- importance_pvalues(rf_model_impurity, method = "altmann", formula = rang_form,
                                       data=data_rang)
    
    
    vip_impurity <- vip_imp_pval %>%
      data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      rename(Importance = importance)
    vip_permutation <- vip_permute_pval %>%
      data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      rename(Importance = importance)
    
    writeLines("Done with importance p-values")
    # Get shap values
    # For classification problem, we're using probability trees
    if(shapley) {
      writeLines("Computing variable importance - shapley")
      predictors <- data_rang %>%
        select(!all_of(y_column)) %>% 
        as.data.frame()
      
      predictors_test <- tuned_model$train.rf %>%
        select(!all_of(y_column)) %>%
        as.data.frame()
      
      #pfun <- function(object) predict(object, data = tuned_model$test.rf)$predictions
      class_options <- levels(data_rang[[y_column]])
      
      pfun <- function(object, newdata) {
        require(ranger)
        predict(object, data = newdata)$predictions[,whichClass]
      }
      
      for( whichClass in seq_along(class_options)) {
        categ_name <- class_options[whichClass]
        pfun <- function(object, newdata) {
          require(ranger)
          predict(object, data = newdata)$predictions[,whichClass]
        }
        shap_vals_imp_temp <-  fastshap::explain(rf_model_impurity, X = predictors, 
                                            pred_wrapper = pfun, nsim = 1000,
                                            newdata = predictors_test,
                                            .parallel = T,
                                            .paropts = list(.packages = c("ranger", "fastshap"),
                                                            .verbose = TRUE))
        print(dim(shap_vals_imp_temp))
        names(shap_vals_imp_temp) <- paste0(names(shap_vals_imp_temp), "_", categ_name)
        if(!exists("shap_vals_imp")) {
          shap_vals_imp <- shap_vals_imp_temp
        } else {
          shap_vals_imp <- cbind(shap_vals_imp, shap_vals_imp_temp)
        }
      }
      
      for( whichClass in seq_along(class_options)) {
        categ_name <- class_options[whichClass]
        pfun <- function(object, newdata) {
          require(ranger)
          predict(object, data = newdata)$predictions[,whichClass]
        }
        print(whichClass)
        shap_vals_perm_temp <-  fastshap::explain(rf_model_permutation, X = predictors, 
                                                 pred_wrapper = pfun, nsim = 1000,
                                                 newdata = predictors_test,
                                                 .parallel = T,
                                                 .paropts = list(.packages = c("ranger", "fastshap"),
                                                                 .verbose = TRUE))
        print(dim(shap_vals_perm_temp))
        names(shap_vals_perm_temp) <- paste0(names(shap_vals_perm_temp), "_", categ_name)
        if(!exists("shap_vals_perm")) {
          shap_vals_perm <- shap_vals_perm_temp
        } else {
          shap_vals_perm <- cbind(shap_vals_perm, shap_vals_perm_temp)
        }
      }
    }
    
    
    # # Calculate variable importance
    # vip_impurity <- vip::vi(rf_model_impurity, target = y_column, pred_wrapper = pfun)
    # vip_permutation <- vip::vi(rf_model_permutation, target = y_column, pred_wrapper = pfun)
    
    # Calculate confusion matrix; Take the column with maximum probability as the predicted value
    writeLines("Cacluating Confusion Matrix")
    pred_impure_max_prob <- colnames(pred_impure$predictions)[max.col(pred_impure$predictions)] %>%
      factor(levels = levels(tuned_model[["test.rf"]][[y_column]]))
    
    pred_permute_max_prob <- colnames(pred_permute$predictions)[max.col(pred_permute$predictions)] %>%
      factor(levels = levels(tuned_model[["test.rf"]][[y_column]]))
    
    conf_mat_impure <- confusionMatrix(pred_impure_max_prob, tuned_model[["test.rf"]][[y_column]])
    conf_mat_permute <- confusionMatrix(pred_permute_max_prob, tuned_model[["test.rf"]][[y_column]])
    
    return_list <- 
      list(rf_impurity = rf_model_impurity,
           rf_permutation = rf_model_permutation,
           vip_impurity = vip_impurity,
           vip_permutation = vip_permutation,
           conf_mat_impure = conf_mat_impure,
           conf_mat_permute = conf_mat_permute,
           test_dat = tuned_model$test.rf,
           train_dat = tuned_model$train.rf,
           test_site = tuned_model$test.site)
    
    if(shapley) {
      shap_list <- list( shap_vals_imp = shap_vals_imp,
            shap_vals_perm = shap_vals_perm)
      return_list <- append(return_list, shap_list)
    }
    
    if(remove_hi_cor) {
      hicor_list <- list(hicor = hicor_report)
      return_list <- append(return_list, hicor_list)
    }
    
    if(remove_nzv) {
      nzv_list <- list(nzv = nzv)
      return_list <- append(return_list, nzv_list)
    }
    return(return_list)
  }
  
}

# Random forest tunning and prediction function for continuous Y
run_cont_random_forest <- function(data = data.rf,
                                    y_column = "BetaNTI",
                                    tune_model = TRUE,
                                    tuned_model = NULL,
                                    remove_hi_cor = TRUE,
                                    remove_nzv = TRUE,
                                    shapley = run_shapley,
                                    seed.val = NULL,
                                    rang_n_trees = 10) { # rang_n_trees = multiplier of trees * n_features
  #attach(environment(caret::train)) # necessary to work with  doMPI see: https://github.com/topepo/caret/issues/1017

  # Prepare data
  # Remove na values and assign training and testing data (with proportional groups)
  if(!is.null(seed.val)) {
    set.seed(seed.val)
  }
  
  data <- data %>%
    # Remove NAs
    na.omit()
  
  if(remove_nzv) {
    # identifying near zero variance predictors. see: https://rdrr.io/cran/caret/man/nearZeroVar.html 
    nzv <- data %>%
      select(-!!as.name(y_column)) %>%
      select_if(is.numeric) %>% 
      caret::nearZeroVar(freqCut = 15,
                         uniqueCut = 10,
                         saveMetrics = TRUE)
    
    
    col_to_remove <- nzv %>% filter(nzv) %>%
      rownames()
    
    writeLines(paste("Removing columns with near zero variance. Columns removed: "))
    writeLines(col_to_remove)
    
    data <- data %>% 
      select(!all_of(col_to_remove))
    
  }
  
  if(remove_hi_cor) {
    # Idenfying variables with correlation greater than 0.7
    cor_vals <- data %>% 
      select(-!!as.name(y_column)) %>%
      select_if(is.numeric) %>%
      cor()
    
    hicor <- cor_vals %>%
      caret::findCorrelation(cutoff = 0.9, names = TRUE, exact = TRUE)
    
    cor_vals[lower.tri(cor_vals,diag=TRUE)] <- NA
    
    hicor_report <- cor_vals %>% data.frame() %>% 
      rownames_to_column(var = "Var1") %>% 
      pivot_longer(-Var1, names_to = "Var2", values_to = "cor") %>% 
      filter(!is.na(cor)) %>%
      mutate(removed1 = ifelse(Var1 %in% hicor, TRUE, FALSE),
             removed2 = ifelse(Var2 %in% hicor, TRUE, FALSE),
             removed = ifelse(removed1|removed2, TRUE, FALSE))
    
    writeLines(paste("Removing columns with high correlations. Columns removed: "))
    writeLines(hicor)
    
    data <- data %>% 
      select(!all_of(hicor))
  }
  
  # Partition Data into training and testing sets
  data <- data %>%
    # add a training and testing column
    mutate(traintest = sample(c(0,1), size = n(), replace = T, prob = c(0.8,0.2))) %>%
    ungroup() %>%
    mutate(traintest = ifelse(traintest == 0, "train", "test"))
  
  
  # Separate into training and testing
  train.rf <- data %>%
    filter(traintest == "train") %>%
    select(-traintest, -Site1, -Site2)
  
  test.rf <- data %>%
    filter(traintest == "test") %>%
    select(-traintest, -Site1, -Site2)
  
  test.site <- data %>% 
    filter(traintest == "test") %>%
    select(all_of(y_column), Site1, Site2, traintest)
  
  if(tune_model) {
    # Set up Model tuning controls
    # cross-validation
    cv <- trainControl(method = "cv",
                       number = 5,
                       verboseIter = TRUE)
    
    
    # create hyperparameter grid
    n_features <- ncol(train.rf) - 1
    hyper_grid <- expand.grid(
      mtry = c(sqrt(n_features), floor(n_features * c(.05, .15, .25, .333, .4))), # 5 evenly spaced across range 2 to # predictors, including the default sqrt(p) value
      min.node.size = c(1, 3, 5, 10),
      splitrule = c("variance", "extratrees")
    )
    
    rang_form <- as.formula(paste0(y_column, "~ ."))
    
    writeLines("Training RF model, columns used: ")
    writeLines(names(train.rf))
    
    # Run the random forest model using tuning to find the best parameters
    rf_tuned <-  caret::train(rang_form,
                              data=train.rf,
                              method="ranger",
                              metric = "RMSE",
                              num.trees = rang_n_trees*n_features,
                              tuneGrid = hyper_grid,
                              trControl = cv,
                              importance = 'impurity') 
    
    print(rf_tuned)
    rf_tuned_plot <- plot(rf_tuned) 
    
    tuned_output <- list(tuned_model = rf_tuned,
                         model_tunning = rf_tuned$bestTune,
                         rf_tuned_plot = rf_tuned_plot,
                         test.rf = test.rf,
                         train.rf = train.rf,
                         test.site = test.site)
    if(remove_hi_cor) {
      hicor_list <- list(hicor = hicor_report)
      tuned_output <- append(tuned_output, hicor_list)
    }
    
    if(remove_nzv) {
      nzv_list <- list(nzv = nzv)
      tuned_output <- append(tuned_output, nzv_list)
    }
    writeLines("i'm in tuned output")
    
    return(tuned_output)
  } else {
    # run tuned model
    # requires that tuned_model != NULL
    if(is.null(tuned_model)) {
      stop("tuned_model is not provided")
    }
    writeLines("I'm in hard-code model")
    
    # Clean names for ranger, which is more fussy about that sort of thing
    data_rang <- janitor::clean_names(tuned_model$train.rf, case = "none")
    
    n_features <- ncol(data_rang) - 1
    
    rang_form <- as.formula(paste0(y_column, "~ ."))
    
    writeLines("Running tuned RF model [impurity], columns used: ")
    writeLines(names(train.rf))
    
    rf_model_impurity <-  ranger(rang_form,
                                 data=data_rang,
                                 num.trees = rang_n_trees*n_features,
                                 mtry = tuned_model$bestTune$mtry,
                                 splitrule = tuned_model$bestTune$splitrule,
                                 min.node.size = tuned_model$bestTune$min.node.size,
                                 importance = 'impurity',
                                 seed = seed.val)
    
    writeLines("Running tuned RF model [permutation], columns used: ")
    writeLines(names(train.rf))
    
    rf_model_permutation <-  ranger(rang_form,
                                    data=data_rang,
                                    num.trees = rang_n_trees*n_features,
                                    mtry = tuned_model$bestTune$mtry,
                                    splitrule = tuned_model$bestTune$splitrule,
                                    min.node.size = tuned_model$bestTune$min.node.size,
                                    importance = 'permutation',
                                    seed = seed.val)
    # Measure Goodness of Fit:
    # Predict and evaluate model performance
    pred_impure <- predict(rf_model_impurity, tuned_model$test.rf)
    pred_permute <- predict(rf_model_permutation, tuned_model$test.rf)
    # check <- bind_cols(tuned_model$test.rf, pred) %>% 
    #   rename(pred = last_col())
    
    writeLines("Computing variable importance p-values")
    vip_permute_pval <- importance_pvalues(rf_model_permutation, method = "altmann", formula = rang_form,
                                   data=data_rang)
    
    vip_imp_pval <- importance_pvalues(rf_model_impurity, method = "altmann", formula = rang_form,
                                           data=data_rang)
    
    
    vip_impurity <- vip_imp_pval %>%
      data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      rename(Importance = importance)
    vip_permutation <- vip_permute_pval %>%
      data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      rename(Importance = importance)
    
    # Get shap values
    if(shapley) {
      writeLines("Computing variable importance - shapley")
      predictors <- data_rang %>%
        select(!all_of(y_column)) %>% 
        as.data.frame()
      
      predictors_test <- tuned_model$train.rf %>%
        select(!all_of(y_column)) %>%
        as.data.frame()
      
      #pfun <- function(object) predict(object, data = tuned_model$test.rf)$predictions
      pfun <- function(object, newdata) {
        predict(object, data = newdata)$predictions
      }
      
      shap_vals_imp <-  fastshap::explain(rf_model_impurity, X = predictors, 
                                          pred_wrapper = pfun, nsim = 1000,
                                          newdata = predictors_test, .parallel = T,
                                          .paropts = list(.packages = c("ranger", "fastshap"),
                                                          .verbose = TRUE))
      shap_vals_perm <-  fastshap::explain(rf_model_permutation, X = predictors, 
                                           pred_wrapper = pfun, nsim = 1000,
                                           newdata = predictors_test, .parallel = T,
                                           .paropts = list(.packages = c("ranger", "fastshap"),
                                                           .verbose = TRUE))
    }
    
    
    r.rap <- tuned_model$test.rf[y_column] - pred_permute$predictions
    rmse_model <- sqrt(sum(r.rap^2)/length(r.rap))
    
    #rmse_impure <- RMSE(pred_impure,test.rf[y_column])
    #rmse_permute <- RMSE(pred_impure,test.rf[y_column])
    
    
    return_list <- 
      list(rf_impurity = rf_model_impurity,
           rf_permutation = rf_model_permutation,
           rmse_model = rmse_model,
           pred_impure = pred_impure,
           pred_permute = pred_permute,
           test_dat = tuned_model$test.rf,
           test_site = tuned_model$test.site,
           train_dat = tuned_model$train.rf,
           vip_impurity = vip_impurity,
           vip_permutation = vip_permutation)
    
    if(shapley) {
      shap_list <- list( shap_vals_imp = shap_vals_imp,
                         shap_vals_perm = shap_vals_perm)
      return_list <- append(return_list, shap_list)
    }
    
    if(remove_hi_cor) {
      hicor_list <- list(hicor = hicor_report)
      return_list <- append(return_list, hicor_list)
    }
    
    if(remove_nzv) {
      nzv_list <- list(nzv = nzv)
      return_list <- append(return_list, nzv_list)
    }
    return(return_list)
  }
  
}

#### ====================================================================== ####

# Run random forest models with Biotic Side (aka, function)
#### ====================================================================== ####
# Remove na values and assign training and testing data (with proportional groups)
data.rf <- betanull.lf.diff %>% 
  select(Assembly_Process, Habitat_comp, Depth_diff, T_air_7d_diff, 
         pct_time_below_WTD_21d_diff, Year_diff, 
         contains("path_diff_abund"), 
         Site1, Site2) %>% 
  select(!contains("scl_")) %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>% 
  mutate(Assembly_Process = gsub(" ", ".", Assembly_Process)) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels))  %>% 
  filter(Habitat_comp %in% c("Fen_Fen", "Bog_Bog", "Palsa_Palsa")) %>%
  mutate(Habitat_comp = droplevels(Habitat_comp)) %>%
  filter(Depth_diff < 10) # Only use depth differences in the same depthLumping 
names_translation_data.rf <- data.frame(original = names(data.rf)) 
data.rf <- data.rf %>%
  janitor::clean_names(case = "none")

names_translation_data.rf$new <- names(data.rf)

# Prepare matrices
writeLines("Running all categorical model...")

# Should we retune the model if an already tuned model file exits?
if(!retune_models & file.exists(here(outputs.fp, "all_tune_model.RDS"))) {
  all_tune_model <- readRDS(here(outputs.fp, "all_tune_model.RDS")) 
} else {
  all_tune_model <- run_categ_random_forest(data = data.rf,
                                            y_column = "Assembly_Process")
  saveRDS(all_tune_model, here(outputs.fp, "all_tune_model.RDS"))
}

all_hardcode_model <- run_categ_random_forest(data = data.rf,
                                      tuned_model = all_tune_model,
                                      tune_model = FALSE,
                                      y_column = "Assembly_Process")


all_hardcode_model$vip_permutation
all_hardcode_model$conf_mat_permute

saveRDS(all_hardcode_model, paste0(outputs.fp, "/all_hardcode_model.RDS"))

#### ====================================================================== ####

#### ====================================================================== ####
# palsa
palsa_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp) %>%
  mutate(Assembly_Process = fct_drop(Assembly_Process))

# Should we retune the model
if(!retune_models & file.exists(here(outputs.fp, "palsa_model.RDS"))) {
  palsa_model <- readRDS(here(outputs.fp, "palsa_model.RDS")) 
} else {
  palsa_model <- run_categ_random_forest(data = palsa_data,
                                         y_column = "Assembly_Process")
  
  plot(palsa_model$tuned_model)
  varImp(palsa_model$tuned_model)
  
  saveRDS(palsa_model, here(outputs.fp, "palsa_model.RDS"))
}


writeLines("Running palsa categorical model...")
palsa_hardcode_model <- run_categ_random_forest(data = palsa_data,
                                          tuned_model = palsa_model,
                                          tune_model = FALSE,
                                          rang_n_trees = 50,
                                          y_column = "Assembly_Process")



palsa_hardcode_model$vip_permutation 
palsa_hardcode_model$conf_mat_permute

saveRDS(palsa_hardcode_model, paste0(outputs.fp, "/palsa_hardcode_model.RDS"))

# Bog 
bog_data <- data.rf %>% filter(Habitat_comp == "Bog_Bog") %>%
  select(-Habitat_comp) %>% 
  filter(Assembly_Process != "Heterogenous selection") %>% # only 1 row
  mutate(Assembly_Process = fct_drop(Assembly_Process))

if(!retune_models & file.exists(here(outputs.fp, "bog_model.RDS"))) {
  bog_model <- readRDS(here(outputs.fp, "bog_model.RDS")) 
} else {
  bog_model <- run_categ_random_forest(data = bog_data,
                                       rang_n_trees = 60,
                                       y_column = "Assembly_Process")
  
  plot(bog_model$tuned_model)
  varImp(bog_model$tuned_model)
  
  saveRDS(bog_model, here(outputs.fp, "bog_model.RDS"))
}

writeLines("Running bog categorical model...")
bog_hardcode_model <- run_categ_random_forest(data = bog_data,
                                          tuned_model = bog_model,
                                          tune_model = FALSE,
                                          rang_n_trees = 60,
                                          y_column = "Assembly_Process")
bog_hardcode_model$rf_permutation
bog_hardcode_model$conf_mat_permute

saveRDS(bog_hardcode_model, paste0(outputs.fp, "/bog_hardcode_model.RDS"))

# Fen

fen_data <- data.rf %>% filter(Habitat_comp == "Fen_Fen") %>%
  select(-Habitat_comp) %>% 
  mutate(Assembly_Process = fct_drop(Assembly_Process))

if(!retune_models & file.exists(here(outputs.fp, "fen_model.RDS"))) {
  fen_model <- readRDS(here(outputs.fp, "fen_model.RDS")) 
} else {
  fen_model <- run_categ_random_forest(data = fen_data,
                                       y_column = "Assembly_Process")
  
  plot(fen_model$tuned_model)
  varImp(fen_model$tuned_model)
  
  saveRDS(fen_model, here(outputs.fp, "fen_model.RDS"))
}

writeLines("Running fen categorical model...")

fen_hardcode_model <- run_categ_random_forest(data = fen_data,
                                                tuned_model = fen_model,
                                                tune_model = FALSE,
                                                rang_n_trees = 10,
                                                y_column = "Assembly_Process")


fen_hardcode_model$rf_permutation
fen_hardcode_model$conf_mat_permute
fen_hardcode_model$vip_permutation

saveRDS(fen_hardcode_model, paste0(outputs.fp, "/fen_hardcode_model.RDS"))

#### ====================================================================== ####


# Run random forest models with Deterministic predictor
#### ====================================================================== ####
# Remove na values and assign training and testing data (with proportional groups)
data.rf <- betanull.lf.diff %>% 
  select(BetaNTI, Habitat_comp, Depth_diff, T_air_7d_diff, 
         pct_time_below_WTD_21d_diff, Year_diff, contains("_path_diff_abund"), Site1, Site2) %>% 
  select(BetaNTI, Habitat_comp, Depth_diff, Year_diff, !contains("scl_"), Site1, Site2) %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>% 
  filter(Habitat_comp %in% c("Fen_Fen", "Bog_Bog", "Palsa_Palsa")) %>%
  filter(Depth_diff < 10) %>% # Only use depth differences in the same depthLumping %>%
  janitor::clean_names(case = "none")

# palsa
palsa_det_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp) #%>%
  #select(!contains("path_diff_abund"))

writeLines("Running palsa deterministic model...")

if(!retune_models & file.exists(here(outputs.fp, "palsa_det_model.RDS"))) {
  palsa_det_model <- readRDS(here(outputs.fp, "palsa_det_model.RDS")) 
} else {
  palsa_det_model <- run_cont_random_forest(data = palsa_det_data,
                                            seed.val = 50,
                                            y_column = "BetaNTI",
                                            rang_n_trees = 50)  
  plot(palsa_det_model$tuned_model)
  varImp(palsa_det_model$tuned_model)
  
  saveRDS(palsa_det_model, here(outputs.fp, "palsa_det_model.RDS"))
}

palsa_hardcode_det_model <- run_cont_random_forest(data = palsa_det_data,
                                                tuned_model = palsa_det_model,
                                                tune_model = FALSE,
                                                rang_n_trees = 50,
                                                y_column = "BetaNTI")


palsa_hardcode_det_model$rmse_model
palsa_hardcode_det_model$rf_permutation
palsa_hardcode_det_model$vip_permutation

saveRDS(palsa_hardcode_det_model, paste0(outputs.fp, "/palsa_hardcode_det_model.RDS"))

# Bog
bog_det_data <- data.rf %>% 
  filter(Habitat_comp == "Bog_Bog") %>%
  select(-Habitat_comp)

writeLines("Running bog deterministic model...")

if(!retune_models & file.exists(here(outputs.fp, "bog_det_model.RDS"))) {
  bog_det_model <- readRDS(here(outputs.fp, "bog_det_model.RDS")) 
} else {
  bog_det_model <- run_cont_random_forest(data = bog_det_data,
                                            seed.val = 50,
                                            y_column = "BetaNTI",
                                            rang_n_trees = 50)  
  plot(bog_det_model$tuned_model)
  varImp(bog_det_model$tuned_model)
  
  saveRDS(bog_det_model, here(outputs.fp, "bog_det_model.RDS"))
}

bog_hardcode_det_model <- run_cont_random_forest(data = bog_det_data,
                                             tuned_model = bog_det_model,
                                             tune_model = FALSE, 
#                                             remove_hi_cor = FALSE,
                                             rang_n_trees = 50,
                                             seed.val = 50,
                                             y_column = "BetaNTI")


bog_hardcode_det_model$rmse_model
bog_hardcode_det_model$rf_permutation
bog_det_model$tuned_model$finalModel

saveRDS(bog_hardcode_det_model, paste0(outputs.fp, "/bog_hardcode_det_model.RDS"))

# Fen
fen_det_data <- data.rf %>% filter(Habitat_comp == "Fen_Fen") %>%
  select(-Habitat_comp)

writeLines("Running fen deterministic model...")


if(!retune_models & file.exists(here(outputs.fp, "fen_det_model.RDS"))) {
  fen_det_model <- readRDS(here(outputs.fp, "fen_det_model.RDS")) 
} else {
  fen_det_model <- run_cont_random_forest(data = fen_det_data,
                                            seed.val = 50,
                                            y_column = "BetaNTI",
                                            rang_n_trees = 50)  
  plot(fen_det_model$tuned_model)
  varImp(fen_det_model$tuned_model)
  
  saveRDS(fen_det_model, here(outputs.fp, "fen_det_model.RDS"))
}

fen_hardcode_det_model <- run_cont_random_forest(data = fen_det_data,
                                             tuned_model = fen_det_model,
                                             tune_model = FALSE,
                                             rang_n_trees = 50,
                                             y_column = "BetaNTI")


fen_hardcode_det_model$rmse_model
fen_hardcode_det_model$rf_permutation
fen_det_model$tuned_model$finalModel

saveRDS(fen_hardcode_det_model, paste0(outputs.fp, "/fen_hardcode_det_model.RDS"))

#### ====================================================================== ####


# Run random forest models with Stochastic predictor
#### ====================================================================== ####
# Remove na values and assign training and testing data (with proportional groups)
data.rf <- betanull.lf.diff %>% 
  select(RCBC.nona, Habitat_comp, Depth_diff, T_air_7d_diff, 
         pct_time_below_WTD_21d_diff, Year_diff, contains("_path_diff_abund"), Site1, Site2) %>% 
  select(RCBC.nona, Habitat_comp, Depth_diff, Year_diff, !contains("scl_"), Site1, Site2) %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>% 
  filter(Habitat_comp %in% c("Fen_Fen", "Bog_Bog", "Palsa_Palsa")) %>%
  filter(Depth_diff < 10) %>% # Only use depth differences in the same depthLumping %>%
  janitor::clean_names(case = "none")


# Palsa
palsa_stoch_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp)

writeLines("Running palsa stochastic model...")

if(!retune_models & file.exists(here(outputs.fp, "palsa_stoch_model.RDS"))) {
  palsa_stoch_model <- readRDS(here(outputs.fp, "palsa_stoch_model.RDS")) 
} else {
  palsa_stoch_model <- run_cont_random_forest(data = palsa_stoch_data,
                                            seed.val = 50,
                                            y_column = "RCBC_nona",
                                            rang_n_trees = 50)  
  plot(palsa_stoch_model$tuned_model)
  varImp(palsa_stoch_model$tuned_model)
  
  saveRDS(palsa_stoch_model, here(outputs.fp, "palsa_stoch_model.RDS"))
}

palsa_hardcode_stoch_model <- run_cont_random_forest(data = palsa_stoch_data,
                                               tuned_model = palsa_stoch_model,
                                               tune_model = FALSE,
                                               rang_n_trees = 50,
                                               y_column = "RCBC_nona")


palsa_hardcode_stoch_model$rmse_model
palsa_hardcode_stoch_model$rf_permutation
palsa_stoch_model$tuned_model$finalModel

saveRDS(palsa_hardcode_stoch_model, paste0(outputs.fp, "/palsa_hardcode_stoch_model.RDS"))

# Bog
bog_stoch_data <- data.rf %>% filter(Habitat_comp == "Bog_Bog") %>%
  select(-Habitat_comp)

writeLines("Running bog stochastic model...")

if(!retune_models & file.exists(here(outputs.fp, "bog_stoch_model.RDS"))) {
  bog_stoch_model <- readRDS(here(outputs.fp, "bog_stoch_model.RDS")) 
} else {
  bog_stoch_model <- run_cont_random_forest(data = bog_stoch_data,
                                              seed.val = 50,
                                              y_column = "RCBC_nona",
                                              rang_n_trees = 50)  
  plot(bog_stoch_model$tuned_model)
  varImp(bog_stoch_model$tuned_model)
  
  saveRDS(bog_stoch_model, here(outputs.fp, "bog_stoch_model.RDS"))
}

bog_hardcode_stoch_model <- run_cont_random_forest(data = bog_stoch_data,
                                             tuned_model = bog_stoch_model,
                                             tune_model = FALSE,
                                             y_column = "RCBC_nona")


bog_hardcode_stoch_model$rmse_model
bog_hardcode_stoch_model$rf_permutation
bog_stoch_model$tuned_model$finalModel

saveRDS(bog_hardcode_stoch_model, paste0(outputs.fp, "/bog_hardcode_stoch_model.RDS"))

# Fen
fen_stoch_data <- data.rf %>% 
  filter(Habitat_comp == "Fen_Fen") %>%
  select(-Habitat_comp)

writeLines("Running fen stochastic model...")


if(!retune_models & file.exists(here(outputs.fp, "fen_stoch_model.RDS"))) {
  fen_stoch_model <- readRDS(here(outputs.fp, "fen_stoch_model.RDS")) 
} else {
  fen_stoch_model <- run_cont_random_forest(data = fen_stoch_data,
                                            seed.val = 50,
                                            y_column = "RCBC_nona",
                                            rang_n_trees = 50)  
  plot(fen_stoch_model$tuned_model)
  varImp(fen_stoch_model$tuned_model)
  
  saveRDS(fen_stoch_model, here(outputs.fp, "fen_stoch_model.RDS"))
}

fen_hardcode_stoch_model <- run_cont_random_forest(data = fen_stoch_data,
                                             tuned_model = fen_stoch_model,
                                             tune_model = FALSE,
                                             rang_n_trees = 50,
                                             y_column = "RCBC_nona")


fen_hardcode_stoch_model$rmse_model
fen_hardcode_stoch_model$rf_permutation
fen_stoch_model$tuned_model$finalModel

saveRDS(fen_hardcode_stoch_model, paste0(outputs.fp, "/fen_hardcode_stoch_model.RDS"))

#### ====================================================================== ####


#### Cleanup parallel
#### ====================================================================== ####

if(server == TRUE) {
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
