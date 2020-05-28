###########################################################################################
# Script to perform cross validation regardless of algorithm choice.
###########################################################################################
cross_validation <- function(drug_source, views_source, view_combination, algorithm, standardize_any=F, standardize_response=F,
                             parameters, k_fold=5, random=NULL) {
  
  # ****************
  # scripts
  source("R/GLMs.R")
  source("R/mgaussian.R")
  source("R/BEMKL.R")
  source("R/L21_R.R")

  # ****************
  # General variables:
  names_view <- names(view_combination) # Need data type 
  N <- nrow(drug_source) # number of observations
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  P <- length(view_combination) # number of views
  
  # Output variables: 
  output <- list()
  
  # prepare for k-fold CV
  rand_obs_ix <- sample(N) # prepare index of randomized observatios
  folds <- cut(seq(1,N),breaks=k_fold,labels=FALSE) #divide in k folds
  folds_ix <- lapply(1:k_fold, function(i){ 
    rand_obs_ix[which(folds==i,arr.ind=TRUE)]
  })
  
  #if random is set to a number, a randomised cross-validation, selecting each time
  #20% of the samples is performed and repeated the specified amount of times
  if (!is.null(random)){
    k_fold = random
  }
  
  # ****************
  # k-fold CV:
  #if random is set to a number, a randomised cross-validation, selecting each time
  #20% of the samples is performed and repeated the specified amount of times
  
  for (k in 1:k_fold){
    cat("\n---\niteration", k,":", "\n")
    
    # indices of observations used for learning (training) and prediction (test)
    if (is.null(random)){
      learning_indices = setdiff(1:N, folds_ix[[k]]);
      prediction_indices = folds_ix[[k]];
    }else{
      prediction_indices = sample(1:N, round(N*20/100))
      learning_indices = setdiff(1:N, prediction_indices)
    }
    
    cat("Model -->", algorithm, "\n")
    
    if (algorithm %in% c("Lasso", "Elastic_Net", "Logistic_EN")){
      
      # Generalized linear models:
      output[[k]] <- GLMs(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                          standardize_any, standardize_response, parameters = parameters[[algorithm]], iteration = k) 
      
    }else if (algorithm %in% "Multi_Task_EN"){
      
      # Multi-task (Elastic Net):
      output[[k]] <- mgaussian(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                               standardize_any, standardize_response, parameters = parameters[[algorithm]], iteration = k) 
      
    }else if (algorithm %in% "BEMKL"){
      
      # BEMKL:
      output[[k]] <- BEMKL(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                           standardize_any, standardize_response, parameters = parameters[[algorithm]], iteration = k) 
      
    }else if (algorithm %in% "L21"){
      
      # Multi-task (L21):
      output[[k]] <- L21_R(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                           standardize_any, standardize_response, parameters = parameters[[algorithm]], iteration = k) 
    }
  }
  return(output)
}
