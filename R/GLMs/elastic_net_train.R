###########################################################################################
# Script to estimate the hyperparameters, specific for elastic net regression.
###########################################################################################

elastic_net_train <- function(drug_source, views_source, view_combination, standardize_any = F, standardize_response = F,
                                family, parameters, measure_type = "mse", k_fold = 10, parallelize = F, iteration) {
  
  # ****************
  # packages
  library(glmnet)
  library(parallel)
  #library(glmnet,lib.loc= Sys.getenv("R_LIBS")) # cluster
  
  # ****************
  # General variables:
  names_view <- view_combination # Need data type 
  N <- nrow(drug_source) # number of observations
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  P <- length(names_view) # number of views
  alpha_range <- parameters$alpha
  
  # Normalization training set: Z score
  ## View:
  views_source.Z <- views_source
  ## Response: 
  drug_source.Z <- drug_source

  # ****************
  # Hyperparameters estimation: alpha and lambda
  
  # initialize variables:
  # Coefficients and mse values for each iteration
  features.glmnet.1se <- features.glmnet.min <- lowest.mse.cv.1se <- lowest.mse.cv.min <- NULL
  
  # data (mse, sd, coef) for each iteration, for every alpha
  data_task_cv.glmnet <- NULL
  
  # Parameters for each iteration 
  parameters.glmnet.1se <- parameters.glmnet.min <- NULL

  coef_for_each_alpha <- NULL
  MSE_for_each_alpha <- NULL
  sd_for_each_alpha <- NULL
  
  # Internal cross validation: number of folds by default is 10
  assign("fit",lapply(alpha_range, function(X){
    # Same foldid for each alpha, but different foldid per iteration.
    set.seed(3000 + 10 * iteration)
    # Applying internal cross-validation:
    cv.glmnet(x = views_source.Z, 
              y = drug_source.Z, 
              family = family,
              alpha = X, 
              type.measure = measure_type, 
              standardize = FALSE, 
              intercept = TRUE, 
              keep = TRUE,
              nfolds = k_fold,
              parallel = parallelize)}))
  
  names(fit) <- paste0(alpha_range,"-")
  
  # For each iteration, MSE, sd and coefficients for each alpha across range of lambdas. 
  # MSE
  MSE_for_each_alpha <- lapply(fit, function(X) {
    tmp = X$cvm
    names(tmp) = X$lambda
    return(tmp)})
  # SD
  sd_for_each_alpha <- lapply(fit, function(X) {
    tmp = X$cvsd
    names(tmp) = X$lambda
    return(tmp)})
  # COEF
  coef_for_each_alpha <- lapply(fit, function(X) {
    tmp = as.matrix(coef(X, s = X$lambda))
    colnames(tmp) = X$lambda
    return(tmp)})
  # Saving values -->
  data_task_cv.glmnet <- list(MSE_for_each_alpha, sd_for_each_alpha, coef_for_each_alpha)
  names(data_task_cv.glmnet) <- c("MSE","SD","coef")

  # Pair of lambda-alpha values:
  value_MSE_alphas <- NULL
  
  # Find alpha-lambda: 1. min MSE, 2. min MSE + 1SE
  value_MSE_alphas <- unlist(data_task_cv.glmnet$MSE)
  all_alphas = sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), head, 1)
  all_lambdas = sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), function(X) return(X[2]))
  
  ## 1. Lambda-alpha value leading to min MSE
  min_MSE <- value_MSE_alphas[which(value_MSE_alphas == min(value_MSE_alphas,na.rm = TRUE), arr.ind=TRUE)]
  alpha_min_MSE = sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), head, 1)
  lambda_min_MSE = sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), function(X) return(X[2]))
  
  ## 2.  Lambda-alpha value leading leading to 1 SE + min MSE 
  one.se <- data_task_cv.glmnet$SD[[paste0(alpha_min_MSE,"-")]][which(names(data_task_cv.glmnet$SD[[paste0(alpha_min_MSE,"-")]]) == lambda_min_MSE)]
  high <- min_MSE + one.se
  interval_mse <- value_MSE_alphas[value_MSE_alphas < high]
  interval_mse <- interval_mse[which(as.numeric(sapply(strsplit(names(interval_mse), split = "-.", fixed = TRUE), function(X) return(X[2])))
                                   >= as.numeric(lambda_min_MSE))]
  one.se_MSE <- interval_mse[which(interval_mse == max(interval_mse))]
  lambda_1se_MSE <- sapply(strsplit(names(one.se_MSE), split = "-.", fixed = TRUE), function(X) return(X[2]))
  alpha_1se_MSE <- sapply(strsplit(names(one.se_MSE), split = "-.", fixed = TRUE), head, 1)
 
  # Saving parameters for each iteration:
  parameters.glmnet.1se <- list(alpha_1se_MSE, lambda_1se_MSE)  
  parameters.glmnet.min <- list(alpha_min_MSE, lambda_min_MSE)
  names(parameters.glmnet.1se) <- names(parameters.glmnet.min) <- c("alpha", "lambda")

  # 3. Selecting coefficients for both pairs of alpha-lambda
  # Alpha-lambda min
  coef_min <- as.matrix(data_task_cv.glmnet$coef[[paste0(alpha_min_MSE,"-")]][,lambda_min_MSE])
  features.glmnet.min <- coef_min ; names(features.glmnet.min) <- rownames(coef_min)
  # Alpha-lambda 1se
  coef_1se <- as.matrix(data_task_cv.glmnet$coef[[paste0(alpha_1se_MSE,"-")]][,lambda_1se_MSE])
  features.glmnet.1se <- coef_1se ; names(features.glmnet.1se) <- rownames(coef_1se)
  # mse.1se and mse.min from internal cross validation for iteration k
  lowest.mse.cv.1se = one.se_MSE
  lowest.mse.cv.min = min_MSE
  
  # Per task
  Coef_cv.glmnet <- list(features.glmnet.1se, features.glmnet.min) ; names(Coef_cv.glmnet) <- c("1se.mse","min.mse")
  MSE_cv.glmnet <- list(lowest.mse.cv.1se, lowest.mse.cv.min) ; names(MSE_cv.glmnet) <- c("1se.mse","min.mse")
  Hyperparameters_cv.glmnet <-  list(parameters.glmnet.1se, parameters.glmnet.min) ; names(Hyperparameters_cv.glmnet) <- c("1se.mse","min.mse")
  
  ## Stability Selection (it may take time if the number of replications is chosen very large and the number of
  ## core is not chosen high enough but it depends on your computer)
  p <- ncol(views_source.Z)
  q <- ncol(drug_source.Z)
  n <- nrow(drug_source.Z)
  
  stabsel.glmnet <- function(i) {
    b_sort <- sort(sample(1:(n * q), floor((n * q) / 2)))
    resultat_glmnet <- glmnet(views_source.Z[b_sort, ],
                              drug_source.Z[b_sort],
                              family = "gaussian",
                              alpha = Hyperparameters_cv.glmnet[["min.mse"]][["alpha"]],
                              lambda = Hyperparameters_cv.glmnet[["min.mse"]][["lambda"]]
    )
    ind_glmnet <- which(as.matrix(resultat_glmnet$beta) != 0)
    return(tabulate(ind_glmnet, (p * q)))
  }
  X <- views_source.Z
  Y <- drug_source.Z
  nb_repli = 1000
  nb.cores = 4
  
  res.cum <- Reduce("+", mclapply(1:nb_repli, stabsel.glmnet, mc.cores = nb.cores))

  freq <- res.cum / nb_repli
  if (is.null(colnames(Y))) {
    colnames(Y) <- 1:ncol(Y)
  }
  if (is.null(colnames(X))) {
    colnames(X) <- 1:ncol(X)
  }
  Freqs <- cbind(rep(colnames(Y), each = p), rep(colnames(X), q), as.data.frame(freq))
  names(Freqs) <- c("Names_of_Y", "Names_of_X", "frequency")
  
  return(list(cv.glmnet.data = data_task_cv.glmnet,
         cv.glmnet.hyperparameters = Hyperparameters_cv.glmnet,
         cv.glmnet.features = Coef_cv.glmnet,
         cv.glmnet.mse = MSE_cv.glmnet,
         cv.glmnet.freq = Freqs[2:3]))
}