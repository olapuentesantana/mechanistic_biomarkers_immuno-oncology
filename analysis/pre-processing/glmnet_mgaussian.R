# #########################################################################################################
# Script to understand how glmnet multi-gaussian works
# #########################################################################################################

# setwd("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop")
# X <- read.csv("R/L21_regularization_EN/python_tmp/pathways_immunecells_SKCM.csv", sep = ",",header = T)
# Y <- read.csv("R/L21_regularization_EN/python_tmp/drugs_SKCM.csv", sep = ",",header = T)
load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new_remove_all_genes/DataViews_no_filter_GBM.RData")
X <- DataViews.no_filter$pathways
load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new/ImmuneResponse_no_filter_GBM_matrix_format.RData")
Y <- ImmuneResponse.no_filter

library(glmnet)

mas.mea.learning.X = colMeans(as.matrix(X), na.rm = T)
mas.std.learning.X = colSds(as.matrix(X), na.rm = T)

X <- sweep(X, 2,   mas.mea.learning.X, FUN = "-")
X <- sweep(X, 2,   mas.std.learning.X, FUN = "/")

mas.mea.learning.Y = colMeans(as.matrix(Y), na.rm = T)
mas.std.learning.Y = colSds(as.matrix(Y), na.rm = T)

Y <- sweep(Y, 2,   mas.mea.learning.Y, FUN = "-")
Y <- sweep(Y, 2,   mas.std.learning.Y, FUN = "/")

# Internal cross validation: number of folds by default is 10
# Same foldid for each alpha, but different foldid per iteration.
iteration <- 1
set.seed(3000 + 10 * iteration)

# Applying internal cross-validation:
fit <- cv.glmnet(x = as.matrix(X), 
                 y = as.matrix(Y), 
                 family = "mgaussian",
                 alpha = 1, 
                 type.measure = "mse", 
                 standardize = FALSE, 
                 intercept = TRUE, 
                 keep = TRUE,
                 nfolds = 5,
                 parallel = T)


coef(fit)

a <- cor(DataViews.no_filter$pathways)

