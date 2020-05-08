# #########################################################################################################
# Script to understand how L21 regularization works
# #########################################################################################################

# setwd("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop")
# X <- read.csv("R/L21_regularization_EN/python_tmp/pathways_immunecells_SKCM.csv", sep = ",",header = T)
# Y <- read.csv("R/L21_regularization_EN/python_tmp/drugs_SKCM.csv", sep = ",",header = T)
load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new_remove_all_genes/DataViews_no_filter_GBM.RData")
X <- DataViews.no_filter$pathways
load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new/ImmuneResponse_no_filter_GBM_matrix_format.RData")
Y <- ImmuneResponse.no_filter

library(RMTL)
mas.mea.learning.X = colMeans(as.matrix(X), na.rm = T)
mas.std.learning.X = colSds(as.matrix(X), na.rm = T)

X <- sweep(X, 2,   mas.mea.learning.X, FUN = "-")
X <- sweep(X, 2,   mas.std.learning.X, FUN = "/")

mas.mea.learning.Y = colMeans(as.matrix(Y), na.rm = T)
mas.std.learning.Y = colSds(as.matrix(Y), na.rm = T)

Y <- sweep(Y, 2,   mas.mea.learning.Y, FUN = "-")
Y <- sweep(Y, 2,   mas.std.learning.Y, FUN = "/")

colnames(ImmuneResponse.no_filter)
colnames(X)
data <- list()

data$X[[1]] <- as.matrix(X)
# data$X[[2]] <- as.matrix(X)
# data$X[[3]] <- as.matrix(X)
# data$X[[4]] <- as.matrix(X)
# data$X[[5]] <- as.matrix(X)
# data$X[[6]] <- as.matrix(X)
# data$X[[7]] <- as.matrix(X)
# data$X[[8]] <- as.matrix(X)

data$Y[[1]] <- as.matrix(Y[,"IFny"])
# data$Y[[2]] <- as.matrix(Y[,"ExpandedImmune"])
# data$Y[[3]] <- as.matrix(Y[,"T_cell_inflamed"])
# data$Y[[4]] <- as.matrix(Y[,9])
# data$Y[[5]] <- as.matrix(Y[,10])
# data$Y[[6]] <- as.matrix(Y[,11])
# data$Y[[7]] <- as.matrix(Y[,5])
# data$Y[[8]] <- as.matrix(Y[,2])

alpha <- c(1.24046956e-01, 1.15686607e-01, 1.07889717e-01, 1.00618311e-01,
           9.38369734e-02, 8.75126755e-02, 8.16146140e-02, 7.61140622e-02,
           7.09842292e-02, 6.62001297e-02, 6.17384625e-02, 5.75774968e-02,
           5.36969662e-02, 5.00779704e-02, 4.67028828e-02, 4.35552648e-02,
           4.06197857e-02, 3.78821481e-02, 3.53290181e-02, 3.29479605e-02,
           3.07273783e-02, 2.86564558e-02, 2.67251066e-02, 2.49239239e-02,
           2.32441348e-02, 2.16775579e-02, 2.02165631e-02, 1.88540344e-02,
           1.75833356e-02, 1.63982776e-02, 1.52930886e-02, 1.42623856e-02,
           1.33011485e-02, 1.24046956e-02, 1.15686607e-02, 1.07889717e-02,
           1.00618311e-02, 9.38369734e-03, 8.75126755e-03, 8.16146140e-03,
           7.61140622e-03, 7.09842292e-03, 6.62001297e-03, 6.17384625e-03,
           5.75774968e-03, 5.36969662e-03, 5.00779704e-03, 4.67028828e-03,
           4.35552648e-03, 4.06197857e-03, 3.78821481e-03, 3.53290181e-03,
           3.29479605e-03, 3.07273783e-03, 2.86564558e-03, 2.67251066e-03,
           2.49239239e-03, 2.32441348e-03, 2.16775579e-03, 2.02165631e-03,
           1.88540344e-03, 1.75833356e-03, 1.63982776e-03, 1.52930886e-03,
           1.42623856e-03, 1.33011485e-03, 1.24046956e-03, 1.15686607e-03,
           1.07889717e-03, 1.00618311e-03, 9.38369734e-04, 8.75126755e-04,
           8.16146140e-04, 7.61140622e-04, 7.09842292e-04, 6.62001297e-04,
           6.17384625e-04, 5.75774968e-04, 5.36969662e-04, 5.00779704e-04,
           4.67028828e-04, 4.35552648e-04, 4.06197857e-04, 3.78821481e-04,
           3.53290181e-04, 3.29479605e-04, 3.07273783e-04, 2.86564558e-04,
           2.67251066e-04, 2.49239239e-04, 2.32441348e-04, 2.16775579e-04,
           2.02165631e-04, 1.88540344e-04, 1.75833356e-04, 1.63982776e-04,
           1.52930886e-04, 1.42623856e-04, 1.33011485e-04, 1.24046956e-04)
datar <- Create_simulated_data(Regularization="L21", type="Regression")
# L21 built-in cross validation
assign("fit",lapply(alpha_range, function(X){
  
  # Same foldid for each alpha, but different foldid per iteration.
  set.seed(3000 + 10 * iteration)
  
  # Applying internal cross-validation:
  cvfit <- cvMTL(data$X , data$Y, type="Regression", Regularization="L21", Lam1_seq=alpha, Lam2=0,
                 nfolds = 5, parallel = TRUE, ncores = 2, opts=list(init=0,  tol=10^-3, maxIter=1500))
  
  # Train
  model <- MTL(data$X, data$Y, type="Regression", Regularization="L21",
               Lam1=cvfit$Lam1.min, Lam2 =cvfit$Lam2 ,
               opts=list(init=0,  tol=10^-3, maxIter=1500))
  
  return(list(cv = list(Lam1_seq = cvfit$Lam1_seq, cvm = cvfit$cvm), model = model))
  
}))

names(fit) <- paste0(alpha_range,"-")






