elastic_net_test <- function(x.test, coef){
  
  Intercept = as.matrix(coef)[1]
  coef = coef[-1, ,drop = FALSE]
  # match features properly
  pos = na.omit(match(colnames(x.test), rownames(coef)))
  coef = coef[pos,, drop = FALSE]
  
  if  (length(coef) > 1){
    Slope <-  coef
    rownames(Slope) <- gsub(" ","",rownames(Slope)) # In case of Dorothea is needed, do not affect the other features
    fit.pred <- t(as.numeric(Intercept) + t(as.numeric(Slope)) %*% t(as.matrix(x.test[,rownames(Slope)])))
    
  }else{
    Slope <- 0
    fit.pred <- as.numeric(Intercept) 
  }
  return(fit.pred)
}