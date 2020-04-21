L21_test <- function(x.test, coef.matrix){
  
  Intercept <- as.matrix(coef.matrix$min.mse)[1,]
  coef <- as.matrix(coef.matrix$min.mse)[-1, ,drop = FALSE]
  
  # Combine views
  x.test.combo <- do.call(cbind, lapply(1:length(x.test), function(x){tmp = x.test[[x]]}))
  
  # match features properly
  pos <- na.omit(match(colnames(x.test.combo), rownames(coef)))
  coef <- coef[pos,, drop = FALSE]
  
  if (length(coef) > 1){
    Slope <-  coef
    rownames(Slope) <- gsub(" ","",rownames(Slope)) # In case of Dorothea is needed, do not affect the other features
    fit.pred <- t(matrix(as.matrix(Intercept), nrow = ncol(Slope), ncol = nrow(x.test.combo))
                  + t(Slope) %*% t(as.matrix(x.test.combo[,rownames(Slope)])))
  }else{
    Slope <- 0
    fit.pred <- matrix(as.matrix(Intercept), nrow = ncol(Slope), ncol = nrow(x.test[[1]]))
  }
  return(fit.pred)
}