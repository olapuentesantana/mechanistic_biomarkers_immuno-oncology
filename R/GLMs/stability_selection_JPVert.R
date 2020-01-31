stability_selection <- function(x,y,nbootstrap=100,nsteps=20,alpha=0.2,plotme=FALSE)
{
  # Stability selection in the spirit of Meinshausen&Buhlman
  # JP Vert, 14/9/2010
  
  # x is the n*p design matrix, y the n*1 variable to predict
  # x should be normalized to unit variance per column before calling this function (if you want to)
  # the result is a score (length p) for each feature, the probability that each feature is selected during the first nsteps steps of the Lasso path when half of the samples are used and the features are reweigthed by a random weight uniformaly sampled in [alpha,1]. This probability is estimated by nbootstrap bootstrap samples
  require(lars)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  freq <- matrix(0,nsteps+1,p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    xs <- t(t(x)*runif(p,alpha,1))
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    # run the randomized lasso on each sample and check which variables are selected
    r <- lars(xs[i1,],y[i1],max.steps=nsteps,normalize=FALSE)
    freq <- freq + abs(sign(coef.lars(r)))
    r <- lars(xs[i2,],y[i2],max.steps=nsteps,normalize=FALSE)
    freq <- freq + abs(sign(coef.lars(r)))
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="LARS iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
}