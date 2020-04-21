###########################################################################################
# Script to create ObtainIPS function:
# Compute ImmunoPhenoScore from log2(TPM+1) values
# Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
# and expression values (i.e. log2(TPM+1)) for each sample in the other columns
###########################################################################################

ObtainIPS <- function(gene_expr,....){

  # ************
  # Literature info: read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  IPSG_read <- read.table("./data/raw_data/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)

  # function: calculate Immunophenoscore
  ipsmap <- function (x) {
    if (x<=0) {
      ips<-0
    } else {
      if (x>=3) {
        ips<-10
      } else {
        ips<-round(x*10/3, digits=0)
      }
    }
    return(ips)
  }
  
  # ****************
  # calculation
  gene_expression <- gene_expr
  sample_names <- colnames(gene_expression)
  
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  IPSG <- IPSG_read
  unique_ips_genes <- as.vector(unique(IPSG$NAME))
  
  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL
  
  # Gene names in expression file
  GVEC <- row.names(gene_expression)
  # Genes names in IPS genes file
  VEC <- as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind <- which(is.na(match(VEC,GVEC)))
  # List genes missing or differently named
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }
  # for (x in 1:length(ind)) {
  #   print(IPSG[ind,])
  # }
  
  for (i in 1:length(sample_names)) {
    GE <- gene_expression[[i]]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    for (gen in unique_ips_genes) {
      MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k <- k+1
    }
    WG <- MIG*WEIGHT
    MHC[i] <- mean(WG[1:10])
    CP[i] <- mean(WG[11:20])
    EC[i] <- mean(WG[21:24])
    SC[i] <- mean(WG[25:26])
    AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }
  names(IPS) <- sample_names
  cat("IPS score computed","\n")
  return(data.frame(IPS = IPS, check.names = FALSE))
}