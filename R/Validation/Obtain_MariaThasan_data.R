BiocManager::install(c("biomaRt",
           "circlize",
           "ComplexHeatmap",
           "corrplot",
           "DESeq2",
           "dplyr",
           "DT",
           "edgeR",
           "ggplot2",
           "limma",
           "lsmeans",
           "reshape2",
           "spatstat",
           "survival",
           "plyr"))

install.packages("Downloads/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)
library(IMvigor210CoreBiologies)

# Preprocessed data
# Transcriptome wide gene expression data
# To load a CountDataSet object called ‘cds’, type:
data(cds)
# This CountDataSet object contains raw counts for all genes as well as basic feature and all sample annotations
# reported in the manuscript. These can be accessed like this:
head(counts(cds))
head(fData(cds))
head(pData(cds))
