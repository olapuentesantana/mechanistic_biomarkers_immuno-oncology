preprocess_gene_names_to_match_LRpairs_genes <- function(RNA.tpm){
  
  # ****************
  # data 
  load("../data/all_genes_ICB_proxies.RData")
  load("../data/gene_names_updated_LRpairs.RData")
  all_LRpairs <- unique(gene.names_LRpairs$Input)
  
  # redefine gene names to match transcripts in TCGA data
  genes_kept <- intersect(rownames(RNA.tpm), all_LRpairs)
  genes_left <- setdiff(all_LRpairs, rownames(RNA.tpm))
  
  sub_genes_left <- gene.names_LRpairs[gene.names_LRpairs$Input %in% setdiff(genes_left, all_genes_to_remove),]
  sub_genes_left <- subset(sub_genes_left, Match.type == "Previous symbol")
  RNA.tpm_post <- RNA.tpm
  
  pick_1 <- na.omit(match(sub_genes_left$Approved.symbol, rownames(RNA.tpm_post)))
  pick_2 <- na.omit(match(rownames(RNA.tpm_post), sub_genes_left$Approved.symbol))
  rownames(RNA.tpm_post)[pick_1] <- sub_genes_left$Input[pick_2]
  
  return(RNA.tpm_post)
}
