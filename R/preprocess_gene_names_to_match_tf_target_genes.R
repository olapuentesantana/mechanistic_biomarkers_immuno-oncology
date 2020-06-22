preprocess_gene_names_to_match_tf_target_genes <- function(RNA.tpm){
  
  # ****************
  # data 
  load("../data/all_genes_ICB_proxies.RData")
  # Pathway responsive genes updates (using multi-symbol checker)
  load("../data/gene_names_updated_tfs.RData")
  all_regulated_transcripts <- unique(gene.names_tfs$Input)
  
  # redefine gene names to match transcripts in TCGA data
  genes_kept <- intersect(rownames(RNA.tpm), all_regulated_transcripts)
  genes_left <- setdiff(all_regulated_transcripts, rownames(RNA.tpm))
  
  sub_genes_left <- gene.names_tfs[gene.names_tfs$Input %in% setdiff(genes_left, all_genes_to_remove),]
  sub_genes_left <- subset(sub_genes_left, Match.type == "Previous symbol")
  RNA.tpm_post <- RNA.tpm
  
  pick_1 <- na.omit(match(sub_genes_left$Approved.symbol, rownames(RNA.tpm_post)))
  pick_2 <- na.omit(match(rownames(RNA.tpm_post), sub_genes_left$Approved.symbol))
  rownames(RNA.tpm_post)[pick_1] <- sub_genes_left$Input[pick_2]
  
  return(RNA.tpm_post)
}

