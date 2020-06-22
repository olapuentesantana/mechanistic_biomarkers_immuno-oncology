preprocess_gene_names_to_match_pathways_genes <- function(RNA.raw_counts){

  # ****************
  # data 
  load("../data/all_genes_ICB_proxies.RData")
  # Tf regulons target genes updates (using multi-symbol checker)
  load("../data/gene_names_updated_pathways.RData")
  pathway_responsive_genes <- unique(gene.names_pathways$Input)
  
  # redefine gene names to match transcripts in TCGA data
  genes_kept <- intersect(rownames(RNA.raw_counts), pathway_responsive_genes)
  genes_left <- setdiff(pathway_responsive_genes, rownames(RNA.raw_counts))
  
  sub_genes_left <- gene.names_pathways[gene.names_pathways$Input %in% setdiff(genes_left, all_genes_to_remove),]
  sub_genes_left <- subset(sub_genes_left, Match.type == "Previous symbol")
  RNA.raw_counts_post <- RNA.raw_counts
  
  pick_1 <- na.omit(match(sub_genes_left$Approved.symbol, rownames(RNA.raw_counts_post)))
  pick_2 <- na.omit(match(rownames(RNA.raw_counts_post), sub_genes_left$Approved.symbol))
  rownames(RNA.raw_counts_post)[pick_1] <- sub_genes_left$Input[pick_2]
  
  return(RNA.raw_counts_post)
}

  