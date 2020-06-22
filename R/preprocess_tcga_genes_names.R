preprocess_tcga_genes_names <- function(gene.expr){
  
  # ****************
  # packages
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(AnnotationDbi))
  
  # Rows with non-valid HGNC symbols were removed.
  if (any(grep("\\?",rownames(gene.expr)))) {
    HGNC_symbol <- sapply(strsplit(rownames(gene.expr), split = "|", fixed = T), head, 1)
    gene.expr <-  gene.expr[-which(HGNC_symbol == "?"),]
    HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol == "?")]
  }
  if (any(grep("\\|",rownames(gene.expr)))) {
    # Keep entrez id for each gene
    entrez_id <- sapply(strsplit(rownames(gene.expr), split = "|", fixed = T),tail, 1)
    rownames(gene.expr) <- sapply(strsplit(rownames(gene.expr), split = "|", fixed = T), head, 1)
    
    # From Entrez id to symbols
    Symbols_from_id <- AnnotationDbi::select(org.Hs.eg.db, keys=entrez_id, columns='SYMBOL', keytype='ENTREZID')

    # TCGA hugo symbols
    pos_1 <- match(Symbols_from_id[!is.na(Symbols_from_id$SYMBOL),"ENTREZID"], entrez_id)
    rownames(gene.expr)[pos_1] <- Symbols_from_id[!is.na(Symbols_from_id$SYMBOL),"SYMBOL"]
    
    pos_2 <- match(Symbols_from_id[is.na(Symbols_from_id$SYMBOL),"ENTREZID"], entrez_id)
    # write.csv(rownames(gene.expr)[pos_2],
    #           file = paste0("../data/multi_symbol_checker/gene_names_to_be_checked_missing_symbol_from_id.csv"),
    #           quote = FALSE, row.names = FALSE, col.names = FALSE)
    gene.names_corrected <- read.csv(file = "../data/multi_symbol_checker/hgnc-symbol-check_missing_symbol_from_id.csv",
                                     header = TRUE, row.names = NULL, skip = 1)
    sub_gene.names_corrected <- subset(gene.names_corrected, Match.type %in% "Previous symbol")
    pos_3 <- match(sub_gene.names_corrected$Input, sapply(strsplit(rownames(gene.expr), split = "|", fixed = T), head, 1))
    rownames(gene.expr)[pos_3] <- sub_gene.names_corrected$Approved.symbol
  }
  
  # Rows corresponding to the same HGNC symbol were averaged.
  gene.symbol <- rownames(gene.expr)
  if(anyDuplicated(gene.symbol) != 0){
    idx <- which(duplicated(gene.symbol) == TRUE)
    dup_genes <- gene.symbol[idx]
    for (ii in dup_genes){
      gene.expr[which(gene.symbol %in% ii)[1],] <- colMeans(gene.expr[which(gene.symbol %in% ii),])
      gene.expr <- gene.expr[-which(gene.symbol %in% ii)[2],]
      gene.symbol <- gene.symbol[-which(gene.symbol %in% ii)[2]]
    }
  }
  rownames(gene.expr) <- gene.symbol
  
  return(gene.expr)
}