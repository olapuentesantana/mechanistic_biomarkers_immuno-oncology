#########################################################################################################
# Script to find the subnetwork of Cytokines pairs within the network of Ligand-Receptor pairs
#########################################################################################################

# ****************
# working directory
setwd("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data")

# Ligand - Receptor Pairs for the TME
#############################################################################################################
# ************
# Input data:
# Average CAGE expression data for human protein-coding genes in the 144 human primary cells
Expr.Lig_Rec <- read.delim(file = "Ligand_Receptors_Rdata/ExpressionLigRec.txt")
# Average CAGE expression data for all ligands and receptors in the 144 human primary cells
Expr_genes <- read.delim(file = "Ligand_Receptors_Rdata/ExpressionGenes.txt")
# Curated ligand-receptor pairs in Ramilowski database (as of April 2015)
Pairs.L.R <- read.delim(file = "Ligand_Receptors_Rdata/PairsLigRec.txt")

# load("data/Ligand_Receptors_Rdata/LR.pairs.Ramilowski.RData")

# ************
# Processing

# Filter by cell types commonly present in the TME (24 cell types) [source: Maisa's report]
TME.cell_types <- c("Mature.Adipocyte","Adipocyte.Omental","CD19..B.cells", "CD4..T.cells",
                    "CD4.CD25.CD45RA..naive.regulatory.T.cells","CD4.CD25.CD45RA..memory.regulatory.T.cells",       
                    "CD4.CD25.CD45RA..naive.conventional.T.cells", "CD4.CD25.CD45RA..memory.conventional.T.cells",
                    "CD8..T.cells", "Dendritic.Monocyte.Immature.derived", "Dendritic.Plasmacytoid", 
                    "Endothelial.Lymphatic","Endothelial.Microvascular", "Fibroblast.Lymphatic", "Fibroblast.Skin.Normal",
                    "Macrophage.Monocyte.derived", "Mast.cells", "Mast.cells.stimulated","CD14..Monocytes",
                    "CD14.CD16..Monocytes", "CD14.CD16..Monocytes.1", "CD14.CD16..Monocytes.2", "NK.cells", "Neutrophils")

Expr.Lig_Rec <- Expr.Lig_Rec[,c(1,2,3,4,match(TME.cell_types, colnames(Expr.Lig_Rec)))]
# We selected only literature supported evidence pairs to build this table
Pairs.L.R <- subset(Pairs.L.R, Pair.Evidence == "literature supported")
Expr.Lig_Rec <- subset(Expr.Lig_Rec, Pair.Evidence == "literature supported")

# Exclude those ligands and receptors that were expressed by a cell type 
# but not paired to another ligand or receptor in the cell network.
Expr.L <- subset(Expr.Lig_Rec, Type == "ligand")
Expr.R <- subset(Expr.Lig_Rec, Type == "receptor")

# Built intermediate table: each sender cell type and its ligand,
# each receiver cell type and its receptor, both with their corresponding expression.
rownames(Expr.L) <- Expr.L$ApprovedSymbol
rownames(Expr.R) <- Expr.R$ApprovedSymbol
Expr.L <- Expr.L[,-c(1,2,3,4)]
Expr.R <- Expr.R[,-c(1,2,3,4)]

Expr.L <- Expr.L[!is.na(rowSums(Expr.L)),]
Expr.R <- Expr.R[!is.na(rowSums(Expr.R)),]

keep_pairs <- intersect(which(Pairs.L.R$Ligand.ApprovedSymbol %in% rownames(Expr.L)), which(Pairs.L.R$Receptor.ApprovedSymbol %in% rownames(Expr.R)))


Receptors_for_Ligands <- Pairs.L.R[match(keep_ligands, Pairs.L.R$Ligand.ApprovedSymbol), c(1,2,4)]
Ligands_for_Receptors <- Pairs.L.R[match(keep_receptors, Pairs.L.R$Receptor.ApprovedSymbol), c(1,2,4)]

LR.pairs <- Pairs.L.R[keep_pairs, c(1,2,4)]
LR.pairs <- LR.pairs[!duplicated(LR.pairs$Pair.Name),]

Intermediate.table <- do.call(rbind, lapply(colnames(Expr.L), function(cell_type){
  
  Intermediate.table <- do.call(rbind, lapply(1:nrow(LR.pairs), function(L_pos){
  
           L <- LR.pairs$Ligand.ApprovedSymbol[L_pos]
           R <- LR.pairs$Receptor.ApprovedSymbol[L_pos]

          tmp.table <- data.frame(Sender = cell_type, 
                                  Ligand = L, 
                                  Expr_Ligand = as.numeric(Expr.L[L,cell_type]),
                                  Receiver = colnames(Expr.R), 
                                  Receptor = R,
                                  Expr_Receptor= as.numeric(Expr.R[R,]))
          
    return(tmp.table)
  }))
  return(Intermediate.table)
}))


# The ligands and receptors were only taken into account if their
# expression was >= 10 TPM in the given cell
anyNA(Intermediate.table) # FALSE
Intermediate.table <- Intermediate.table[which(Intermediate.table$Expr_Ligand >= 10 & Intermediate.table$Expr_Receptor >= 10),]

# Expressed molecules per cell type: Pairs that ocurred most often (576 times) agreed with Maisa's report
# tmp <- subset(Intermediate.table, Ligand == "VIM" & Receptor == "CD44")
unique(tmp$Ligand)
unique(tmp$Receptor)

# Total number of unique pairs
LR.pairs_network <- unique(paste0(Intermediate.table$Ligand,"_",Intermediate.table$Receptor))

# Create filtered Ramilowsky dataset in order to obtain cytokines subnetwork only for the TME.
filtered.LR.pairs.Ramilowski <- LR.pairs[which(LR.pairs$Pair.Name %in% LR.pairs_network),]

# Cytokine Pairs for the TME
#############################################################################################################

# ***************
# iTALK database
load("Ligand_Receptors_Rdata/iTALK_database.rda")
cytokine.iTALK <- subset(database, Classification == "cytokine")
CYTOKINE.pairs_iTALK <- cytokine.iTALK$Pair.Name
CYTOKINES_iTALK <- unique(c(sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), head, 1), sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), tail, 1)))
# write.csv(CYTOKINES_iTALK, file = "unique_cytokines_iTALK_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_iTALK.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check.csv", header = TRUE, row.names = NULL, skip = 1)
# Matching with Ramilowski L-R pairs
CYTOKINE.pairs_iTALK <- filtered.LR.pairs.Ramilowski[na.exclude(which(filtered.LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_iTALK.corrected$Approved.symbol)),
                                            "Pair.Name", drop = FALSE]
CYTOKINE.pairs_iTALK <- unique(CYTOKINE.pairs_iTALK$Pair.Name)

# ***************
# ImmuneXpresso database
database.ImmuneXpresso <- read.csv("Ligand_Receptors_Rdata/ImmuneXpressoResults.csv", header = TRUE, sep = ',')
CYTOKINES_ImmuneXpresso <- unique(database.ImmuneXpresso$Label)
# write.csv(CYTOKINES_ImmuneXpresso, file = "unique_cytokines_ImmuneXpresso_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_ImmuneXpresso.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check-2.csv", header = TRUE, row.names = NULL, skip = 1)

# Matching with Ramilowski L-R pairs
CYTOKINE.pairs_ImmuneXpresso <- filtered.LR.pairs.Ramilowski[na.exclude(which(filtered.LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_ImmuneXpresso.corrected$Approved.symbol)), 
                                                    "Pair.Name", drop = FALSE]
CYTOKINE.pairs_ImmuneXpresso <- CYTOKINE.pairs_ImmuneXpresso$Pair.Name

# ***************
# CellPhoneDB database
CellPhoneDB.database <- read.csv("Ligand_Receptors_Rdata/cellphonedb_protein_curated.csv", header = TRUE, sep = ",")
keep <- unique(grep("Cytokine", CellPhoneDB.database$secreted_desc, fixed = FALSE))

CYTOKINES_CellPhoneDB <- CellPhoneDB.database[keep,]
CYTOKINES_CellPhoneDB <- sapply(strsplit(CYTOKINES_CellPhoneDB$protein_name, split = "_", fixed = T), head, 1)
# write.csv(CYTOKINES_CellPhoneDB, file = "unique_cytokines_CellPhoneDB_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_CellPhoneDB.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check-3.csv", header = TRUE, row.names = NULL, skip = 1)

# Matching with Ramilowski
CYTOKINE.pairs_CellPhoneDB <- filtered.LR.pairs.Ramilowski[na.exclude(which(filtered.LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_CellPhoneDB.corrected$Approved.symbol)),
                                                  "Pair.Name", drop = FALSE]

CYTOKINE.pairs_CellPhoneDB <- unique(CYTOKINE.pairs_CellPhoneDB$Pair.Name)

# ***************
# Intersecting cytokine -cytokine receptor pairs in the three databases
CYTOKINE.pairs_subnetwork <- Reduce(intersect, list(CYTOKINE.pairs_iTALK, CYTOKINE.pairs_ImmuneXpresso, CYTOKINE.pairs_CellPhoneDB))

# Intersecting cytokines in the three databases
CYTOKINES_subnetwork <- Reduce(intersect, list(CYTOKINES_iTALK.corrected$Approved.symbol,
                                               CYTOKINES_ImmuneXpresso.corrected$Approved.symbol,
                                               CYTOKINES_CellPhoneDB.corrected$Approved.symbol))

# ***************
# Comparing with Maisa externship intersection
we_have_more <- setdiff(CYTOKINE.pairs_subnetwork, colnames(DataViews.no_filter$CYTOKINEpairs))
setdiff(colnames(DataViews.no_filter$CYTOKINEpairs), CYTOKINE.pairs_subnetwork)

# ***************
# save list of cytokines pairs derived from intersection between the three databases
# save(CYTOKINE.pairs_subnetwork, file = "Ligand_Receptors_Rdata/CYTOKINE_pairs_subnetwork_only_TME.RData")
# save(LR.pairs_network, file = "Ligand_Receptors_Rdata/LR_pairs_network_only_TME.RData")



