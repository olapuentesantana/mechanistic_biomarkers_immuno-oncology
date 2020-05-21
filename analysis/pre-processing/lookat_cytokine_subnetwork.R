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
Intermediate.table_Oscar <- Intermediate.table[which(Intermediate.table$Expr_Ligand >= 10 & Intermediate.table$Expr_Receptor >= 10),]
Intermediate.table_Maisa <- read.csv("Ligand_Receptors_Rdata/intermediate_table_Maisa.csv", header = TRUE)

# Expressed molecules per cell type: Pairs that ocurred most often (576 times) agreed with Maisa's report
# tmp <- subset(Intermediate.table, Ligand == "VIM" & Receptor == "CD44")
# unique(tmp$Ligand)
# unique(tmp$Receptor)

# Total number of unique pairs
LR.pairs_network_Oscar <- unique(paste0(Intermediate.table_Oscar$Ligand,"_",Intermediate.table_Oscar$Receptor))
LR.pairs_network_Maisa <- unique(paste0(Intermediate.table_Maisa$ligands,"_",Intermediate.table_Maisa$receptors))

# Need to redefine some gene names (that's the reason for the different unique L-R pairs)
# My simbols my be unapdated
setdiff(LR.pairs_network_Maisa, LR.pairs_network_Oscar)
setdiff(LR.pairs_network_Oscar, LR.pairs_network_Maisa)
genes_that_need_updated <- unique(unlist(strsplit(setdiff(LR.pairs_network_Oscar, LR.pairs_network_Maisa), split = "_")))
write.csv(genes_that_need_updated, file = "Ligand_Receptors_Rdata/genes_from_intermediate_table_to_update.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_that_need_updated_corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check-4.csv", header = TRUE, row.names = NULL, skip = 1)

for (X in 1:nrow(genes_that_need_updated_corrected)){
  
  LR.pairs_network_Oscar <- gsub(genes_that_need_updated_corrected$Input[X],genes_that_need_updated_corrected$Approved.symbol[X],
                                 LR.pairs_network_Oscar, fixed = F)
  
}
setdiff(LR.pairs_network_Maisa, LR.pairs_network_Oscar)
setdiff(LR.pairs_network_Maisa, LR.pairs_network_Oscar)

# Create filtered Ramilowsky dataset in order to obtain cytokines subnetwork only for the TME.
# Use LR.pairs.Ramilowski.RData from Maisa (it containts updated hugo symbols)
load("Ligand_Receptors_Rdata/LR.pairs.Ramilowski.RData")
filtered.LR.pairs.Ramilowski <- LR.pairs.Ramilowski[which(LR.pairs.Ramilowski$Pair.Name %in% LR.pairs_network_Oscar),c(1,2,4)]

# Cytokine Pairs for the TME
#############################################################################################################

# ***************
# iTALK database
load("Ligand_Receptors_Rdata/iTALK_database.rda")
cytokine.iTALK <- subset(database, Classification == "cytokine")
CYTOKINE.pairs_iTALK <- cytokine.iTALK$Pair.Name
CYTOKINES_iTALK <- unique(c(sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), head, 1), sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), tail, 1)))

# Obtaining cytokine pairs from iTALK database
keep_pairs <- intersect(which(cytokine.iTALK$Ligand.ApprovedSymbol %in% CYTOKINES_iTALK), 
                        which(cytokine.iTALK$Receptor.ApprovedSymbol %in% CYTOKINES_iTALK))

CYTOKINE.pairs_iTALK <- cytokine.iTALK[na.exclude(keep_pairs),"Pair.Name", drop = FALSE]
CYTOKINE.pairs_iTALK <- unique(CYTOKINE.pairs_iTALK$Pair.Name)

CYTOKINES_iTALK <- unique(c(sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_", fixed = TRUE), head, 1),
                            sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_", fixed = TRUE), tail, 1)))

# Take corrected HGNC symbols:
#write.csv(CYTOKINES_iTALK, file = "Ligand_Receptors_Rdata/unique_cytokines_iTALK_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_iTALK.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check.csv", header = TRUE, row.names = NULL, skip = 1)
CYTOKINES_iTALK.corrected <- subset(CYTOKINES_iTALK.corrected, Match.type %in% c("Approved symbol", "Previous symbol", "Alias symbol"))
CYTOKINES_iTALK.corrected <- CYTOKINES_iTALK.corrected[!duplicated(CYTOKINES_iTALK.corrected$Input),]

# Obtain iTALK cytokine pairs:
for (X in 1:nrow(CYTOKINES_iTALK.corrected)){
  
  CYTOKINE.pairs_iTALK <- gsub(CYTOKINES_iTALK.corrected$Input[X],CYTOKINES_iTALK.corrected$Approved.symbol[X],
                                 CYTOKINE.pairs_iTALK, fixed = F)
  
}

# ***************
# ImmuneXpresso database: 
# file to be downloaded: interactions betweeen selected cells and cytokines under certain context
database.ImmuneXpresso <- read.csv("Ligand_Receptors_Rdata/ImmuneXpressoResults.csv", header = TRUE, sep = ',')
CYTOKINES_ImmuneXpresso <- unique(database.ImmuneXpresso$Cytokine.Ontology.Label)

# Take corrected HGNC symbols:
#write.csv(CYTOKINES_ImmuneXpresso, file = "Ligand_Receptors_Rdata/unique_cytokines_ImmuneXpresso_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_ImmuneXpresso.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check-2.csv", header = TRUE, row.names = NULL, skip = 1)
CYTOKINES_ImmuneXpresso.corrected <- subset(CYTOKINES_ImmuneXpresso.corrected, Match.type %in% c("Approved symbol", "Previous symbol", "Alias symbol"))
CYTOKINES_ImmuneXpresso.corrected <- CYTOKINES_ImmuneXpresso.corrected[!duplicated(CYTOKINES_ImmuneXpresso.corrected$Input),]

# Obtain ImmuneExpresso cytokine pairs:
# Matching with Ramilowski L-R pairs (only ligands)
keep_pairs <- which(LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_ImmuneXpresso.corrected$Approved.symbol)
                        # which(LR.pairs.Ramilowski$Receptor.ApprovedSymbol %in% CYTOKINES_ImmuneXpresso.corrected$Approved.symbol))

CYTOKINE.pairs_ImmuneXpresso <- LR.pairs.Ramilowski[na.exclude(keep_pairs),"Pair.Name", drop = FALSE]
CYTOKINE.pairs_ImmuneXpresso <- unique(CYTOKINE.pairs_ImmuneXpresso$Pair.Name)

# ***************
# CellPhoneDB database
CellPhoneDB.database_Oscar <- read.csv("Ligand_Receptors_Rdata/cellphonedb_protein_curated.csv", header = TRUE, sep = ",")
# CellPhoneDB.database_Maisa <- read.csv("Ligand_Receptors_Rdata/Maisa_download/protein_complex_cellphonedb.csv", header = TRUE, sep = ",")
keep_Oscar <- unique(c(grep("Cytokine", CellPhoneDB.database_Oscar$secreted_desc, fixed = FALSE), grep("Cytokine", CellPhoneDB.database_Oscar$receptor_desc, fixed = FALSE)))
keep_Maisa <- unique(c(grep("Cytokine", CellPhoneDB.database_Maisa$secreted_desc, fixed = FALSE)))

CYTOKINES_CellPhoneDB <- CellPhoneDB.database_Oscar[keep_Oscar,]
CYTOKINES_CellPhoneDB <- sapply(strsplit(CYTOKINES_CellPhoneDB$protein_name, split = "_", fixed = T), head, 1)

# Take corrected HGNC symbols:
# write.csv(CYTOKINES_CellPhoneDB, file = "Ligand_Receptors_Rdata/unique_cytokines_CellPhoneDB_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_CellPhoneDB.corrected <- read.csv(file = "Ligand_Receptors_Rdata/hgnc-symbol-check-3.csv", header = TRUE, row.names = NULL, skip = 1)
CYTOKINES_CellPhoneDB.corrected <- subset(CYTOKINES_CellPhoneDB.corrected, Match.type %in% c("Approved symbol", "Previous symbol", "Alias symbol"))
CYTOKINES_CellPhoneDB.corrected <- CYTOKINES_CellPhoneDB.corrected[!duplicated(CYTOKINES_CellPhoneDB.corrected$Input),]

# Obtain CellPhoneDB cytokine pairs:
# Matching with Ramilowski
keep_pairs <- unique(c(which(LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_CellPhoneDB.corrected$Approved.symbol), 
                       which(LR.pairs.Ramilowski$Receptor.ApprovedSymbol %in% CYTOKINES_CellPhoneDB.corrected$Approved.symbol)))

CYTOKINE.pairs_CellPhoneDB <- LR.pairs.Ramilowski[na.exclude(keep_pairs),"Pair.Name", drop = FALSE]
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
# filtered for L-R pairs only in the TME:
CYTOKINE.pairs_subnetwork <- CYTOKINE.pairs_subnetwork[which(CYTOKINE.pairs_subnetwork %in% LR.pairs_network_Oscar)]
#save(CYTOKINE.pairs_subnetwork, file = "Ligand_Receptors_Rdata/CYTOKINE_pairs_subnetwork_only_TME.RData")
# save(LR.pairs_network, file = "Ligand_Receptors_Rdata/LR_pairs_network_only_TME.RData")



