#########################################################################################################
# Script to find the subnetwork of Cytokines pairs within the network of Ligand-Receptor pairs
#########################################################################################################

# ****************
# working directory
setwd("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data")

# ***************
# iTALK database
load("iTALK_database.rda")
cytokine.iTALK <- subset(database, Classification == "cytokine")
CYTOKINE.pairs_iTALK <- cytokine.iTALK$Pair.Name
CYTOKINES_iTALK <- unique(c(sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), head, 1), sapply(strsplit(CYTOKINE.pairs_iTALK, split = "_"), tail, 1)))
# write.csv(CYTOKINES_iTALK, file = "unique_cytokines_iTALK_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_iTALK.corrected <- read.csv(file = "hgnc-symbol-check.csv", header = TRUE, row.names = NULL, skip = 1)
# Matching with Ramilowski L-R pairs
CYTOKINE.pairs_iTALK <- LR.pairs.Ramilowski[na.exclude(which(LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_iTALK.corrected$Approved.symbol)),
                                            "Pair.Name", drop = FALSE]
CYTOKINE.pairs_iTALK <- unique(CYTOKINE.pairs_iTALK$Pair.Name)

# ***************
# ImmuneXpresso database
database.ImmuneXpresso <- read.csv("ImmuneXpressoResults.csv", header = TRUE, sep = ',')
load("LR.pairs.Ramilowski.RData")
CYTOKINES_ImmuneXpresso <- unique(database.ImmuneXpresso$Label)
# write.csv(CYTOKINES_ImmuneXpresso, file = "unique_cytokines_ImmuneXpresso_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_ImmuneXpresso.corrected <- read.csv(file = "hgnc-symbol-check-2.csv", header = TRUE, row.names = NULL, skip = 1)

# Matching with Ramilowski L-R pairs
CYTOKINE.pairs_ImmuneXpresso <- LR.pairs.Ramilowski[na.exclude(which(LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_ImmuneXpresso.corrected$Approved.symbol)), 
                                                    "Pair.Name", drop = FALSE]
CYTOKINE.pairs_ImmuneXpresso <- CYTOKINE.pairs_ImmuneXpresso$Pair.Name

# ***************
# CellPhoneDB database
CellPhoneDB.database <- read.csv("cellphonedb_protein_curated.csv", header = TRUE, sep = ",")
keep <- unique(grep("Cytokine", CellPhoneDB.database$secreted_desc, fixed = FALSE))

CYTOKINES_CellPhoneDB <- CellPhoneDB.database[keep,]
CYTOKINES_CellPhoneDB <- sapply(strsplit(CYTOKINES_CellPhoneDB$protein_name, split = "_", fixed = T), head, 1)
# write.csv(CYTOKINES_CellPhoneDB, file = "unique_cytokines_CellPhoneDB_input.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
CYTOKINES_CellPhoneDB.corrected <- read.csv(file = "hgnc-symbol-check-3.csv", header = TRUE, row.names = NULL, skip = 1)

# Matching with Ramilowski
CYTOKINE.pairs_CellPhoneDB <- LR.pairs.Ramilowski[na.exclude(which(LR.pairs.Ramilowski$Ligand.ApprovedSymbol %in% CYTOKINES_CellPhoneDB.corrected$Approved.symbol)),
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
# save(CYTOKINE.pairs_subnetwork, file = "CYTOKINE_pairs_subnetwork.RData")



