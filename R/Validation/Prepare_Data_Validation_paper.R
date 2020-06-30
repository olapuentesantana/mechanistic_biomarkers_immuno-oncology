#########################################################################################################
# Script to generate DataViews from external datasets:

# External datasets --> 
  ## Auslander, Gide, Hugo, Kim, Liu, Riaz, Zhao
#########################################################################################################

# ****************
# working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/R")

# ****************
# scritps 
source("compute.pathways.scores.R")
source("compute.TF.activity.R")
source("compute.LR.pairs.R")
# Pre-processing tcga gene names to match pathway genes
source("preprocess_gene_names_to_match_pathway_genes.R")
# Pre-processing tcga gene names to match pathway genes
source("preprocess_gene_names_to_match_tf_target_genes.R")
# Pre-processing tcga gene names to match LR pairs genes
source("preprocess_gene_names_to_match_LRpairs_genes.R")

# ****************
# External datasets
Datasets.names <-  dir("../data/Validation/Francesca", full.names = F, recursive = F)
DataViews.test <- list()
# Labels.test <- list()
DataViews.test <- do.call(c, lapply(Datasets.names, function(dataset){
  
  # Samples information
  load(paste0("../data/Validation/Francesca/", dataset, "/", dataset, "_sampleinfo.rdata"))
  
  cat("\n",dataset,": \n")
  
  # Raw counts values
  gene.count <- read.table(file = paste0("../data/Validation/Francesca/", dataset, "/", dataset, "_gene_count.txt"),  
                           sep = "\t", header = TRUE, row.names = 1)
  
  # Tpm values
  gene.tpm <- read.table(file = paste0("../data/Validation/Francesca/", dataset, "/", dataset, "_gene_tpm.txt"),  
                         sep = "\t", header = TRUE, row.names = 1)
  
  # Cell fractions
  cell_fractions <- as.data.frame(t(read.table(file = paste0("../data/Validation/Francesca/", dataset, "/", dataset, "_cell_fractions.txt"),  
                               sep = "\t", header = TRUE, row.names = 1)))
  # # LRpairs
  # LRpairs <- read.table(file = paste0("../data/Validation/Francesca/", dataset, "/", tolower(dataset), ".txt"),  
  #                        sep = "\t", header = TRUE, row.names = 2)
  # # Remove three annotation columns
  # LRpairs <- as.data.frame(LRpairs[,-c(1,2,3)])
  # 
  # # CYTOKINEpairs
  # CYTOKINEpairs <- read.table(file = paste0("../data/Validation/Francesca/", dataset, "/", tolower(dataset), "_cytokine.txt"),  
  #                        sep = "\t", header = TRUE, row.names = 2)
  # # Remove three annotation columns
  # CYTOKINEpairs <- as.data.frame(CYTOKINEpairs[,-c(1,2,3)])
  
  #--------------------------------------------------------------------#
  # Filter samples according to clinical data (pre-treatment samples only) 
  #--------------------------------------------------------------------#
  if(dataset == "Riaz"){
    
    pre.sampleinfo <- subset(sampleinfo, Time == "Pre" & Response != "UNK")
    order.pre <- match(pre.sampleinfo$SRA, colnames(gene.count))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions
    cell_fractions.pre <- cell_fractions[,match(pre.sampleinfo$SRA, colnames(cell_fractions))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, match(pre.sampleinfo$SRA, colnames(LRpairs))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, match(pre.sampleinfo$SRA, colnames(CYTOKINEpairs))]
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, match(pre.sampleinfo$SRA, colnames(gene.tpm))]
    }
    
    label.pre <- data.frame(Sample = colnames(gene.tpm.pre), label = pre.sampleinfo$Response)
  }
  
  if(dataset == "Hugo"){
    
    samples.available <- intersect(sampleinfo$patient_id_s, colnames(cell_fractions))
    # Two samples for Pt27 in RNAseq data, but not in cell-fractions *
    # SRR3184292 is an early on-treatment sample (Pt16)
    pre.sampleinfo <- subset(sampleinfo, patient_id_s != "Pt16")
    
    order.pre <- match(rownames(pre.sampleinfo), colnames(gene.count))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions
    cell_fractions.pre <- cell_fractions[, match(pre.sampleinfo$patient_id_s, colnames(cell_fractions))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    colnames(cell_fractions.pre) <- rownames(pre.sampleinfo)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, match(rownames(pre.sampleinfo), colnames(LRpairs))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, match(rownames(pre.sampleinfo), colnames(CYTOKINEpairs))]
    
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, match(rownames(pre.sampleinfo), colnames(gene.tpm))]
    }
    label.pre <- data.frame(Sample = rownames(pre.sampleinfo), label = pre.sampleinfo$anti_pd_1_response_s)
  }
  
  if(dataset == "Auslander"){
    
    sampleinfo$Time <- as.numeric(as.vector(sampleinfo$Time))
    pre.sampleinfo <- subset(sampleinfo, Time <= 0 )
    
    order.pre <- match(pre.sampleinfo$SRA, colnames(gene.count))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions
    cell_fractions.pre <- cell_fractions[,match(pre.sampleinfo$SRA, colnames(cell_fractions))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, match(pre.sampleinfo$SRA, colnames(LRpairs))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, match(pre.sampleinfo$SRA, colnames(CYTOKINEpairs))]
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, match(pre.sampleinfo$SRA, colnames(gene.tpm))]
    }
    label.pre <- data.frame(Sample = pre.sampleinfo$SRA, label = pre.sampleinfo$Response)
  }
  
  if(dataset == "Gide"){
    
    clinical_PD1 <- read.csv(file = "../data/Validation/Francesca/Gide/1-s2.0-S1535610819300376-mmc2.csv", 
                             sep = ",", header = T, skip = 2)
    clinical_PD1$Patient.no. <- paste0("PD1_",clinical_PD1$Patient.no., "_", clinical_PD1$RNA.Sequencing)
    
    clinical_PD1_CTLA4 <- read.csv(file = "../data/Validation/Francesca/Gide/1-s2.0-S1535610819300376-mmc3.csv",
                                   sep = ",", header = T, skip = 2)
    
    clinical_PD1_CTLA4$RNA.Sequencing[38] <- "PRE" # mistake in clinical data
    clinical_PD1_CTLA4$Patient.no. <- paste0("ipiPD1_",clinical_PD1_CTLA4$Patient.no., "_", clinical_PD1_CTLA4$RNA.Sequencing)
    
    clinical_PD1.pre <- subset(clinical_PD1, RNA.Sequencing %in% c("PRE", "PRE ", "PRE and EDT"))
    clinical_PD1_CTLA4.pre <- subset(clinical_PD1_CTLA4, RNA.Sequencing %in% c("PRE", "PRE ", "PRE and EDT"))
    
    tmp.response.pre <- data.frame(Alias =  c(clinical_PD1.pre$Patient.no.,
                                              clinical_PD1_CTLA4.pre$Patient.no.),
                                   Response = c(as.vector(clinical_PD1.pre$Best.RECIST.response), 
                                                as.vector(clinical_PD1_CTLA4.pre$Best.RECIST.response)))
    
    tmp.response.pre$Alias <- sapply(strsplit(as.character(tmp.response.pre$Alias), split = " ", fixed = T), head, 1)
    
    clinical_PD1.on <- subset(clinical_PD1, RNA.Sequencing %in% c("EDT", "PRE and EDT"))
    clinical_PD1_CTLA4.on <- subset(clinical_PD1_CTLA4, RNA.Sequencing %in% c("EDT", "PRE and EDT"))
    
    tmp.response.on <- data.frame(Alias =  c(clinical_PD1.on$Patient.no., 
                                             clinical_PD1_CTLA4.on$Patient.no.),
                                  Response = c(as.vector(clinical_PD1.on$Best.RECIST.response), 
                                               as.vector(clinical_PD1_CTLA4.on$Best.RECIST.response)))
    
    tmp.response.on$Alias <- as.character(tmp.response.on$Alias)
    these <- which(sapply(strsplit(tmp.response.on$Alias, "_",  fixed = T), tail, 1) == "PRE and EDT")
    tmp.response.on$Alias[these] <- paste0(sapply(strsplit(tmp.response.on$Alias[these], split = "PRE ", fixed = T), head, 1),
                                    sapply(strsplit(tmp.response.on$Alias[these], split = " ", fixed = T), tail, 1))
    
    tmp.response <- rbind(tmp.response.pre, tmp.response.on)
    
    sampleinfo$Response <- rep("?", times = nrow(sampleinfo))
    sampleinfo$Alias <- as.character(sampleinfo$Alias)
    
    sampleinfo <- sampleinfo[match(tmp.response$Alias, sampleinfo$Alias),] 
    sampleinfo$Response <- tmp.response$Response
    
    pre.sampleinfo <- subset(sampleinfo, sapply(strsplit(as.character(Alias), split = "_", fixed = T), tail, 1) ==  "PRE" )
    
    order.pre <- match(rownames(pre.sampleinfo), colnames(gene.count))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions
    cell_fractions.pre <- cell_fractions[,match(rownames(pre.sampleinfo), colnames(cell_fractions))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, match(rownames(pre.sampleinfo), colnames(LRpairs))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, match(rownames(pre.sampleinfo), colnames(CYTOKINEpairs))]
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, match(rownames(pre.sampleinfo), colnames(gene.tpm))]
    }
   
    label.pre <- data.frame(Sample = rownames(pre.sampleinfo), label = pre.sampleinfo$Response)
  }
  
  if(dataset == "Liu"){
    
    gene.count <- as.data.frame(t(gene.count))
    gene.tpm <- as.data.frame(t(gene.tpm))
    
    # Complete set of pre-treatment samples follow this criteria
    pre.sampleinfo <- subset(sampleinfo, daysBiopsyAfterIpiStart == "na" & numPriorTherapies == 0 &
                               biopsyContext..1.Pre.Ipi..2.On.Ipi..3.Pre.PD1..4.On.PD1. %in% c(1,3))
                             
    order.pre <- na.omit(match(rownames(pre.sampleinfo), colnames(gene.count)))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions (*Ask francesca for new computation of all samples (n=121))
    cell_fractions.pre <- cell_fractions[,na.omit(match(rownames(pre.sampleinfo), colnames(cell_fractions)))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, na.omit(match(rownames(pre.sampleinfo), colnames(LRpairs)))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, na.omit(match(rownames(pre.sampleinfo), colnames(CYTOKINEpairs)))]
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, na.omit(match(rownames(pre.sampleinfo), colnames(gene.tpm)))]
    }
    
    label.pre <- data.frame(Sample = rownames(pre.sampleinfo)[match(colnames(gene.count.pre), rownames(pre.sampleinfo))]
                            , label = pre.sampleinfo$BR[match(colnames(gene.count.pre), rownames(pre.sampleinfo))])
  }
  
  if(dataset == "Kim"){ # n = 45 (CR=3, PR=9, PD=18, SD=15)
    
    pre.sampleinfo <- subset(sampleinfo, Assay.Type == "RNA-Seq")
    
    # By hand, supp table 2 Kim et al., Nat Med 2018
    infopaper <- data.frame(Id = paste0("PB-16-", c("045", "068", "031", "005", "010", "013", "019", "020", "021", "038",
                                                    "044", "048", "053", "059", "063", "066", "006", "016", "022", "026","034",
                                                    "037", "040", "042", "047", "049", "052", "054", "056", "058", "062", "064", "069",
                                                    "003", "014", "001", "002", "004", "011", "015", "018", "023", "024", "025", "029",
                                                    "030", "032", "035", "039", "041", "043", "051", "055", "057", "060", "067", "009")),
                            Response = c("CR", "CR", "CR", 
                                         "PR", "PR", "PR",  "PR", "PR", "PR", "PR", "PR", "PR", "PR", "PR", "PR",
                                         "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD",
                                         "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "PD", "PD", "PD",
                                         "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD",
                                         "PD", "PD", "PD", "PD", "PD", "PD", "PD"))
    
    pre.sampleinfo <- pre.sampleinfo[na.omit(match(infopaper$Id, pre.sampleinfo$title)),]
    pre.sampleinfo$Response <- infopaper$Response[match(pre.sampleinfo$title, infopaper$Id)]
    
    order.pre <- match(rownames(pre.sampleinfo), colnames(gene.count))
    gene.count.pre <- gene.count[, order.pre]
    
    # cell fractions (*Ask francesca for new computation of all samples (n=121))
    cell_fractions.pre <- cell_fractions[,match(rownames(pre.sampleinfo), colnames(cell_fractions))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, match(rownames(pre.sampleinfo), colnames(LRpairs))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, match(rownames(pre.sampleinfo), colnames(CYTOKINEpairs))]
    
    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, match(rownames(pre.sampleinfo), colnames(gene.tpm))]
    }
    
    label.pre <- data.frame(Sample = rownames(pre.sampleinfo), label = pre.sampleinfo$Response)
  }
  
  if(dataset == "Zhao"){

    clinical_PD1 <- read.csv(file = "../data/Validation/Francesca/Zhao/GBM_PD1_samples.csv",
                             sep = ",", header = T)

    clinical_PD1.pre <- subset(clinical_PD1, pre.or.post == "pre" & RNA.seq. == "Yes")
    
    clinical_PD1.pre$DNA.Sampel.Name <- gsub(" ","", as.character(clinical_PD1.pre$DNA.Sampel.Name), fixed = T)
    sampleinfo$Sample.Name <- gsub(" ","", as.character(sampleinfo$Sample.Name), fixed = T)
    
    pre.sampleinfo <- sampleinfo[na.omit(match(clinical_PD1.pre$DNA.Sampel.Name, sampleinfo$Sample.Name)),]

    order.pre <- na.omit(match(rownames(pre.sampleinfo), colnames(gene.count)))
    gene.count.pre <- gene.count[, order.pre]
    
    pre.sampleinfo.aval <- pre.sampleinfo[na.omit(match(colnames(gene.count), rownames(pre.sampleinfo))),]
    clinical_PD1.pre.aval <- clinical_PD1.pre[na.omit(match(pre.sampleinfo.aval$Sample.Name, 
                                                            clinical_PD1.pre$DNA.Sampel.Name)),]
    # cell fractions 
    cell_fractions.pre <- cell_fractions[,na.omit(match(rownames(pre.sampleinfo.aval), colnames(cell_fractions)))]
    rownames(cell_fractions.pre)[9] <- "T_cells_regulatory_Tregs"
    rownames(cell_fractions.pre) <-  gsub(".","_", rownames(cell_fractions.pre), fixed = T)
    
    # # ligand-receptor pairs
    # LRpairs.pre <- LRpairs[, na.omit(match(rownames(pre.sampleinfo.aval), colnames(LRpairs)))]
    # 
    # # cytokine pairs
    # CYTOKINEpairs.pre <- CYTOKINEpairs[, na.omit(match(rownames(pre.sampleinfo.aval), colnames(CYTOKINEpairs)))]

    if (identical(colnames(gene.count), colnames(gene.tpm))){
      gene.tpm.pre <- gene.tpm[, order.pre]
    } else{
      gene.tpm.pre <- gene.tpm[, na.omit(match(rownames(pre.sampleinfo), colnames(gene.tpm)))]
    }
    
    label.pre <- data.frame(Sample = rownames(pre.sampleinfo.aval), label = clinical_PD1.pre.aval$Response)
  }
    
  # ---------------------------------------------------------------------------------- #
  # Pathways activity computation
  # ---------------------------------------------------------------------------------- #
  
  # Correct gene names in order to match correctly pathway responsive genes
  gene.count.pre_processed_v2 <- preprocess_gene_names_to_match_pathways_genes(gene.count.pre)
  
  # Computation of Pathway scores (input matrix [genes, samples], ouput matrix [sample, pathways])
  Pathway_activities <- compute.pathways.scores(RNA.raw_counts = gene.count.pre_processed_v2,
                                                remove.genes.ICB_proxies=TRUE)

  # ---------------------------------------------------------------------------------- #
  # TF activity computation
  # ---------------------------------------------------------------------------------- #

  # Correct gene names in order to match correctly tf target genes
  gene.tpm.pre_processed_v2 <- preprocess_gene_names_to_match_tf_target_genes(gene.tpm.pre)
   
  # Computation of TF activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  TF_activities <- compute.TF.activity(RNA.tpm = gene.tpm.pre_processed_v2, 
                                       remove.genes.ICB_proxies=TRUE)

  # ---------------------------------------------------------------------------------- #
  # Immune cells fractions
  # ---------------------------------------------------------------------------------- #
  
  # Immune cells fractions (ouput matrix [sample, Immunecells])
  immunecells <- data.frame(t(cell_fractions.pre))

  # ---------------------------------------------------------------------------------- #
  # Ligand-Receptor and cytokine pairs computation
  # ---------------------------------------------------------------------------------- #
    
  # Correct gene names in order to match correctly pathway responsive genes
  gene.tpm.pre_processed_v2 <- preprocess_gene_names_to_match_LRpairs_genes(gene.tpm.pre)
  
  # Computation of L-R and cytokine pairs (ouput matrix [sample, LRpairs])
  LRpairs.data <- compute.LR.pairs(RNA.tpm = gene.tpm.pre_processed_v2, 
                                   remove.genes.ICB_proxies=TRUE, compute.cytokines.pairs=TRUE)

  LRpairs <- as.data.frame(LRpairs.data$LRpairs)
  CYTOKINEpairs <- data.frame(LRpairs.data$CYTOKINEpairs)

  # Remove NA values in DataViews LRpairs and CYTOKINE pairs (feature-wise)
  if (anyNA(LRpairs)){
    LR_sum <- apply(LRpairs, 2, sum)
    LRpairs <- LRpairs[,!is.na(LR_sum)]
    Cyt_sum <- apply(CYTOKINEpairs, 2, sum)
    CYTOKINEpairs <- CYTOKINEpairs[,!is.na(Cyt_sum)]
  }

  # ****************
  # data
  DataViews.test[[dataset]] <- list(Pathways = data.frame(Pathway_activities$scores),
                                    ImmuneCells = immunecells,
                                    TFs = data.frame(TF_activities$scores),
                                    LRpairs = LRpairs,
                                    CYTOKINEpairs = CYTOKINEpairs,
                                    Transcript = data.frame(log2(t(gene.tpm.pre)+1)))
    
    # Labels.test[[dataset]] <- label.pre
  
  return(DataViews.test)
    
}))

All.DataViews.test <- DataViews.test
# All.Labels.test <- Labels.test
#save(All.DataViews.test, file = "../data/PanCancer_draft_v1/Validation/All_DataViews_test_pre.RData")
# save(All.Labels.test, file = "../data/PanCancer_draft_v1/Validation/All_Labels_test_pre.RData")

# --------------------------- #
# Join SKCM datasets
# --------------------------- #

# Hugo, Liu, Riaz --> PD-1
# Gide, Auslander --> PD-1 & CTLA4

load("../data/PanCancer_draft_v1/Validation/All_DataViews_test_pre.RData")
load("../data/PanCancer_draft_v1/Validation/All_Labels_test_pre.RData")

# All validation datasets
Datasets.names <-  dir("../data/Validation/Francesca", full.names = F, recursive = F)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM")

# Only validation SKCM datasets
filter.cancer <- which(names(Datasets.names) == "SKCM")
Datasets.names.SKCM <- Datasets.names[filter.cancer]

# Initialize new comb.datasets
All.DataViews.test[["comb.Hugo_Liu_Riaz"]] <- vector("list", length(names(All.DataViews.test$Auslander)))
names(All.DataViews.test[["comb.Hugo_Liu_Riaz"]]) <- names(All.DataViews.test$Auslander)
All.DataViews.test[["comb.Gide_Auslander"]] <- vector("list", length(names(All.DataViews.test$Auslander)))
names(All.DataViews.test[["comb.Gide_Auslander"]]) <- names(All.DataViews.test$Auslander)

All.Labels.test[["comb.Hugo_Liu_Riaz"]] <- vector("list", length(names(All.Labels.test$Auslander)))
names(All.Labels.test[["comb.Hugo_Liu_Riaz"]]) <- names(All.Labels.test$Auslander)
All.Labels.test[["comb.Gide_Auslander"]] <- vector("list", length(names(All.Labels.test$Auslander)))
names(All.Labels.test[["comb.Gide_Auslander"]]) <- names(All.Labels.test$Auslander)

# Hugo, Liu, Riaz --> PD-1
All.DataViews.test[["comb.Hugo_Liu_Riaz"]] <- lapply(names(All.DataViews.test$Auslander), function(view){
  
  tmp.var <- intersect(colnames(All.DataViews.test[["Hugo"]][[view]]), colnames(All.DataViews.test[["Riaz"]][[view]]))
  tmp.var <- intersect(colnames(All.DataViews.test[["Liu"]][[view]]), tmp.var)
  
  All.DataViews.test[["comb.Hugo_Liu_Riaz"]][[view]] <- rbind(All.DataViews.test[["Hugo"]][[view]][,tmp.var], 
                                                              All.DataViews.test[["Liu"]][[view]][,tmp.var],
                                                              All.DataViews.test[["Riaz"]][[view]][,tmp.var])
  
})
names(All.DataViews.test[["comb.Hugo_Liu_Riaz"]]) <- names(All.DataViews.test$Auslander)

All.Labels.test[["comb.Hugo_Liu_Riaz"]] <- rbind(All.Labels.test[["Hugo"]], 
                                                 All.Labels.test[["Liu"]],
                                                 All.Labels.test[["Riaz"]])

# Gide, Auslander --> PD-1 & CTLA4
All.DataViews.test[["comb.Gide_Auslander"]] <- lapply(names(All.DataViews.test$Auslander), function(view){
  
  tmp.var <- intersect(colnames(All.DataViews.test[["Gide"]][[view]]), colnames(All.DataViews.test[["Auslander"]][[view]]))
  All.DataViews.test[["comb.Gide_Auslander"]][[view]] <- rbind(All.DataViews.test[["Gide"]][[view]][,tmp.var], 
                                                               All.DataViews.test[["Auslander"]][[view]][,tmp.var])
  
})
names(All.DataViews.test[["comb.Gide_Auslander"]]) <- names(All.DataViews.test$Auslander)

All.Labels.test[["comb.Gide_Auslander"]] <- rbind(All.Labels.test[["Gide"]], 
                                                  All.Labels.test[["Auslander"]])
                                                        
#save(All.DataViews.test, file = "../data/PanCancer_draft_v1/Validation/All_DataViews_test_pre.RData")
#save(All.Labels.test, file = "../data/PanCancer_draft_v1/Validation/All_Labels_test_pre.RData")

