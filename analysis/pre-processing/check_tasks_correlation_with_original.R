# #########################################################################################################
# Script to check correlation between computed tasks and original scores
# #########################################################################################################

# ****************
# working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/")

# ****************
# functions 
panel.cor <- function(x, y, digits=2, font.cor = 1, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x,y)$estimate
  p <- cor.test(x,y)$p.value
  txt_r <- format(r, digits=digits)
  txt_p <- format(p, scientific = TRUE, digits=digits)
  txt <- paste("cor=", txt_r, "\np=", txt_p, sep="")
  
  if(txt_r >= 0.7 & txt_p >= 0.05) font.cor <- 2
  
  text(0.5, 0.5, txt, cex = 1, font = font.cor)
}

panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.8, col.smooth = "#A1A1A1", ...) {
  
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  # abline(stats::lm(y ~ x),  col = col.smooth, ...)
  abline(a=0, b=1,  col = col.smooth, ...)
}

# ****************
# Select cancer type
## no filter
load("analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

# ****************
# Original data
## CYT
CYT.paper <- read.csv("data/raw_data_tcga/cytolyticActivity.csv", row.names = 1, header = T)
CYT.paper <- CYT.paper[,"Cytolytic.Activity", drop = FALSE]
CYT.paper <- CYT.paper[!is.na(CYT.paper$Cytolytic.Activity), , drop = FALSE]

## IS
IS.paper <- read.csv("data/raw_data_tcga/immuneSignature.csv", row.names = 1, header = T, skip = 1)
IS.paper <- IS.paper[-c(9082:9083),2, drop = FALSE]
IS.paper <- IS.paper[!is.na(IS.paper$Immune_signature_score), , drop = FALSE]

## IPS
IPS.paper <- read.csv("data/raw_data_tcga/patientsAll_IPS_literature.csv", row.names = 1, header = TRUE, sep = "\t")
IPS.paper <- IPS.paper[,"ips_ctla4_neg_pd1_neg", drop = FALSE]

## IMPRES
IMPRES.paper <- read.csv("data/raw_data_tcga/IMPRES_TCGA.csv", row.names = 1, header = T)
IMPRES.paper <- IMPRES.paper[,2, drop = FALSE]

## RohIS

## Chemokine

## IS_Davoli
Davoli.paper <- read.csv("data/raw_data_tcga/IS_Proliferation_Davoli_TCGA.csv", header = T, skip = 5, row.names = 1)
IS_Davoli.paper <- Davoli.paper[,"Immune.Signature.Score", drop = FALSE]
IS_Davoli.paper <- IS_Davoli.paper[!is.na(IS_Davoli.paper[,1]), , drop = FALSE]

## Proliferation
Proliferation.paper <- Davoli.paper[,"CellCycle.Signature.Score", drop = FALSE]
Proliferation.paper <- Proliferation.paper[!is.na(Proliferation.paper[,1]), , drop = FALSE]

## IFny
IFny.paper <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
    TIDE.table <- read.table(file = paste0("data/data_processed_TIDE/output_TIDE_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
    IFny.paper <- TIDE.table[keep.samples.no_filter, "IFNG", drop = FALSE]
    return(IFny.paper)
}))

## ExpandedImmune

## T_cell_inflamed

## TIDE
TIDE.paper <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  TIDE.table <- read.table(file = paste0("data/data_processed_TIDE/output_TIDE_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
  TIDE.paper <- TIDE.table[keep.samples.no_filter, "TIDE", drop = FALSE]
  return(TIDE.paper)
}))

## MSI
MSI.paper <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  TIDE.table <- read.table(file = paste0("data/data_processed_TIDE/output_TIDE_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
  MSI.paper <- TIDE.table[keep.samples.no_filter, "MSI.Score", drop = FALSE]
  return(MSI.paper)
}))

# Collect computed tasks over all cancer types
comb_tasks <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load previous data (remove all)
  load(paste0("data/PanCancer/", Cancer,"/new/ImmuneResponse_no_filter_", Cancer,"_matrix_format.RData"))  
  tmp_data <- as.matrix(ImmuneResponse.no_filter)
  
  return(tmp_data)
}))

all.computed.scores <- as.data.frame(comb_tasks)
all.computed.scores <- all.computed.scores[,-c(2,5,6,10,11)]

# Dataframe with all original scores
original.scores <- list(CYT = CYT.paper, 
                        IPS = IPS.paper,
                        IMPRES = IMPRES.paper, 
                        IS_Davoli = IS_Davoli.paper, 
                        Proliferation = Proliferation.paper,
                        IFny = IFny.paper,
                        TIDE = TIDE.paper,
                        MSI = MSI.paper)

all.original.scores.list <- list()
all.original.scores.list <- lapply(names(original.scores), function(score){
  # Match samples
  tmp <- intersect(substr(rownames(all.computed.scores), 1, 12), substr(rownames(original.scores[[score]]), 1, 12))
  tmp.data <- original.scores[[score]][tmp, , drop = FALSE]
  all.original.scores.list[[score]] <- tmp.data
})
names(all.original.scores.list) <- names(original.scores)

all.computed.scores.list <- list()
all.computed.scores.list <- lapply(colnames(all.computed.scores), function(score){
  # Match samples
  tmp <- intersect(substr(rownames(all.computed.scores), 1, 12), substr(rownames(original.scores[[score]]), 1, 12))
  tmp.data <- all.computed.scores[tmp,score , drop = FALSE]
  all.computed.scores.list[[score]] <- tmp.data
})
names(all.computed.scores.list) <- colnames(all.computed.scores)

# ****************
# Check correlation

data <- data.frame(original = all.original.scores.list[[8]], computed = all.computed.scores.list[[8]])
pairs( ~ ., data = data, upper.panel = panel.cor,lower.panel = panel.lm,
       cex.labels = 1)

cor.matrix <- cor(data)
cor.sig <- cor.mtest(data)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor.matrix, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.8, number.cex = 0.75, #Text label color and rotation
         # Combine with significance
         p.mat = cor.sig$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)