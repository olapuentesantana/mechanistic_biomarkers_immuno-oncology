# #########################################################################################################
# Script to check correlation between tasks: proxy's of response to ICBs
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
load("TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)


comb_tasks <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load previous data (remove all)
  load(paste0("/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/",
              Cancer,"/new/ImmuneResponse_no_filter_", Cancer,"_matrix_format.RData"))  
  
  tmp_data <- as.matrix(ImmuneResponse.no_filter)
  
  return(tmp_data)
  
}))

data <- as.data.frame(comb_tasks)

pairs( ~ . , data = data, upper.panel = panel.cor,lower.panel = panel.lm,
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
