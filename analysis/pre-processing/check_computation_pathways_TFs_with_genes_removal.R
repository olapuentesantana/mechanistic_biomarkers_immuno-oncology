# #########################################################################################################
# Script to visualize the data
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

# ****************
# views
views <- c(pathways = 'gaussian', #1
           Protall = 'gaussian', #2
           immunecells = 'gaussian', #3
           TFs = 'gaussian', #4
           transcript = 'gaussian', #5
           sTIL = 'gaussian', #6
           LRpairs = 'gaussian', #7
           CYTOKINEpairs = 'gaussian')  #8) 

# **************** 
# Select data to examine

## Pathways ## 
view_combinations <- list(views[1])

## TFs ## 
#view_combinations <- list(views[4])

# Check
comb_remove_IS_CYT_IPS_IMPRES_genes <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load previous data (remove IS,CYT,IPS,IMPRES)
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/repository/mechanistic_biomarkers_immuno-oncology/data/PanCancer/",
              Cancer,"/new/DataViews_no_filter_", Cancer,".RData"))  
  a <- DataViews.no_filter[[names(view_combinations[[1]])]]
  
  tmp_data <- as.matrix(a)

  return(tmp_data)
  
}))

comb_remove_all_genes <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load previous data (remove all)
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/repository/mechanistic_biomarkers_immuno-oncology/data/PanCancer/",
              Cancer,"/new_remove_all_genes/DataViews_no_filter_", Cancer,".RData"))  
  b <- DataViews.no_filter[[names(view_combinations[[1]])]]
  
  tmp_data <- as.matrix(b)

  return(tmp_data)
  
}))


comb_keep_all_genes <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load previous data (keep_all)
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/repository/mechanistic_biomarkers_immuno-oncology/data/PanCancer/",
              Cancer, "/new_keep_all_genes/DataViews_no_filter_",Cancer,".RData"))  
  c <- DataViews.no_filter[[names(view_combinations[[1]])]]
  
  tmp_data <- as.matrix(c)
  
  return(tmp_data)
  
}))

keep_common_features <- colnames(comb_remove_all_genes)

data <- data.frame(remove_IS_CYT_IPS_IMPRES_genes = as.vector(comb_remove_IS_CYT_IPS_IMPRES_genes[,keep_common_features]),
                   remove_all_genes = as.vector(comb_remove_all_genes[,keep_common_features]),
                   keep_all_genes = as.vector(comb_keep_all_genes[,keep_common_features]))

pairs( ~ . , data = data, upper.panel = panel.cor,lower.panel = panel.lm,
       cex.labels = 1)



