# Plot correlations including p-value, correlation value and correlation line 
panel.cor <- function(x, y, digits=2, cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x,y)$estimate
  p <- cor.test(x,y)$p.value
  txt_r <- format(r, digits=digits)
  txt_p <- format(p, scientific = TRUE, digits=digits)
  txt <- paste("cor=", txt_r, "\npVal=", txt_p, sep="")
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "#A1A1A1", ...) {
  
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  # abline(stats::lm(y ~ x),  col = col.smooth, ...)
  abline(a=0, b=1,  col = col.smooth, ...)
}

# Check correlation between features
pairs( ~ ., data = DataViews$PROGENy, upper.panel = panel.cor,lower.panel = panel.lm) 