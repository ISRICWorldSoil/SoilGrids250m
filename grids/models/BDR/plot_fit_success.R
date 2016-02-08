## Plotting success of fitting BDTICM using log-log scale

library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(plotKML)
library(R.utils)

gunzip("RF_fit_BDTICM.csv.gz")
fit <- read.csv("RF_fit_BDTICM.csv")
pfun <- function(x,y, ...){ 
         panel.hexbinplot(x,y, ...)  
         panel.abline(0,1,lty=1,lw=2,col="black") 
}
d.meas <- min(fit$predicted, na.rm=TRUE)
out <- fit$predicted+ifelse(d.meas==0, 1, d.meas)
meas <- fit$observed+ifelse(d.meas==0, 1, d.meas)
hexbinplot(out~meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main="Absolute depth to bedrock (cm)", xlab="measured", ylab="predicted (machine learning)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=25, density=40, xlim=range(meas, na.rm=TRUE), ylim=range(meas, na.rm=TRUE), panel=pfun, colorcut=c(0,0.005,0.01,0.03,0.07,0.15,0.25,0.5,1))