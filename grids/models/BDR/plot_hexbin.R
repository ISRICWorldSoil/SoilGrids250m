## Plotting results of cross-validation:
## Tom.Hengl@isric.org

plot_hexbin <- function(j, breaks, main, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1), pal=R_pal[["bpy_colors"]][1:18], in.file, log.plot){
  out.file = paste0("plot_CV_", t.vars[j], ".png")
  if(!file.exists(out.file)){
    load(in.file)
    assign("m", get(paste0("CV_", t.vars[j]))[[1]])
    png(file = out.file, res = 150, width=850, height=850, type="cairo")
    if(log.plot==TRUE){
      d.meas <- min(m$Observed, na.rm=TRUE)
      out <- m$Predicted+ifelse(d.meas==0, 1, d.meas)
      meas <- m$Observed+ifelse(d.meas==0, 1, d.meas)
      lim <- range(breaks)+ifelse(d.meas==0, 1, d.meas)
      plt <- hexbinplot(out~meas, colramp=colorRampPalette(pal), main=main, xlab="measured", ylab="predicted (SoilGrids250m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=25, density=40, xlim=lim, ylim=lim, panel=pfun, colorcut=colorcut)
    } else {
      lim <- range(breaks)
      plt <- hexbinplot(m$Predicted~m$Observed, colramp=colorRampPalette(pal), main=main, xlab="measured", ylab="predicted (SoilGrids250m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, xlim=lim, ylim=lim, asp=1, xbins=25, density=40, panel=pfun, colorcut=colorcut)
    }
    print(plt)
    dev.off()
  }
}