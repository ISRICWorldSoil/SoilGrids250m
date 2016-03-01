library(plyr)
setwd("~/big/soildepth2/code")
t.vars <- c("BDTICM", "BDRICM", "BDRLOG")
plt.names <- c("Absolute depth to bedrock", "Censored depth to bedrock", "Occurrence of the R horizon")
vn1 <- NULL
for(i in 1 :3)
{
  vn1 <- c(vn1, read.table(paste0("../imp_",t.vars[i],".txt"),as.is=T)[,1])
 
}
vn1 <-sort(unique(vn1))
vn1
vn2 <- c("DTB by Pelletier", "Land surface elevation","Surface roughness", "MODIS NIR band 4 Apr.", "MODIS NIR band 4 May", "MODIS NIR band 4 Oct.", "Latitude", "MODIS MIR band 7 Apr.", "MODIS MIR band 7 May",  "MODIS MIR band 7 Oct.", "precipitation Jan.", "precipitation May", "precipitation Jul.", "precipitation Oct.", "precipitation Nov.", "Total annual precipitation", "Slope", "MODIS daytime LST Mar.", "MODIS daytime LST Apr.", "MODIS daytime LST May" , "MODIS daytime LST Jun.", "MODIS daytime LST Sep.",  "MODIS daytime LST Oct.",   "MODIS daytime LST Nov.","Wetness index", "MRVBF", "Valley depth", "PWV May to Jun.", "PWV Nov. to Dec.")
vn <- data.frame(variable=vn1,vn=vn2)
#"Multiresolution Index of Valley Bottom Flatness"
#"Precipitable Water Vapor"
bitmap(paste0("../pics/Fig_imp.tif"), width = 7.48, height = 8, units = "in", res =300, type = "tiffcrle", pointsize =11)
#bitmap(paste0("../pics/Fig_imp.png"), width = 7.48, height = 8, units = "in", res =1000, type = "pngmono", pointsize =11)
par(mfrow=c(2, 2), mar=c(4, 1, 3, 1))
for(i in 1 :3)
{
  imp <- read.table(paste0("../imp_",t.vars[i],".txt"))
  imp <-join(imp, vn, type ="left")
  row.names(imp) <- imp[,5]
  imp$variable <-NULL
  imp$vn <-NULL
  imp <- as.matrix(imp)
  
  ord <-  rev(order(imp[,2], decreasing=TRUE))
  dotchart(imp[ord,2], xlab="Scaled importance", ylab="", xlim=c(0, max(imp[,2])),main=plt.names[i])

#"residual sum of squares"
}
dev.off()



