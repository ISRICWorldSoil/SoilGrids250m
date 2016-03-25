rm(list = ls(all = TRUE))
library(sp)
library(maps)
library(grDevices)
library(maptools)
library(lattice)
#par(resetPar())
# global define
#gdal.dir = shortPathName("C:\\ms4w")
#gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)
a.dir <- "/data/shang009/big"# dir of the project
m.dir <- paste0(a.dir, "/soildepth2")
setwd(m.dir)
source(paste0(a.dir, "/soildepth/code/head/functions.r"))


resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}
Partmp <- resetPar() 
par(Partmp)
soilwell <- read.csv("ov.BDR_SoilGrids250m.csv")
soilwell <- soilwell[,1:8]
simulated <- subset(soilwell, soilwell$SOURCEDB== "Simulated")
soilwell <- subset(soilwell, soilwell$SOURCEDB!= "Simulated")
soilwell2 <- subset(soilwell, soilwell$BDRICM<250)


####EDA for BDTICM
bitmap(paste0("./pics/", "Hist_BDT.tif"), width = 7.48, height = 4, units = "in", res =1000, type = "tiffcrle", pointsize =11)
#bitmap(paste0("./pics/", "Hist_BDT.png"), width = 7.48, height = 4, units = "in", res =1000, type = "pngmono", pointsize =11)
par(mar = c(4,4,0,2), mfcol=c(1,2))
hist(soilwell$BDTICM, main = "", xlim =  quantile(soilwell$BDTICM, probs = c(0, 0.99), na.rm = T), breaks =1000 , xlab = "Absolute Depth to bedrock (cm)", col="grey")
hist(log1p(soilwell$BDTICM), main = "", xlab = "Log-transformed absolute depth to bedrock (cm)", col="grey")
dev.off()
par(mar = c(5.1,4.1,4.1,2.1), mfcol=c(1,1))
summary(soilwell$BDTICM)
length(soilwell$BDTICM)- sum(is.na(soilwell$BDTICM))
#1586581
####EDA for BDRICM
bitmap(paste0("./pics/", "Hist_BDR.tif"), width = 7.48, height = 4, units = "in", res =1000, type = "tiffcrle", pointsize =11)
#bitmap(paste0("./pics/", "Hist_BDR.png"), width = 7.48, height = 4, units = "in", res =1000, type = "pngmono", pointsize =11)
par(mar = c(4,4,0,0), mfcol=c(1,2))
hist(soilwell2[soilwell2$SOURCEDB!= "Wells",]$BDRICM, main = "", xlab = "Censored depth to bedrock from soils (cm)", breaks= 20 , col="grey" )
hist(soilwell2[soilwell2$SOURCEDB == "Wells",]$BDRICM, main = "", xlab = "Censored depth to bedrock from wells (cm)", breaks= 20, col="grey" )
dev.off()
par(mar = c(5.1,4.1,4.1,2.1), mfcol=c(1,1))
summary(soilwell2$BDRICM)
sum(soilwell$BDRICM==250)

sum(!is.na(soilwell2[soilwell2$SOURCEDB!= "Wells",]$BDRICM))


sum(!is.na(soilwell2[soilwell2$SOURCEDB== "Wells",]$BDRICM))

