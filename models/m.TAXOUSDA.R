## Fit model for TAXOUSDA

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(dplyr)
library(snowfall)
library(rgdal)
library(nnet)
library(psych)

source("extract.equi7t3.R")
source("wrapper.predict_c.R")
load("/data/covs/equi7t3.rda")
des <- read.csv("/data/covs/SoilGrids250m_COVS250m.csv")
load("TAXOUSDA.pnts.rda")
ov <- extract.equi7t3(x=TAXOUSDA.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40)
## TAKES >20 MINS FOR 30k points
str(ov)
ov$LATWGS84 <- as.numeric(sapply(paste(ov$LOC_ID), function(x){strsplit(x, "_")[[1]][2]}))
## 30,270 profiles
write.csv(ov, file="ov.TAXOUSDA_SoilGrids250m.csv")
library(utils)
zip("ov.TAXOUSDA_SoilGrids250m.csv", zipfile="ov.TAXOUSDA_SoilGrids250m.zip")
save(ov, file="ov.TAXOUSDA.rda")
pr.lst <- des$WORLDGRIDS_CODE
formulaString.USDA = as.formula(paste('TAXOUSDA.f ~ LATWGS84 + ', paste(pr.lst, collapse="+")))
formulaString.USDA
## TAKES > 20 mins to fit...
## can it be parallelized?
m_TAXOUSDA <- nnet::multinom(formulaString.USDA, ov, MaxNWts = 9000)
# groups "Anthrepts", "Ustox" are empty
str(fitted(m_TAXOUSDA))
## ?? points
head(signif(fitted(m_TAXOUSDA),3))
## goodness of fit:
cout.m <- as.factor(paste(predict(m_TAXOUSDA, newdata=ov, na.action = na.pass)))
cf <- mda::confusion(cout.m, as.character(ov[,"TAXOUSDA.f"]))
## remove missing classes:
a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
c.kappa = psych::cohen.kappa(cf[a,b])
ac <- sum(diag(cf))/sum(cf)*100
message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
## 38%
message(paste("Map purity:", signif(ac, 3)))
## 37%
save(m_TAXOUSDA, file="m_TAXOUSDA.rda")

## Alternative models??
#mrf_TAXOUSDA <- randomForestSRC::rfsrc(formulaString.USDA, ov)

## predict for sample locations:
dir.lst <- list.dirs("/data/covs")[-1]
Sys.chmod(gsub("covs", "predicted", dir.lst), "777")
wrapper.predict_c(i="NA_060_036", varn="TAXOUSDA", gm=m_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted")
wrapper.predict_c(i="OC_087_063", varn="TAXOUSDA", gm=m_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted")
wrapper.predict_c(i="EU_051_012", varn="TAXOUSDA", gm=m_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted")
## plot in GE:
plotKML.env(convert="convert", show.env=FALSE)
x <- readGDAL("/data/predicted/NA_060_036/TAXOUSDA_Xeralfs_NA_060_036.tif")
kml(x, file.name="TAXOUSDA_Xeralfs_NA_060_036.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_NA_060_036.png")
x <- readGDAL("/data/predicted/OC_087_063/TAXOUSDA_Xeralfs_OC_087_063.tif")
kml(x, file.name="TAXOUSDA_Xeralfs_OC_087_063.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_OC_087_063.png")
x <- readGDAL("/data/predicted/EU_051_012/TAXOUSDA_Aqualfs_EU_051_012.tif")
kml(x, file.name="TAXOUSDA_Aqualfs_EU_051_012.kml", folder.name="Aqualfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Aqualfs_EU_051_012.png")

## run in parallel:


