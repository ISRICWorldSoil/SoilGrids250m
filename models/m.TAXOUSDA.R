## Fit model for TAXOUSDA
## By Tom.Hengl@isric.org 
## Contributions by Mario Guevara <mguevara@udel.edu>

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(dplyr)
library(snowfall)
library(rgdal)
library(nnet)
library(randomForest)
library(kknn)
library(psych)
library(plotKML)
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
## Color legend:
col.legend <- read.csv("TAXOUSDA_legend.csv")
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)

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

## Alternative models
mrf_TAXOUSDA <- randomForest(formulaString.USDA, ov) ## very fast!
mrf_TAXOUSDA
save(mrf_TAXOUSDA, file="mrf_TAXOUSDA.rda")
#ss <- sample(1:nrow(ov), size=1000)
mkk_TAXOUSDA <- kknn(formulaString.USDA, train=ov, test=ov, distance=1, kernel="triangular")
save(mkk_TAXOUSDA, file="mkk_TAXOUSDA.rda")

## create dirs:
dir.lst <- list.dirs("/data/covs")[-1]
x <- lapply(gsub("covs", "predicted", dir.lst), dir.create, recursive=TRUE)
## predict for sample locations:
wrapper.predict_c(i="NA_060_036", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="../covs", out.path="../predicted", col.legend=col.legend)
wrapper.predict_c(i="OC_087_063", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="../covs", out.path="../predicted", col.legend=col.legend)
wrapper.predict_c(i="EU_051_012", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="../covs", out.path="../predicted", col.legend=col.legend)

## plot in GE:
x <- readGDAL("/data/predicted/NA_060_036/TAXOUSDA_Xeralfs_NA_060_036.tif")
kml(x, file.name="TAXOUSDA_Xeralfs_NA_060_036.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_NA_060_036.png")
x <- readGDAL("/data/predicted/NA_060_036/TAXOUSDA_NA_060_036.tif")
x$cl <- col.legend[match(x$band1, col.legend$Number),"Group"]
x$cl <- as.factor(paste(x$cl))
pal <- col.legend[match(levels(x$cl), col.legend$Group),"COLOR"]
#raster::image(raster(x["cl"]), col=pal)
kml(x, file.name="TAXOUSDA_NA_060_036.kml", folder.name="SoilGrids: TAXOUSDA", subfolder.name="Predicted", colour=cl, colour_scale=pal, plot.legend=FALSE, raster_name="TAXOUSDA_NA_060_036.png")
x <- readGDAL("/data/predicted/OC_087_063/TAXOUSDA_Xeralfs_OC_087_063.tif")
kml(x, file.name="TAXOUSDA_Xeralfs_OC_087_063.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_OC_087_063.png")
x <- readGDAL("/data/predicted/EU_051_012/TAXOUSDA_Aqualfs_EU_051_012.tif")
kml(x, file.name="TAXOUSDA_Aqualfs_EU_051_012.kml", folder.name="Aqualfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Aqualfs_EU_051_012.png")
x <- readGDAL("/data/predicted/EU_051_012/TAXOUSDA_EU_051_012.tif")
x$cl <- col.legend[match(x$band1, col.legend$Number),"Group"]
x$cl <- as.factor(paste(x$cl))
pal <- col.legend[match(levels(x$cl), col.legend$Group),"COLOR"]
#raster::image(raster(x["cl"]), col=pal)
kml(x, file.name="TAXOUSDA_EU_051_012.kml", folder.name="SoilGrids: TAXOUSDA", subfolder.name="Predicted", colour=cl, colour_scale=pal, plot.legend=FALSE, raster_name="TAXOUSDA_EU_051_012.png")

## run all predictions in parallel:
pr.dirs <- basename(list.dirs("/data/predicted")[-1])
## 2594 dirs
sfInit(parallel=TRUE, cpus=35)
sfExport("wrapper.predict_c", "m_TAXOUSDA", "mrf_TAXOUSDA", "col.legend")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(plyr)
sfLibrary(nnet)
sfLibrary(randomForest)
x <- sfLapply(pr.dirs, fun=function(x){ try( wrapper.predict_c(x, varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend) ) } )

## Create mosaicks:
x <- lapply(paste0("/data/GEOG/", names(equi7t3)), dir.create, recursive=TRUE)
mosiack.equi7t3 <- function(i, j){
  out.tif <- paste0('/data/GEOG/', j, '/TAXOUSDA_', i, '_', j, '_250m_ll.tif')
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("TAXOUSDA_", i, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"bilinear\" -ot \"Byte\" -dstnodata \"255\" -tr 0.002083333 0.002083333 -co \"COMPRESS=DEFLATE\"'))
  }
}
## Bounding boxes need to be set manually?
#NA: -150 12 -50 60

for(j in names(equi7t3)){
  sfInit(parallel=TRUE, cpus=35)
  sfExport("mosiack.equi7t3", "j", "gdalbuildvrt", "gdalwarp")
  x <- sfLapply(m_TAXOUSDA$lev, fun=mosiack.equi7t3, j=j)
  sfStop()
}

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("TAXOUSDA_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
