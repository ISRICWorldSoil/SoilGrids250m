## Fit models for USDA GreatGroups, texture classes, and several soil properties using the NASIS points (ca 350,000)
## Code by Tom.Hengl@isric.org / points prepared by Amanda Rachmaran () and Travis Nauman (tnauman@usgs.gov)

setwd("/data/NASIS")
load(".RData")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
library(caret)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger") ## version to deal with Memory problems
library(ranger)
library(xgboost)
library(nnet)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(raster)
library(plotKML)
library(GSIF)
library(randomForestSRC)
library(parallel)
options(rf.cores=detectCores(), mc.cores=detectCores())

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "gdal_translate"
  gdalwarp =  "gdalwarp"
  gdalbuildvrt = "gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
}

source("/data/models/saveRDS_functions.R")
source("functions_NASIS.R")

## legends
des <- read.csv("SoilGrids_USA48_Covs100m.csv")
LNDCOV6.leg <- read.csv("/data/GEOG/NASIS/SoilGrids_USA48_LandCover.csv")
PMTGSS7.leg <- read.csv("/data/GEOG/NASIS/SoilGrids_USA48_gSSURGO_pmaterial.csv")
DRNGSS7.leg <- read.csv("/data/GEOG/NASIS/SoilGrids_USA48_gSSURGO_drainage.csv")
PVEGKT6.leg <- read.csv("/data/GEOG/NASIS/potential_vegetation_legend.csv")
COUNTY6.leg <- read.csv("/data/GEOG/NASIS/Counties_legend.csv")
#GEOUSG6.leg <- read.csv("/data/GEOG/NASIS/geology_legend.csv")
GESUSG6.leg <- read.csv("/data/GEOG/NASIS/surfacegeology_legend.csv")

## Training points (taxonomy, texture classes):
TAX_gg.pnts = readRDS("TAX_gg.pnts.rds")
str(TAX_gg.pnts@data)
NASISpscs.pnts = readRDS("NASISpscs.pnts.rds")
NASISpscs.pnts$LOC_ID = paste("ID", NASISpscs.pnts$xwgs84, NASISpscs.pnts$ywgs84, sep="_")

## Soil properties NRCS NCSS:
system("7za e Geo_Peds07272016.csv.gz")
NCSS_peds = read.csv("Geo_Peds07272016.csv")
#summary(as.factor(NCSS_peds$geoid))
NCSS_peds$LOC_ID <- paste("ID", NCSS_peds$long, NCSS_peds$lat, sep="_")
t.props = c("clay", "sand", "soc", "bd", "n_tot") 
## clay, sand, soil organic carbon, total nitrogen, bulk density
lapply(t.props, function(i){summary(NCSS_peds[,i])})
NCSS_peds = NCSS_peds[,c("pedon_key","LOC_ID","long","lat","hzn_top","hzn_bot",t.props)]
str(NCSS_peds)

## Unique locations:
pnts = data.frame(LOC_ID=unique(c(TAX_gg.pnts$LOC_ID, NASISpscs.pnts$LOC_ID, NCSS_peds$LOC_ID)), LONWGS84=NA, LATWGS84=NA)
pnts$LONWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][2]}))
pnts$LATWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][3]}))
pnts = pnts[pnts$LONWGS84< -59 & pnts$LONWGS84> -127 & pnts$LATWGS84< 52 & pnts$LATWGS84> 24,]
str(pnts)
## 305,124 
coordinates(pnts) <- ~ LONWGS84 + LATWGS84
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")

## Tiling system:
obj <- GDALinfo("/data/GEOG/NASIS/COUNTY6.tif")
tile.lst <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=TRUE)
proj4string(tile.lst)
tile.tbl <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## 6300 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
writeOGR(tile.pol, "tiles50km.shp", "tiles50km", "ESRI Shapefile")
plot(spTransform(pnts, CRS(proj4string(tile.lst))), pch="+")
lines(tile.lst, col="red")
## Overlay tiles and admin units (fully parallelized):
system(paste('gdal_translate /data/GEOG/NASIS/COUNTY6.tif /data/USA48/COUNTY6.sdat -ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
system(paste0(saga_cmd, ' -c=48 shapes_grid 2 -GRIDS=\"/data/USA48/COUNTY6.sgrd\" -POLYGONS=\"tiles50km.shp\" -PARALLELIZED=1 -RESULT=\"ov_COUNTY6_tiles50km.shp\"'))
ov_COUNTY6 = readOGR("ov_COUNTY6_tiles50km.shp", "ov_COUNTY6_tiles50km")
summary(sel.t <- !ov_COUNTY6@data[,"COUNTY6..MA"]==0)
## 3390 tiles with values
ov_COUNTY6 = ov_COUNTY6[sel.t,]
str(ov_COUNTY6@data)

## Prepare covariates as tiles:
t.sel = as.character(ov_COUNTY6$ID)
new.dirs <- paste0("/data/NASIS/covs100m/T", t.sel)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
x <- lapply(gsub("covs100m", "predicted100m", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## run in parallel:
covs.lst = as.character(des$WORLDGRIDS_CODE[-which(des$WORLDGRIDS_CODE=="COUNTY6")])
#make_RDS_tiles(i=as.numeric(t.sel)[1], tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg)

#rds.lst = list.files(path="./covs100m", pattern=glob2rx("*.rds"), recursive = TRUE, full.names = TRUE)
#unlink(rds.lst)
#make_RDS_tiles(i=1922, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)
## TAKES 5-7 HOURS TO PREPARE
sfInit(parallel=TRUE, cpus=25)
sfExport("make_RDS_tiles", "tile.tbl", "covs.lst", "LNDCOV6.leg", "PMTGSS7.leg", "DRNGSS7.leg", "PVEGKT6.leg", "COUNTY6.leg", "GESUSG6.leg")
sfLibrary(plyr)
sfLibrary(rgdal)
sfLibrary(randomForestSRC)
sfLibrary(sp)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( make_RDS_tiles(i, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg) ) } )
sfStop()
save.image()

## OVERLAY (takes ca 15 mins):
ov <- extract.tiled(x=pnts, tile.pol=tile.pol, path="/data/NASIS/covs100m", ID="ID", cpus=48)
ov$PMTGSS7 = as.factor(ov$PMTGSS7)
ov$DRNGSS7 = as.factor(ov$DRNGSS7)
str(ov)
## "PMTGSS7" has too many levels
## http://stackoverflow.com/questions/14921805/convert-a-factor-to-indicator-variables
PMTGSS7.mat = data.frame(model.matrix(~PMTGSS7, ov))
str(PMTGSS7.mat)
PMTGSS7.mat$X.Intercept. = NULL
xsum.PMTGSS7 = sapply(PMTGSS7.mat, sum, na.rm=TRUE)
## 10 categories have <10 observations
sel.rm.PMTGSS7 = names(PMTGSS7.mat)[which(xsum.PMTGSS7<6)]
## Land cover:
summary(ov$LNDCOV6)
LNDCOV6.mat = data.frame(model.matrix(~LNDCOV6, ov))
LNDCOV6.mat$X.Intercept. = NULL
xsum.LNDCOV6 = sapply(LNDCOV6.mat, sum, na.rm=TRUE)
sel.rm.LNDCOV6 = names(LNDCOV6.mat)[which(xsum.LNDCOV6<6)]
## Potential vegetation:
PVEGKT6.mat = data.frame(model.matrix(~PVEGKT6, ov))
PVEGKT6.mat$X.Intercept. = NULL
xsum.PVEGKT6 = sapply(PVEGKT6.mat, sum, na.rm=TRUE)
sel.rm.PVEGKT6 = names(PVEGKT6.mat)[which(xsum.PVEGKT6<6)]
## Drainage class:
summary(ov$DRNGSS7)
DRNGSS7.mat = data.frame(model.matrix(~DRNGSS7, ov))
DRNGSS7.mat$X.Intercept. = NULL
xsum.DRNGSS7 = sapply(DRNGSS7.mat, sum, na.rm=TRUE)
sel.rm.DRNGSS7 = names(DRNGSS7.mat)[which(xsum.DRNGSS7<6)]

## Merge all indicators with other covariates
ov.fs = plyr::join_all(list(ov, 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PMTGSS7.mat))]), PMTGSS7.mat[,names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(LNDCOV6.mat))]), LNDCOV6.mat[,names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PVEGKT6.mat))]), PVEGKT6.mat[,names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(DRNGSS7.mat))]), DRNGSS7.mat[,names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)]])))
## join sampled values and covariates
ovA_NCSS_peds = plyr::join(NCSS_peds, ov.fs, by="LOC_ID", type="left")
ovA_TAX_gg.pnts = plyr::join(as.data.frame(TAX_gg.pnts), ov.fs, by="LOC_ID", type="left")
ovA_NASISpscs.pnts = plyr::join(as.data.frame(NASISpscs.pnts), ov.fs, by="LOC_ID", type="left")
saveRDS.gz(ovA_NCSS_peds, file="ovA_NCSS_peds.rds")
saveRDS.gz(ovA_TAX_gg.pnts, file="ovA_TAX_gg.pnts.rds")
saveRDS.gz(ovA_NASISpscs.pnts, file="ovA_NASISpscs.pnts.rds")

## covariates - WE DO NOT USE ANY FACTORS ONLY INDICATORS
pr.lst <- des$WORLDGRIDS_CODE
## clean up list of covariates
pr.lst = c(paste(pr.lst[-c(grep(pattern = "GESUSG6",pr.lst), grep(pattern = "PMTGSS7",pr.lst), grep(pattern = "PVEGKT6",pr.lst), grep(pattern = "LNDCOV6",pr.lst), grep(pattern = "DRNGSS7",pr.lst), grep(pattern = "COUNTY6",pr.lst), grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst), grep(pattern = glob2rx("N??PR5"),pr.lst), grep(pattern = glob2rx("N??PRI5"),pr.lst), grep(pattern = glob2rx("V??PRI5"),pr.lst))]), names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)], names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)], names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)], names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)])
## 205 covs

## FIT MODELS FOR FACTORS:
t.fvars = c("soiltype", "textype")
t.names = list("USDA Great Group", "SPCS classes")
names(t.names) = t.fvars
ov.lst = list(ovA_TAX_gg.pnts, ovA_NASISpscs.pnts)
names(ov.lst) = t.fvars
save.image()
## fit models in a loop:
for(j in t.fvars){
  if(!file.exists(paste0("mRF_", j,".rds"))){
    formulaString = as.formula(paste(j, ' ~ ', paste(pr.lst, collapse="+")))
    x = ov.lst[[j]][,all.vars(formulaString)] ## sample.int(nrow(ov.lst[[j]]), 2e4)
    x[,j] = droplevels(x[,j])
    x = x[complete.cases(x),]
    mRF = rfsrc(formulaString, data=x, ntree = 500, membership=TRUE) ## ov.lst[[j]][sample.int(nrow(ov.lst[[j]]), 2e4),all.vars(formulaString)]
    #mRF = ranger::ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, probability=TRUE) ## respect.unordered.factors=TRUE
    #save(formulaString, x, file="classification.RData")
    ## https://www.r-bloggers.com/on-ranger-respect-unordered-factors/
    gc(); gc()
    cat("Results of model fitting 'randomForest':\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"))
    cat("\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    cat(paste("Variable:", t.names[[j]]), file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    cat("\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    sink(file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE, type="output")
    cat("\n Random forest model:", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    print(mRF)
    cat("\n Variable importance:\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    s.mRF = rfsrc(formulaString, data=x[sample.int(nrow(x), 5e4),], ntree = 500)
    xl <- vimp(s.mRF)
    print(data.frame(xl$importance)[,1], max.levels=length(all.vars(formulaString)))
    write.csv(data.frame(xl$importance), file=paste0("vimp_mRF_", j,".csv"))
    #xl <- as.list(ranger::importance(mRF))
    #print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
    sink()
    saveRDS.gz(mRF, file=paste0("mRF_", j,".rds"))
  }
}

## test prediction:
## Check:
m = readRDS("./covs100m/T1049/T1049.rds")
str(m@data)
spplot(m["LNDCOV6"])
#spplot(m["PMTGSS7"])
#spplot(m["PMTGSS7.f"])
spplot(m["EX2MOD5"])
#spplot(m["DRNGSS7"])
mRF = readRDS.gz("mRF_soiltype.rds")
m$PMTGSS7 = factor(m$PMTGSS7, levels=unique(PMTGSS7.leg$pmaterial_class_f))
m$DRNGSS7 = factor(m$DRNGSS7, levels=unique(DRNGSS7.leg$drainage_class))
m.PVEGKT6 = data.frame(model.matrix(~PVEGKT6-1, m@data))
m.LNDCOV6 = data.frame(model.matrix(~LNDCOV6-1, m@data))
m.PMTGSS7 = data.frame(model.matrix(~PMTGSS7-1, m@data))
m.DRNGSS7 = data.frame(model.matrix(~DRNGSS7-1, m@data))
#mRF$xvar.names[which(!mRF$xvar.names %in% names(cbind(m@data, m.PVEGKT6, m.LNDCOV6, m.PMTGSS7, m.DRNGSS7)))]
m@data = cbind(m@data, m.PVEGKT6, m.LNDCOV6, m.PMTGSS7, m.DRNGSS7)[,mRF$xvar.names]
soiltype.RF = predict.rfsrc(mRF, m@data, na.action="na.impute", importance=FALSE, membership=TRUE)
m1 = m[1]
m1@data = data.frame(soiltype.RF$predicted)*100
plot(stack(m1[30:40]))
for(j in 1:ncol(m1)){
  out <- paste0("/data/NASIS/predicted100m/T1049/TAXgg_", names(m1)[j], "_T1049.tif")
  if(!file.exists(out)){
    writeGDAL(m1[j], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

## FIT MODELS FOR NUMERIC VARIABLES:


