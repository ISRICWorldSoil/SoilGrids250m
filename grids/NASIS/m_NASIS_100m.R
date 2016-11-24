## Fit models for USDA GreatGroups, texture classes, and several soil properties using the NASIS points (ca 350,000 with classification and 45,000 points with soil properties)
## Code by Tom.Hengl@isric.org / points prepared by Amanda Rachmaran (a.m.ramcharan@gmail.com) and Travis Nauman (tnauman@usgs.gov)

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "ROCR", "randomForest", "R.utils", "plyr", "parallel", "psych", "mda", "dismo", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "mxnet", "doParallel", "caret", "plotKML")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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
#library(randomForestSRC)
#options(rf.cores=detectCores(), mc.cores=detectCores())
library(parallel)
library(doParallel)
library(mxnet)

source("/data/models/saveRDS_functions.R")
source("/data/models/wrapper.predict_cs.R")
source("functions_NASIS.R")
levsGG <- read.csv("TAXOUSDA_GreatGroups.csv")

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

## Grid definition
r <- raster("/data/USA48/elev48i0100a.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
NODATA_value = -32768
te = paste(extent(r)[c(1,3,2,4)], collapse = " ")

## check if all layers are ready:
des <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_Covs100m.csv")
s = raster::stack(paste0('/data/GEOG/NASIS/covs100m/', des$WORLDGRIDS_CODE, ".tif"))
#str(s[1])

## legends
des <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_Covs100m.csv")
LNDCOV6.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_LandCover.csv")
PMTGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_pmaterial.csv")
DRNGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_drainage.csv")
PVEGKT6.leg <- read.csv("/data/GEOG/NASIS/covs100m/potential_vegetation_legend.csv")
COUNTY6.leg <- read.csv("/data/GEOG/NASIS/covs100m/Counties_legend.csv")
#GEOUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/geology_legend.csv")
GESUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/surfacegeology_legend.csv")

## Training points (taxonomy, texture classes):
TAX_gg.pnts = readRDS("TAX_gg.pnts.rds")
str(TAX_gg.pnts@data)
NASISpscs.pnts = readRDS("NASISpscs.pnts.rds")
NASISpscs.pnts$LOC_ID = paste("ID", NASISpscs.pnts$xwgs84, NASISpscs.pnts$ywgs84, sep="_")

## Soil properties NRCS NCSS:
#system("7za e Geo_Peds07272016.csv.gz")
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
obj <- GDALinfo("/data/GEOG/NASIS/covs100m/COUNTY6.tif")
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
system(paste('gdal_translate /data/GEOG/NASIS/covs100m/COUNTY6.tif /data/USA48/COUNTY6.sdat -ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
system(paste0(saga_cmd, ' -c=48 shapes_grid 2 -GRIDS=\"/data/USA48/COUNTY6.sgrd\" -POLYGONS=\"tiles50km.shp\" -PARALLELIZED=1 -RESULT=\"ov_COUNTY6_tiles50km.shp\"'))
ov_COUNTY6 = readOGR("ov_COUNTY6_tiles50km.shp", "ov_COUNTY6_tiles50km")
summary(sel.t <- !ov_COUNTY6@data[,"COUNTY6..MA"]==-32767)
## 3392 tiles with values
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

## Clean up:
#rds.lst = list.files(path="./covs100m", pattern=glob2rx("*.rds"), recursive = TRUE, full.names = TRUE)
#unlink(rds.lst)
## Test:
#make_RDS_tiles(i=1922, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)

## TAKES ca 16 HOURS TO PREPARE
sfInit(parallel=TRUE, cpus=48) ## 25
sfExport("make_RDS_tiles", "tile.tbl", "covs.lst", "LNDCOV6.leg", "PMTGSS7.leg", "DRNGSS7.leg", "PVEGKT6.leg", "COUNTY6.leg", "GESUSG6.leg")
sfLibrary(plyr)
sfLibrary(rgdal)
#sfLibrary(randomForestSRC)
sfLibrary(ranger)
sfLibrary(sp)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( make_RDS_tiles(i, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg) ) } )
sfStop()
save.image()

## Clean up empty dirs
pr.dirs <- basename(dirname(list.files(path="/data/NASIS/covs100m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
## 3397
pr.dirs.c <- list.dirs("/data/NASIS/covs100m")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
#x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})

## Gap-filled Parent material and drainage maps:
for(j in c("PMTGSS7","DRNGSS7")){
  if(!file.exists(paste0('/data/GEOG/NASIS/covs100m/', j, '_f.tif'))){
    pm.tif = list.files(path="/data/NASIS/covs100m", pattern=paste0(j, "_predicted"), recursive = TRUE, full.names = TRUE)
    cat(pm.tif, sep="\n", file="my_liste.txt")
    system(paste0(gdalbuildvrt, ' -input_file_list my_liste.txt ', j, '_f.vrt'))
    system(paste0(gdalwarp, ' ', j, '_f.vrt /data/GEOG/NASIS/covs100m/', j, '_f.tif -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000 -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
  }
}
save.image()

## OVERLAY (takes ca 15 mins):
ov <- extract.tiled(x=pnts, tile.pol=tile.pol, path="/data/NASIS/covs100m", ID="ID", cpus=48)
ov$PMTGSS7 = as.factor(ov$PMTGSS7)
ov$DRNGSS7 = as.factor(ov$DRNGSS7)
str(ov)
## "PMTGSS7" has too many levels for RF models
## http://stackoverflow.com/questions/14921805/convert-a-factor-to-indicator-variables
PMTGSS7.mat = data.frame(model.matrix(~PMTGSS7-1, ov))
xsum.PMTGSS7 = sapply(PMTGSS7.mat, sum, na.rm=TRUE)
## 10 categories have <10 observations
sel.rm.PMTGSS7 = names(PMTGSS7.mat)[which(xsum.PMTGSS7<6)]
## Land cover:
summary(ov$LNDCOV6)
LNDCOV6.mat = data.frame(model.matrix(~LNDCOV6-1, ov))
xsum.LNDCOV6 = sapply(LNDCOV6.mat, sum, na.rm=TRUE)
sel.rm.LNDCOV6 = names(LNDCOV6.mat)[which(xsum.LNDCOV6<6)]
## Potential vegetation:
PVEGKT6.mat = data.frame(model.matrix(~PVEGKT6-1, ov))
xsum.PVEGKT6 = sapply(PVEGKT6.mat, sum, na.rm=TRUE)
sel.rm.PVEGKT6 = names(PVEGKT6.mat)[which(xsum.PVEGKT6<6)]
## Drainage class:
summary(ov$DRNGSS7)
DRNGSS7.mat = data.frame(model.matrix(~DRNGSS7-1, ov))
xsum.DRNGSS7 = sapply(DRNGSS7.mat, sum, na.rm=TRUE)
sel.rm.DRNGSS7 = names(DRNGSS7.mat)[which(xsum.DRNGSS7<6)]

## Merge all indicators with other covariates
ov.fs = plyr::join_all(list(ov, 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PMTGSS7.mat))]), PMTGSS7.mat[,names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(LNDCOV6.mat))]), LNDCOV6.mat[,names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PVEGKT6.mat))]), PVEGKT6.mat[,names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(DRNGSS7.mat))]), DRNGSS7.mat[,names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)]])))
## remove layers with too many missing values?
c.check = sapply(ov.fs, function(i){sum(is.na(i))})
names(ov.fs)[c.check>4000]
ov.fs$WETGSS7 = ifelse(is.na(ov.fs$WETGSS7), 0, ov.fs$WETGSS7)
summary(ov.fs$WETGSS7)

## join sampled values and covariates
ovA_NCSS_peds = plyr::join(NCSS_peds, ov.fs, by="LOC_ID", type="left")
ovA_TAX_gg.pnts = plyr::join(as.data.frame(TAX_gg.pnts), ov.fs, by="LOC_ID", type="left")
ovA_NASISpscs.pnts = plyr::join(as.data.frame(NASISpscs.pnts), ov.fs, by="LOC_ID", type="left")
## Regression matrix per variable:
saveRDS.gz(ovA_NCSS_peds, file="ovA_NCSS_peds.rds")
saveRDS.gz(ovA_TAX_gg.pnts, file="ovA_TAX_gg.pnts.rds")
saveRDS.gz(ovA_NASISpscs.pnts, file="ovA_NASISpscs.pnts.rds")
#write.csv(ovA_TAX_gg.pnts[1:200,], "ovA_TAX_gg_pnts_example.csv")

## covariates - USE ONLY INDICATORS (NO FACTOR-TYPE predictors)
## See also: 'respect.unordered.factors=TRUE' https://www.r-bloggers.com/on-ranger-respect-unordered-factors/
pr.lst <- des$WORLDGRIDS_CODE
## tidy up list of covariates for soil-class mapping
prC.lst = c(paste(pr.lst[-c(grep(pattern = "GESUSG6",pr.lst), grep(pattern = "PMTGSS7",pr.lst), grep(pattern = "PVEGKT6",pr.lst), grep(pattern = "LNDCOV6",pr.lst), grep(pattern = "DRNGSS7",pr.lst), grep(pattern = "COUNTY6",pr.lst), grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst), grep(pattern = glob2rx("N??PRI5"),pr.lst), grep(pattern = glob2rx("V??PRI5"),pr.lst))]), names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)], names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)], names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)], names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)])
str(prC.lst)
## 214 covs
prT.lst = c(prC.lst, paste(pr.lst[c(grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst))]))
str(prT.lst)
## '"DRNGSS7NULL"'? 

## Clean up:
rds.lst = list.files(path="./covs100m", pattern=glob2rx("cT*.rds"), recursive = TRUE, full.names = TRUE)
unlink(rds.lst)
#make_newdata(i=1922, in.path="/data/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prT.lst)
## Make prediction locations in parallel (TAKES 1.5 HRS):
sfInit(parallel=TRUE, cpus=48)
sfExport("t.sel", "make_newdata", "PMTGSS7.leg", "DRNGSS7.leg", "prT.lst")
sfLibrary(plyr)
sfLibrary(sp)
x <- sfClusterApplyLB(as.numeric(t.sel), fun=function(i){ try( make_newdata(i, in.path="/data/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prT.lst) ) } ) 
sfStop()

#############################################
## FIT MODELS FOR FACTORS
#############################################

t.fvars = c("soiltype", "textype")
t.names = list("USDA Great Group", "SPCS classes")
names(t.names) = t.fvars
save.image()
unlink(paste0("mRF_", t.fvars,".rds"))

## fit models for factor variables
## TAKES ca 40 mins per variable and >190GB RAM
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
Nsub = 8e3
for(j in t.fvars){
  if(!file.exists(paste0("mRF_", j,".rds"))){
    formulaString = as.formula(paste(j, ' ~ ', paste(prC.lst, collapse="+")))
    if(j=="soiltype"){ x = ovA_TAX_gg.pnts[,all.vars(formulaString)]  }
    if(j=="textype"){ x = ovA_NASISpscs.pnts[,all.vars(formulaString)]  }
    x = x[complete.cases(x),] ## 325,566 for 'soiltype'
    x[,j] = droplevels(x[,j])
    rf.tuneGrid <- expand.grid(mtry = seq(5,round(sqrt(ncol(x))*2),by=5))
    if(!file.exists(paste0("t.mRF_", j,".rds"))){
      ## estimate Mtry using cross-validation (subsample original RM)
      xs = x[sample.int(nrow(x), Nsub),]
      xs[,j] = droplevels(xs[,j])
      cl <- makeCluster(48)
      registerDoParallel(cl)
      t.mrfX <- caret::train(formulaString, data=xs, method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid, num.trees=85)
      stopCluster(cl)    
      saveRDS.gz(t.mrfX, paste0("t.mRF_", j,".rds"))
    } else {
      t.mrfX = readRDS.gz(paste0("t.mRF_", j,".rds"))
    }
    ## How Many Trees in a Random Forest? 65-120 is more than enough (http://dx.doi.org/10.1007/978-3-642-31537-4_13)
    gc(); gc()
    mRF = ranger::ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, probability=TRUE, num.trees=85, mtry=t.mrfX$bestTune$mtry) ## optimize based on caret
    saveRDS.gz(mRF, file=paste0("mRF_", j,".rds"))
    cat("Results of model fitting 'randomForest':\n", file=paste0("NASIS_resultsFit_", j, ".txt"))
    cat("\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    cat(paste("Variable:", t.names[[j]]), file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    cat("\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    sink(file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE, type="output")
    cat("\n Random forest model:", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    print(mRF)
    cat("\n Variable importance:\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    xl <- as.list(ranger::importance(mRF))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:150]])))
    sink()
  }
}
closeAllConnections()
gc(); gc()
#object.size(mRF)/1e6

#############################################
## PREDICTIONS - CLASSES
#############################################

# del.lst <- list.files(path="/data/NASIS/predicted100m", pattern=glob2rx("^TAXgg*.rds$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst) 
# del.lst <- list.files(path="/data/NASIS/predicted100m", pattern=glob2rx("^PSCS*.tif$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)
# del.lst <- list.files(path="/data/NASIS/predicted100m", pattern=glob2rx("^TAXgg*.tif$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)

## Predict classes in parallel:
#mRF = readRDS.gz("mRF_soiltype.rds")
mRF = readRDS.gz("mRF_textype.rds")
#str(mRF)
levs = attr(mRF$predictions, "dimnames")[[2]]
## 291 class for soiltype
## 65 for textype
#varn = "TAXgg"
varn = "PSCS"

out.path = "/data/NASIS/predicted100m"
## Test predictions:
#make_newdata(i=1923, in.path="/data/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prC.lst)
#predict_factor_tile(i=1922, mRF, varn, levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m")

## run predictions in sequence:
#for(i in as.numeric(t.sel)){ if(!file.exists(paste0("/data/NASIS/predicted100m/T", i, "/", varn, "_T", i, ".tif"))){ try( predict_factor_tile(i, mRF, varn, levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m") ) } }

## run predictions in parallel 30 HRS = 3600 tiles * 1 min / 2 tiles / 60 mins:
n.cores = round(200/(object.size(mRF)/1e9))
#cl <- makeCluster(ifelse(n.cores>48, 48, n.cores), type="FORK")
cl <- makeCluster(30, type="FORK")
x = parLapply(cl, as.numeric(t.sel), fun=function(i){ if(!file.exists(paste0("/data/NASIS/predicted100m/T", i, "/", varn, "_T", i, ".tif"))){ try( predict_factor_tile(i, mRF, varn, levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", mc.cores=1) ) } } )
stopCluster(cl)
gc(); gc()

## mosaic tiles:
levs = list.files(path="./predicted100m/T5445", pattern=glob2rx(paste0("^",varn,"_*_T*.tif$")))
levs = sapply(basename(levs), function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_100m", "varn", "te")
out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_100m(x, in.path="/data/NASIS/predicted100m", varn=varn, te=te) )})
sfStop()

## Clean up emty dirs:
pr.dirs <- basename(dirname(list.files(path="/data/NASIS/predicted100m", pattern=glob2rx("PSCS_T*.tif$"), recursive = TRUE, full.names = TRUE)))
## 3392 - 3379 = 13 empty dirs
pr.dirs.c = list.dirs("/data/NASIS/predicted100m")[-1]
pr.dirsT <- list.dirs("/data/NASIS/covs100m")[-1]
selD <- which(!basename(pr.dirs.c) %in% basename(pr.dirsT))
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})
save.image()

#############################################
## FIT MODELS FOR NUMERIC VARIABLES
#############################################

str(ovA_NCSS_peds)
xs = lapply(t.props, function(i){summary(ovA_NCSS_peds[,i])})
names(xs) = t.props
ovA_NCSS_peds$DEPTH = ovA_NCSS_peds$hzn_top + (ovA_NCSS_peds$hzn_bot - ovA_NCSS_peds$hzn_top)/2
## 7 standard depths
breaks = c(0, rowMeans(data.frame(c(0,5,15,30,60,100), c(5,15,30,60,100,200))), 200, 4500)
ovA_NCSS_peds$DEPTH_c = cut(ovA_NCSS_peds$DEPTH, breaks, labels = paste0("sl", 1:8))
summary(ovA_NCSS_peds$DEPTH_c)
hist(ovA_NCSS_peds$DEPTH) ## some profiles are up to 35m deep?
## Estimate SG variables per depth:
sg.props = c("CLYPPT", "SNDPPT", "ORCDRC", "BLDFIE", "")
for(j in 1:length(t.props)){
  if(!t.props[j]=="n_tot"){
    ovA_NCSS_peds[,paste0("sg_", t.props[j])] = ifelse(ovA_NCSS_peds$DEPTH_c=="sl1", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl1_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl2", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl2_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl3", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl3_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl4", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl4_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl5", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl5_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl6", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl6_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl7", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl7_100m")], NA)))))))
  }
}
#xyplot(clay~sg_clay, ovA_NCSS_peds)

## fit models in a loop
## Generic settings for caret:
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
Nsub = 8e3
rf.tuneGrid <- expand.grid(mtry = seq(10,60,by=5))

cl <- makeCluster(48)
registerDoParallel(cl)
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="SPROPS_resultsFit.txt")
for(j in t.props){
  out.rf <- paste0("mRF.", j,".rds")
  if(!file.exists(out.rf)){
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    cat(paste("Variable:", j, file="SPROPS_resultsFit.txt", append=TRUE))
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    if(j == "n_tot"){
      formulaString = as.formula(paste(j, ' ~ DEPTH + ', paste(prC.lst, collapse="+")))
    } else {
      formulaString = as.formula(paste0(j, ' ~ DEPTH + ', paste(prC.lst, collapse="+"), "+ sg_", j))
    }
    x = ovA_NCSS_peds[,all.vars(formulaString)]
    x = x[complete.cases(x),]
    ## Caret training settings (reduce number of combinations to speed up):
    ## optimize mtry parameter:
    if(!file.exists(gsub("mRF","t.mRF",out.rf))){
      t.mrfX <- caret::train(formulaString, data=x[sample.int(nrow(x), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
      saveRDS.gz(t.mrfX, file=gsub("mRF","t.mRF",out.rf))
    } else {
      t.mrfX <- readRDS.gz(gsub("mRF","t.mRF",out.rf))
    }
    ## fit RF model using 'ranger' (fully parallelized)
    ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
    mrfX <- ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
    saveRDS.gz(mrfX, file=out.rf)
    ## Top covariates:
    sink(file="SPROPS_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
    xl <- as.list(ranger::importance(mrfX))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:150]])))
    ## XGBoost
    if(!file.exists(paste0("mGB.",j,".rds"))){
      mGBX <- caret::train(formulaString, data=x, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid)
      saveRDS.gz(mGBX, file=paste0("mGB.",j,".rds"))
      ## save also binary model for prediction purposes:
      xgb.save(mGBX$finalModel, paste0("Xgb.",j))
    } else {
      mGBX <- readRDS.gz(paste0("mGB.",j,".rds"))
    }
    importance_matrix <- xgb.importance(mGBX$coefnames, model = mGBX$finalModel)
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(mGBX)
    cat("\n XGBoost variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:150,])
    cat("--------------------------------------\n", file="SPROPS_resultsFit.txt", append=TRUE)
    sink()
  }
}
rm(mrfX); rm(mGBX)
stopCluster(cl); closeAllConnections()

#############################################
## PREDICTIONS - NUMERIC VARIABLES
#############################################

z.min <- as.list(c(0,0,0,50,0))
names(z.min) = t.props
z.max <- as.list(c(100,100,710,3500,580))
names(z.max) = t.props
type.lst <- c("Byte","Byte","Int16","Int16","Int16")
names(type.lst) = t.props
mvFlag.lst <- c(255, 255, -32768, -32768, -32768)
names(mvFlag.lst) = t.props
sg.var <- c("CLYPPT", "SNDPPT", "ORCDRC", "BLDFIE", "")
names(sg.var) = t.props

## Run per property
## TAKES ABOUT 10hrs for RF models and about 2hrs for XGB models
for(j in t.props){
  ## it appears that 'parallel' package conflicts with snowfall package so needs to be turned off:
  detach("package:snow", unload=TRUE)
  detach("package:snowfall", unload=TRUE)
  ##``Error in UseMethod("sendData") : 
  ##  no applicable method for 'sendData' applied to an object of class "SOCK0node"''
  if(j=="pH"|j=="soc"){ multiplier = 10 }
  if(j %in% c("n_tot")){ multiplier = 100 }
  if(j %in% c("bd")){ multiplier = 1000 }
  if(j %in% c("clay","sand")){ multiplier = 1 }
  ## Random forest predictions:
  gm = readRDS.gz(paste0("mRF.", j,".rds"))
  gm1.w = 1/gm$prediction.error
  cpus = unclass(round((256-50)/(3.5*(object.size(gm)/1e9))))
  cl <- makeCluster(ifelse(cpus>35, 35, cpus), type="FORK")
  if(j == "n_tot"){
    x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/NASIS/predicted100m/", x, "/", j,"_", x, "_rf.rds"))){ try( split_predict_n(x, gm, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", split_no=NULL, varn=j, method="ranger", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0(in.path, "/", x, "/c", x,".rds"), SG.col=NULL ) ) } } )
  } else {
    x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/NASIS/predicted100m/", x, "/", j,"_", x, "_rf.rds"))){ try( split_predict_n(x, gm, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", split_no=NULL, varn=j, method="ranger", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0(in.path, "/", x, "/c", x,".rds"), SG.col=paste0(sg.var[j], "_M_sl", 1:7, "_100m") ) ) } } )  
  }
  stopCluster(cl)
  gc(); gc()
  ## XGBoost:
  gm = readRDS.gz(paste0("mGB.", j,".rds"))
  gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
  cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
  cl <- makeCluster(ifelse(cpus>35, 35, cpus), type="FORK")
  if(j == "n_tot"){
    x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/NASIS/predicted100m/", x, "/", j,"_", x, "_xgb.rds"))){ try( split_predict_n(x, gm, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", split_no=NULL, varn=j, method="xgboost", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0(in.path, "/", x, "/c", x,".rds"), SG.col=NULL ) ) } } )
  } else {
    x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/NASIS/predicted100m/", x, "/", j,"_", x, "_xgb.rds"))){ try( split_predict_n(x, gm, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", split_no=NULL, varn=j, method="xgboost", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0(in.path, "/", x, "/c", x,".rds"), SG.col=paste0(sg.var[j], "_M_sl", 1:7, "_100m") ) ) } } )  
  }
  stopCluster(cl)
  rm(gm)
  gc(); gc()
  ## sum up predictions:
  library(snowfall)
  sfInit(parallel=TRUE, cpus=45)
  sfExport("t.sel", "sum_predict_ensemble", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(paste0("T", t.sel), fun=function(x){ try( if(length(list.files(path = paste0("/data/NASIS/predicted100m/", x, "/"), glob2rx(paste0("^", j, "_M_sl*_", x, ".tif$"))))==0){ try( sum_predict_ensemble(x, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", varn=j, num_splits=NULL, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]]) ) } )  } )
  sfStop()
}

## corrupt or missing tiles:
sfInit(parallel=TRUE, cpus=length(t.vars))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("missing.tiles", "t.vars", "pr.dirs")
missing.lst <- sfLapply(t.vars, missing.tiles, pr.dirs=pr.dirs)
sfStop()
names(missing.lst) = t.vars

## mosaic tiles:
levs = list.files(path="./predicted100m/T5445", pattern=glob2rx(paste0("^",varn,"_*_T*.tif$")))
levs = sapply(basename(levs), function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_100m", "varn", "te")
out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_100m(x, in.path="/data/NASIS/predicted100m", varn=varn, te=te) )})
sfStop()
