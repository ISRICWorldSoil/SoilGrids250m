## Fit models for GreatGroups based on NASIS points (ca 350,000)
## By Tom.Hengl@isric.org points prepare by Travis Neuman

setwd("/data/NASIS")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
library(caret)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger", ref="forest_memory") ## version to deal with Memory problems
library(ranger)
library(nnet)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(plotKML)
library(GSIF)
library(rgeos)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("/data/models/extract.equi7t3.R")
source("/data/models/saveRDS_functions.R")
source("/data/models/wrapper.predict_cs.R")

load("/data/models/equi7t3.rda")
load("/data/models/equi7t1.rda")
des <- read.csv("/data/models/SoilGrids250m_COVS250m.csv")

## Great groups
unzip("NASIS_L48_gg.zip")
NASISgg.pnts <- readOGR("nasispts16_gg_L48.shp", "nasispts16_gg_L48")
## 304,454 points
summary(NASISgg.pnts$gg)
## Remove smaller classes:
NASISgg.pnts$soiltype <- NA
for(i in levels(NASISgg.pnts$gg)){
  sel <- grep(pattern=i, NASISgg.pnts$gg)
  if(length(sel)>5){
    NASISgg.pnts$soiltype[sel] <- i
  }
}
str(summary(NASISgg.pnts$gg))

## soil texture data:
unzip("nasispts16_pscs_L48.zip")
NASISpscs.pnts <- readOGR("nasispts16_pscs_L48.shp", "nasispts16_pscs_L48")
## 299,487 points
str(NASISpscs.pnts@data)
summary(NASISpscs.pnts$pscs)
NASISpscs.pnts$textype <- NA
for(i in levels(NASISpscs.pnts$pscs)){
  sel <- grep(pattern=i, NASISpscs.pnts$pscs)
  if(length(sel)>5){
    NASISpscs.pnts$textype[sel] <- i
  }
}
NASISpscs.pnts$textype <- as.factor(NASISpscs.pnts$textype)
summary(NASISpscs.pnts$textype)

## OVERLAY AND FIT MODELS:
ov <- extract.equi7(x=NASISgg.pnts, y=des$WORLDGRIDS_CODE, equi7=equi7t3, path="/data/covs", cpus=48)
write.csv(ov, file="ov.NASISgg_SoilGrids250m.csv")
unlink("ov.NASISgg_SoilGrids250m.csv.gz")
gzip("ov.NASISgg_SoilGrids250m.csv")
save(ov, file="ov.NASISgg.rda")
summary(ov$soiltype)
#load("ov.NASISgg.rda")
ov2 <- extract.equi7(x=NASISpscs.pnts, y=des$WORLDGRIDS_CODE, equi7=equi7t3, path="/data/covs", cpus=48)
save(ov2, file="ov.NASISpscs.rda")

## ------------- MODEL FITTING -----------

pr.lst <- des$WORLDGRIDS_CODE
formulaString.USDA = as.formula(paste('soiltype ~ ', paste(pr.lst, collapse="+")))
#formulaString.USDA
ovA <- ov[complete.cases(ov[,all.vars(formulaString.USDA)]),]
ovA$soiltype <- as.factor(paste(ovA$soiltype))

formulaString.pscs = as.formula(paste('textype ~ ', paste(pr.lst, collapse="+")))
ovA2 <- ov2[complete.cases(ov2[,all.vars(formulaString.pscs)]),]
str(ovA2)

## Ranger package (https://github.com/imbs-hl/ranger)
mrfX_NASISgg <- ranger::ranger(formulaString.USDA, ovA, importance="impurity", write.forest=TRUE, probability=TRUE)
mrfX_NASISpscs <- ranger::ranger(formulaString.pscs, ovA2, importance="impurity", write.forest=TRUE, probability=TRUE)

cat("Results of model fitting 'randomForest':\n", file="NASIS_resultsFit.txt")
cat("\n", file="NASIS_resultsFit.txt", append=TRUE)
cat(paste("Variable: USDA Great Group"), file="NASIS_resultsFit.txt", append=TRUE)
cat("\n", file="NASIS_resultsFit.txt", append=TRUE)
sink(file="NASIS_resultsFit.txt", append=TRUE, type="output")
cat("\n Random forest model:", file="NASIS_resultsFit.txt", append=TRUE)
print(mrfX_NASISgg)
cat("\n Variable importance:\n", file="NASIS_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_NASISgg))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
cat("\n", file="NASIS_resultsFit.txt", append=TRUE)
cat(paste("Variable: SPCS classes"), file="NASIS_resultsFit.txt", append=TRUE)
cat("\n", file="NASIS_resultsFit.txt", append=TRUE)
cat("\n Random forest model:", file="NASIS_resultsFit.txt", append=TRUE)
print(mrfX_NASISpscs)
cat("\n Variable importance:\n", file="NASIS_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_NASISpscs))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
sink()

## save objects in parallel:
saveRDS.gz(mrfX_NASISgg, file="mrfX_NASISgg.rds")
saveRDS.gz(mrfX_NASISpscs, file="mrfX_NASISpscs.rds")
save.image()

## ------------- PREDICTIONS -----------

## Predict for the whole of USA:
library(maps)
library(maptools)
usa.m <- map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(usa.m$names, ":"), function(x) x[1])
prj = "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
state = map2SpatialPolygons(usa.m, IDs=IDs)
proj4string(state) = "+proj=longlat +datum=WGS84"
state <- spTransform(state, CRS(proj4string(equi7t1[["NA"]])))
ov.state <- over(y=state, x=equi7t1[["NA"]])
#ov.state <- gIntersection(state, equi7t1[["NA"]], byid = TRUE)
#str(ov.state@data)

new.dirs = unique(paste0("NA_", equi7t1[["NA"]]$TILE[which(!is.na(ov.state))])) #levels(as.factor(ov$equi7))
## 906 tiles
x <- lapply(paste0("./", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## Split models otherwise too large in size:
num_splits=20
#mrfX_NASISgg = readRDS.gz("mrfX_NASISgg.rds")
mrfX_NASISgg_final <- split_rf(mrfX_NASISgg, num_splits)
for(j in 1:length(mrfX_NASISgg_final)){
  gm = mrfX_NASISgg_final[[j]]
  saveRDS.gz(gm, file=paste0("mrfX_NASISgg_", j,".rds"))
}
rm(mrfX_NASISgg); rm(mrfX_NASISgg_final)
gc(); gc()
save.image()

mrfX_NASISpscs = readRDS.gz("mrfX_NASISpscs.rds")
mrfX_NASISpscs_final <- split_rf(mrfX_NASISpscs, num_splits)
for(j in 1:length(mrfX_NASISpscs_final)){
  gm = mrfX_NASISpscs_final[[j]]
  saveRDS.gz(gm, file=paste0("mrfX_NASISpscs_", j,".rds"))
}
rm(mrfX_NASISpscs); rm(mrfX_NASISpscs_final)
gc(); gc()
save.image()

#del.lst <- list.files(path="/data/NASIS", pattern=glob2rx("^TAXgg_*_*_*_*.rds$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

#model.n = "mrfX_NASISgg_"
#varn = "TAXgg"
model.n = "mrfX_NASISpscs_"
varn = "PSCS"
out.path = "/data/NASIS"
for(j in 1:num_splits){
  gm = readRDS.gz(paste0(model.n, j,".rds"))
  cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
  sfInit(parallel=TRUE, cpus=ifelse(cpus>46, 46, cpus))
  sfExport("gm", "new.dirs", "split_predict_c", "j", "varn", "out.path")
  sfLibrary(ranger)
  x <- sfLapply(new.dirs, fun=function(x){ if(length(list.files(path = paste0(out.path, "/", x, "/"), glob2rx("*.rds$")))<j){ try( split_predict_c(x, gm, in.path="/data/covs1t", out.path=out.path, split_no=j, varn=varn) ) } } )  ## , num.threads=5
  sfStop()
  rm(gm)
}

## Test it:
#sum_predict_ranger(i="NA_060_036", in.path="/data/covs1t", out.path="/data/NASIS", varn="TAXgg", num_splits)
sum_predict_ranger(i="NA_060_036", in.path="/data/covs1t", out.path="/data/NASIS", varn="PSCS", num_splits)

## Sum up predictions:
sfInit(parallel=TRUE, cpus=30)
sfExport("new.dirs", "sum_predict_ranger", "num_splits", "varn")
sfLibrary(rgdal)
sfLibrary(plyr)
x <- sfLapply(new.dirs, fun=function(x){ try( sum_predict_ranger(x, in.path="/data/covs1t", out.path="/data/NASIS", varn=varn, num_splits) )  } )
sfStop()


## Create mosaics:
mosaic_tiles_NASIS <- function(j, in.path, varn, tr=0.002083333, r="bilinear", ot="Byte", dstnodata=255, out.path){
  out.tif <- paste0(out.path, varn, "_", j, '_250m_ll.tif')
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"'))
  }
}

levs = list.files(path="./NA_060_036", pattern=glob2rx(paste0("^",varn,"_*_*_*_*.tif$")))
levs = sapply(basename(levs), function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_NASIS", "varn")
out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_NASIS(x, in.path="/data/NASIS/", varn=varn, out.path="/data/NASIS/") )})
sfStop()

save.image()
