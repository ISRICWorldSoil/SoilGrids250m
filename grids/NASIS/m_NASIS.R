## Fit models for GreatGroups based on NASIS points (ca 350,000)
## Code by Tom.Hengl@isric.org / points prepared by Travis Nauman (tnauman@usgs.gov)

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

## Great groups (points prepared by Travis)
unzip("NASIS_L48_gg.zip", junkpaths=TRUE)
NASISgg.pnts <- readOGR("nasispts16_gg_L48.shp", "nasispts16_gg_L48")
## 310,690 points
str(NASISgg.pnts@data)
#summary(NASISgg.pnts$gg)
NASISgg.xy <- plyr::rename(NASISgg.pnts@data[,c("upedonid","xwgs84","ywgs84","grtgrp_nm")], replace=c("upedonid"="SOURCEID", "xwgs84"="LONWGS84", "ywgs84"="LATWGS84", "grtgrp_nm"="TAXNUSDA"))
NASISgg.xy$SOURCEDB = "NASIS"
#plot(NASISgg.xy[,2:3], pch=".")

## NCSS points (prepared by Amanda and Dylan):
load("PedsT05162016.Rda")
Peds.tax <- Peds[!duplicated(Peds$pedon_key),c("pedon_key","long","lat","tax_grtgrp")]
## 34,183 points
Peds.tax <- plyr::rename(Peds.tax, replace=c("pedon_key"="SOURCEID", "long"="LONWGS84", "lat"="LATWGS84", "tax_grtgrp"="TAXNUSDA"))
Peds.tax$SOURCEDB = "NCSS"

## Global compilation of soil classes (various sources):
TAXOUSDA_global <- readRDS("TAXNUSDA.global.rds")
TAXOUSDA_global <- as.data.frame(TAXOUSDA_global[!TAXOUSDA_global$SOURCEDB=="NCSS",])[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNUSDA","TAXOUSDA")]
## 25,933 points

## Combine the three data sets:
all.pnts <- plyr::rbind.fill(NASISgg.xy, Peds.tax, TAXOUSDA_global)
all.pnts$LOC_ID <- paste("ID", all.pnts$LONWGS84, all.pnts$LATWGS84, sep="_") 
summary(duplicated(all.pnts$LOC_ID))
## Select NASIS points and then all other points but which are not duplicates:
selR = all.pnts$LOC_ID[all.pnts$SOURCEDB=="NASIS"]
sel.pnts <- all.pnts[!all.pnts$SOURCEDB=="NASIS",]
## 60,116
all.pnts <- rbind(all.pnts[all.pnts$SOURCEDB=="NASIS",], sel.pnts[which(!sel.pnts$LOC_ID %in% selR),])
## 351,375

## Great_Group names --> we try to locate them in the raw names one by one
levsGG <- read.csv("TAXOUSDA_GreatGroups.csv")
#summary(levsGG$Great_Group)
levs = levels(levsGG$Great_Group)
## 434 levels
strip_s <- function(x){ ifelse(substr(x, nchar(x), nchar(x))=="s", substr(x, 1, nchar(x)-1), x) }
taxn.lst <- list(NULL)
## takes 2 mins:
for(j in 1:length(levs)){
  ## remove "s" if at the end of the class name:
  pat <- strip_s(levs[j]) 
  sel1 <- grep(pat, paste(all.pnts$TAXNUSDA), ignore.case=TRUE) # , fixed=TRUE
  sel2 <- grep(pat, paste(all.pnts$TAXOUSDA), ignore.case=TRUE) # , fixed=TRUE
  sel <- unique(c(sel1, sel2))
  ## bind together
  if(length(sel)>0){
    taxn.lst[[j]] <- all.pnts[sel,]
    taxn.lst[[j]]$gg <- levs[j]
  }
}
TAX_gg.pnts <- do.call(rbind, taxn.lst)
rm(taxn.lst)
TAX_gg.pnts$gg <- as.factor(TAX_gg.pnts$gg)

## Remove smaller classes:
xg = summary(TAX_gg.pnts$gg, maxsum=length(levels(TAX_gg.pnts$gg)))
selg.levs = attr(xg, "names")[xg > 5]
TAX_gg.pnts$soiltype <- as.factor(TAX_gg.pnts$gg)
TAX_gg.pnts$soiltype[which(!TAX_gg.pnts$gg %in% selg.levs)] <- NA
TAX_gg.pnts$soiltype <- droplevels(TAX_gg.pnts$soiltype)
str(summary(TAX_gg.pnts$soiltype, maxsum=length(levels(TAX_gg.pnts$soiltype))))
## On the end total of: 318 classes
coordinates(TAX_gg.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAX_gg.pnts) = NASISgg.pnts@proj4string

## soil texture data:
unzip("NASIS_L48_pscs.zip", junkpaths=TRUE)
NASISpscs.pnts <- readOGR("nasispts16_pscs_L48.shp", "nasispts16_pscs_L48")
## 306,583 points
str(NASISpscs.pnts@data)
xs = summary(NASISpscs.pnts$pscs, maxsum=length(levels(NASISpscs.pnts$pscs)))
sel.levs = attr(xs, "names")[xs > 5]
NASISpscs.pnts$textype <- NASISpscs.pnts$pscs
NASISpscs.pnts$textype[which(!NASISpscs.pnts$pscs %in% sel.levs)] <- NA
NASISpscs.pnts$textype <- droplevels(NASISpscs.pnts$textype)
save.image()


## OVERLAY AND FIT MODELS:
ov <- extract.equi7(x=TAX_gg.pnts, y=des$WORLDGRIDS_CODE, equi7=equi7t3, path="/data/covs", cpus=48)
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

formulaString.pscs = as.formula(paste('textype ~ ', paste(pr.lst, collapse="+")))
ovA2 <- ov2[complete.cases(ov2[,all.vars(formulaString.pscs)]),]
str(ovA2)

## Ranger package (https://github.com/imbs-hl/ranger)
mrfX_NASISgg <- ranger::ranger(formulaString.USDA, ovA, importance="impurity", write.forest=TRUE, probability=TRUE)
mrfX_NASISpscs <- ranger::ranger(formulaString.pscs, ovA2, importance="impurity", write.forest=TRUE, probability=TRUE)
gc()
gc()

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
#library(maps)
#library(maptools)
#usa.m <- map('state', plot=FALSE, fill=TRUE)
#IDs <- sapply(strsplit(usa.m$names, ":"), function(x) x[1])
#prj = "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#state = map2SpatialPolygons(usa.m, IDs=IDs)
#proj4string(state) = "+proj=longlat +datum=WGS84"
#state <- spTransform(state, CRS(proj4string(equi7t1[["NA"]])))
## Most of North America, including Canada, Mexico and the Caribbean:
state = readOGR("EQUI7_T1_NA_Selection.shp", "EQUI7_T1_NA_Selection")
ov.state <- over(y=state, x=equi7t1[["NA"]])
#ov.state <- gIntersection(state, equi7t1[["NA"]], byid = TRUE)
#str(ov.state@data)

new.dirs = unique(paste0("NA_", equi7t1[["NA"]]$TILE[which(!is.na(ov.state))])) #levels(as.factor(ov$equi7))
#new.dirs = c("NA_061_055","NA_090_043")
## 2071 tiles
x <- lapply(paste0("./", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## Split models otherwise too large in size:
num_splits=30
#mrfX_NASISgg = readRDS.gz("mrfX_NASISgg.rds")
mrfX_NASISgg_final <- split_rf(mrfX_NASISgg, num_splits)
for(j in 1:length(mrfX_NASISgg_final)){
  gm = mrfX_NASISgg_final[[j]]
  saveRDS.gz(gm, file=paste0("mrfX_NASISgg_", j,".rds"))
}
rm(mrfX_NASISgg); rm(mrfX_NASISgg_final)
gc(); gc()
save.image()

#mrfX_NASISpscs = readRDS.gz("mrfX_NASISpscs.rds")
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
#del.lst <- list.files(path="/data/NASIS", pattern=glob2rx("^PSCS_*_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/NASIS", pattern=glob2rx("^TAXgg_*_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

#model.n = "mrfX_NASISgg_"
#varn = "TAXgg"
model.n = "mrfX_NASISpscs_"
varn = "PSCS"
out.path = "/data/NASIS"
for(j in 1:num_splits){
  gm = readRDS.gz(paste0(model.n, j,".rds"))
  cpus = unclass(round((256-50)/(3.5*(object.size(gm)/1e9))))
  sfInit(parallel=TRUE, cpus=ifelse(cpus>46, 46, cpus))
  sfExport("gm", "new.dirs", "split_predict_c", "j", "varn", "out.path")
  sfLibrary(ranger)
  x <- sfClusterApplyLB(new.dirs, fun=function(x){ if(length(list.files(path = paste0(out.path, "/", x, "/"), glob2rx("*.rds$")))<j){ try( split_predict_c(x, gm, in.path="/data/covs1t", out.path=out.path, split_no=j, varn=varn) ) } } )  ## , num.threads=5
  sfRemoveAll()
  sfStop()
  closeAllConnections()
  rm(gm)
  gc()
}

## without parallelization:
# for(j in 1:num_splits){
#   gm = readRDS.gz(paste0(model.n, j,".rds"))
#   x = lapply(new.dirs, function(x){ if(length(list.files(path = paste0(out.path, "/", x, "/"), glob2rx("*.rds$")))<j){ try( split_predict_c(x, gm, in.path="/data/covs1t", out.path=out.path, split_no=j, varn=varn) ) } } )
# }
#sum_predict_ranger(i="NA_060_036", in.path="/data/covs1t", out.path="/data/NASIS", varn="TAXgg", num_splits)
#sum_predict_ranger(i="NA_060_036", in.path="/data/covs1t", out.path="/data/NASIS", varn="PSCS", num_splits)

## Sum up predictions:
sfInit(parallel=TRUE, cpus=25)
sfExport("new.dirs", "sum_predict_ranger", "num_splits", "varn")
sfLibrary(rgdal)
sfLibrary(plyr)
x <- sfLapply(new.dirs, fun=function(x){ try( sum_predict_ranger(x, in.path="/data/covs1t", out.path="/data/NASIS", varn=varn, num_splits) )  } )
sfRemoveAll()
sfStop()

## Create mosaics:
mosaic_tiles_NASIS <- function(j, in.path, varn, latlon=FALSE, tr=0.002083333, r="bilinear", ot="Byte", dstnodata=255, out.path){
  out.tif <- paste0(out.path, varn, "_", j, '_250m_ll.tif')
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    if(latlon==TRUE){
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"'))
    } else {
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"'))
    }
  }
}

levs = list.files(path="./NA_060_036", pattern=glob2rx(paste0("^",varn,"_*_*_*_*.tif$")))
levs = sapply(basename(levs), function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_NASIS", "varn")
out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_NASIS(x, in.path="/data/NASIS/", varn=varn, out.path="/data/NASIS/") )})
sfStop()

save.image()
