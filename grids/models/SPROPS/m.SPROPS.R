## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "ranger", "mda", "psych", "stringr", "caret", "plotKML", "maptools", "maps", "stringr", "R.utils", "grDevices")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/models/SPROPS")
load(".RData")
library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
library(devtools)
#install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
library(xgboost) ## xgboost_0.6-4
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger", ref="forest_memory") ## version to deal with Memory problems
library(ranger) ## ranger_0.6.7
library(caret)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(snowfall)
library(utils)
library(plotKML)
library(R.utils)
library(GSIF)
library(parallel)
library(doParallel)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "gdalwarp"
gdalbuildvrt = "gdalbuildvrt"
system("gdal-config --version")
source("/data/models/wrapper.predict_cs.R")
source("/data/models/saveRDS_functions.R")
source('/data/models/mosaick_functions_ll.R')
source('/data/models/extract_tiled.R')
## metadata:
metasd <- read.csv('/data/GEOG/META_GEOTIFF_1B.csv', stringsAsFactors = FALSE)
sel.metasd = names(metasd)[-sapply(c("FileName","VARIABLE_NAME"), function(x){grep(x, names(metasd))})]
## covariates:
des <- read.csv("/data/models/SoilGrids250m_COVS250m.csv")
mask_value <- as.list(des$MASK_VALUE)
names(mask_value) = des$WORLDGRIDS_CODE

## points:
load("/data/profs/SPROPS/SPROPS.pnts.rda") ## spatial locations only
load("/data/profs/SPROPS/all.pnts.rda")

## spatia overlay (20 mins):
tile.pol = rgdal::readOGR("/data/models/tiles_ll_100km.shp", "tiles_ll_100km")
ov <- extract.tiled(x=SPROPS.pnts, tile.pol=tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=48)
#str(ov)
ovA <- join(all.pnts, ov, type="left", by="LOC_ID")
## 807,962 rows

## Data inspection (final checks)
hist(log1p(ovA$CECSUM))
hist(ovA$BLD)
summary(ovA$BLD) ## many missing values
hist(ovA$CRFVOL) ## probably many missing values are 0?
hist(log1p(ovA$ORCDRC))
hist(ovA$PHIHOX)
summary(ovA$PHIKCL) ## also many missing values
summary(ovA$DEPTH.f)
ovA <- ovA[ovA$DEPTH.f >= 0,]
summary(ovA$DEPTH.f>600)
ovA[which(ovA$DEPTH.f>600)[1],1:20]
## mask out "Artic" data?
#ovA <- ovA[!ovA$SOURCEDB %in% c("Artic"),]
summary(ovA$ORCDRC > 80) ## 5% of profiles >8% ORC

## Fill in all missing CRFVOL values (we assume that if it is missing it is most likely "0"):
ovA$CRFVOL.f = ifelse(is.na(ovA$CRFVOL), 0, ovA$CRFVOL)
summary(ovA$CRFVOL.f)

## relationship between ORC and BLD ----
library(scales)
library(hexbin)
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)
  panel.loess(x, y, ..., col = "black",lty=1,lw=2,span=1/18)
}
pal = R_pal[["bpy_colors"]][1:18]
hexbinplot(ovA$BLD~ovA$ORCDRC, colramp=colorRampPalette(pal), xlab="Organic carbon (permille)", ylab="Bulk density (kg/cubic-m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, asp=1, xbins=30, ybins=30, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))
## Fill in gaps in soil BLD - we can make a PedoTransfer function (which includes soil type predictions)
ov.tax = raster::extract(raster("/data/GEOG/TAXOUSDA_250m_ll.tif"), SPROPS.pnts)
dfs = plyr::join(ovA[,c("SOURCEID", "BLD", "ORCDRC", "CLYPPT", "SNDPPT", "PHIHOX", "DEPTH.f")], data.frame(SOURCEID=SPROPS.pnts$SOURCEID, TAXOUSDA=as.factor(ov.tax)))
dfs = dfs[complete.cases(dfs[,c("ORCDRC", "BLD", "CLYPPT", "SNDPPT", "PHIHOX", "DEPTH.f", "TAXOUSDA")]),]
## 98,650 horizons
ind.tax = data.frame(model.matrix(~TAXOUSDA-1, dfs))
dfs_tbl = cbind(dfs, ind.tax)
saveRDS(dfs_tbl, "wosis_tbl.rds")
saveRDS(ind.tax, "ov_taxousda.rds")
#ind.tax = readRDS("ov_taxousda.rds")
#dfs_tbl = readRDS("wosis_tbl.rds")

fm.BLD = as.formula(paste("BLD ~ ORCDRC + CLYPPT + SNDPPT + PHIHOX + DEPTH.f +", paste(names(ind.tax), collapse="+")))
m.BLD_PTF <- ranger(fm.BLD, dfs_tbl, num.trees = 85, importance='impurity')
#ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
#m.BLD_PTF <- caret::train(BLD~log1p(ORCDRC)+log1p(DEPTH.f)+PHIHOX+SNDPPT+CLYPPT, method="lm", ovA, na.action=na.omit, trControl=ctrl)
m.BLD_PTF
## RMSE = 181 kg/m3
## R-square = 54%
saveRDS.gz(m.BLD_PTF, "m_BLD_PTF.rds")
hist(m.BLD_PTF$predictions)
xl <- as.list(ranger::importance(m.BLD_PTF))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
# ORCDRC     1348034076
# DEPTH.f     503196074
# CLYPPT      326637083
# SNDPPT      289527403
# TAXOUSDA53  254542988
# PHIHOX      252207245
# TAXOUSDA84  122500205
# TAXOUSDA26  104138623
# TAXOUSDA21   89730868
# TAXOUSDA80   68314069

## https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/aggregated/TAXOUSDA.sld
## Udalfs
ind.tax.new = ind.tax[which(ind.tax$TAXOUSDA84==1)[1],]
predict(m.BLD_PTF, cbind(data.frame(ORCDRC=11, CLYPPT=22, PHIHOX=6.5, SNDPPT=35, DEPTH.f=5), ind.tax.new))$predictions
## 1522 kg/m3
predict(m.BLD_PTF, cbind(data.frame(ORCDRC=5, CLYPPT=24, PHIHOX=6.8, SNDPPT=28, DEPTH.f=125), ind.tax.new))$predictions
## 1658 kg/m3

## organic soil:
ind.tax.new = ind.tax[which(ind.tax$TAXOUSDA13==1)[1],]
predict(m.BLD_PTF, cbind(data.frame(ORCDRC=220, CLYPPT=12, PHIHOX=5.5, SNDPPT=35, DEPTH.f=5), ind.tax.new))$predictions
## 803 kg/m3
ind.tax.new = ind.tax[which(ind.tax$TAXOUSDA13==1)[1],]
predict(m.BLD_PTF, cbind(data.frame(ORCDRC=320, CLYPPT=8, PHIHOX=5.5, SNDPPT=45, DEPTH.f=10), ind.tax.new))$predictions
## 778 kg/m3

## Simple correction from Kochy et al (2015):
#round((-0.31 * log1p(220/10) + 1.38)*1000)
m.BLD_ls = loess(BLD ~ ORCDRC, ovA, span=1/18)
predict(m.BLD_ls, data.frame(ORCDRC=220))
## 329.2059
sel.BLD = complete.cases(ovA[,all.vars(fm.BLD)[-1]])
## Take weighted average because the RF-based PTF over-estimates BLD for high ORC
BLD.f = (predict(m.BLD_PTF, ovA[sel.BLD,])$predictions + predict(m.BLD_ls, ovA[sel.BLD,]))/2
BLD.f = ifelse(is.nan(BLD.f), NA, BLD.f)
## replace values where missing:
ovA[sel.BLD,"BLD.f"] = BLD.f 
ovA$BLD.f = round(ifelse(is.na(ovA$BLD), ovA$BLD.f, ovA$BLD))
## Soil organic carbon density:
ovA$OCDENS = round(ovA$ORCDRC/1000 * ovA$BLD.f * (100-ovA$CRFVOL.f)/100)
hist(log1p(ovA$OCDENS))
ovA[1920,c("SOURCEID","BLD","ORCDRC","BLD.f","OCDENS")]
ovA[81051,c("SOURCEID","BLD","ORCDRC","BLD.f","OCDENS")]
View(ovA[which(ovA$OCDENS>200),c("SOURCEID","BLD","ORCDRC","BLD.f","OCDENS")])
## Organic carbon density as function of ORC:
hexbinplot(log1p(ovA$OCDENS)~log1p(ovA$ORCDRC), colramp=colorRampPalette(pal), xlab="log - Organic carbon (permille)", ylab="log - Organic carbon density (kg/cubic-m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, xbins=30, ybins=30, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))
## treshold at ca 12%
expm1(4.8)

write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
unlink("ov.SPROPS_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_SoilGrids250m.csv")
save(ovA, file="ovA.rda")
#load("ovA.rda")
## 1.3GB
#load(".RData")
pnts.ll = ovA[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","CRFVOL","SNDPPT","SLTPPT","CLYPPT","BLD","PHIHOX","PHIKCL","ORCDRC","CECSUM","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","DEPTH","HZDTXT","PHICAL","LOC_ID","UHDICM.f","LHDICM.f","DEPTH.f","CRFVOL.f","BLD.f","OCDENS")]
pnts.ll <- pnts.ll[!is.na(pnts.ll$LONWGS84),]
coordinates(pnts.ll) = ~ LONWGS84+LATWGS84
proj4string(pnts.ll) = CRS("+proj=longlat +datum=WGS84")
writeOGR(pnts.ll, "global_soil_points.gpkg", "global_soil_points", "GPKG")

## ------------- MODEL FITTING -----------

t.vars <- c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD.f", "CECSUM", "OCDENS")
#lapply(ovA[,t.vars], quantile, probs=c(0.01,0.5,0.99), na.rm=TRUE)

z.min <- as.list(c(0,20,20,0,0,0,0,50,0,0))
names(z.min) = t.vars
z.max <- as.list(c(800,110,110,100,100,100,100,3500,2200,10000))
names(z.max) = t.vars
## FIT MODELS:
pr.lst <- basename(list.files(path="/data/stacked250m", ".tif"))
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-unlist(sapply(c("QUAUEA3","LCEE10","N11MSD3","B08CHE3","B09CHE3","CSCMCF5","S01ESA4","S02ESA4","S11ESA4","S12ESA4"), function(x){grep(x, pr.lst)}))]
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ DEPTH.f +', paste(pr.lst, collapse="+")))})
#all.vars(formulaString.lst[[1]])
save.image()

## Median value per cov:
cov.m <- lapply(ovA[,all.vars(formulaString.lst[[1]])[-1]], function(x){quantile(x, probs=c(0.01,0.5,0.99), na.rm=TRUE)})
col.m <- as.data.frame(cov.m)
write.csv(col.m, "covs_quantiles.csv")

## cleanup:
#unlink(list.files(pattern=glob2rx("^mrf.*.rds")))
#unlink(list.files(pattern=glob2rx("^t.mrf.*.rds")))
#unlink(list.files(pattern=glob2rx("^Xgb.*")))
#unlink(list.files(pattern=glob2rx("^mgb.*.rds")))
#unlink(list.files(pattern=glob2rx("^mrf.*.rds")))
#unlink(list.files(pattern=glob2rx("^mrfX_*_*.rds")))
#unlink(list.files(pattern=glob2rx("^RF_fit_*.csv.gz")))
#unlink(list.files(pattern=glob2rx("*_resultsFit.txt")))

## sub-sample to speed up model fitting:
## TAKES >2 hrs to fit all models
Nsub <- 1.2e4 
## Caret training settings (reduce number of combinations to speed up):
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:4, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
rf.tuneGrid <- expand.grid(mtry = seq(5,50,by=5))
## Initiate cluster
cl <- makeCluster(48)
doParallel::registerDoParallel(cl)
## Takes 1 hour to fit all models:
for(j in 1:length(t.vars)){
  out.file = paste0(t.vars[j],"_resultsFit.txt")
  if(!file.exists(out.file)){
    cat("Results of model fitting 'randomForest / XGBoost':\n\n", file=out.file)
    cat("\n", file=out.file, append=TRUE)
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file=out.file, append=TRUE)
    cat("\n", file=out.file, append=TRUE)
    out.rf <- paste0("mrf.",t.vars[j],".rds")
    if(!file.exists(out.rf)|!file.exists(paste0("mgb.",t.vars[j],".rds"))){
      LOC_ID <- ovA$LOC_ID
      dfs <- ovA[,all.vars(formulaString.lst[[j]])]
      sel <- complete.cases(dfs)
      dfs <- dfs[sel,]
      ## optimize mtry parameter:
      if(!file.exists(gsub("mrf","t.mrf",out.rf))){
        t.mrfX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
        saveRDS.gz(t.mrfX, file=gsub("mrf","t.mrf",out.rf))
      } else {
        t.mrfX <- readRDS.gz(gsub("mrf","t.mrf",out.rf))
      }
      ## fit RF model using 'ranger' (fully parallelized)
      ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
      if(!file.exists(paste0("mrf.",t.vars[j],".rds"))){
        mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
        saveRDS.gz(mrfX, file=paste0("mrf.",t.vars[j],".rds"))
      } else {
        mrfX <- readRDS.gz(paste0("mrf.",t.vars[j],".rds"))
      }
      ## Top 15 covariates:
      sink(file=out.file, append=TRUE, type="output")
      print(mrfX)
      cat("\n Variable importance:\n", file=out.file, append=TRUE)
      xl <- as.list(ranger::importance(mrfX))
      print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:35]])))
      ## save fitting success vectors:
      fit.df <- data.frame(LOC_ID=LOC_ID[sel], observed=dfs[,1], predicted=predictions(mrfX))
      unlink(paste0("RF_fit_", t.vars[j], ".csv.gz"))
      write.csv(fit.df, paste0("RF_fit_", t.vars[j], ".csv"))
      gzip(paste0("RF_fit_", t.vars[j], ".csv"))
      if(!file.exists(paste0("mgb.",t.vars[j],".rds"))){
        ## fit XGBoost model (uses all points):
        mgbX <- caret::train(formulaString.lst[[j]], data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid) 
        saveRDS.gz(mgbX, file=paste0("mgb.",t.vars[j],".rds"))
        ## save also binary model for prediction purposes:
        xgb.save(mgbX$finalModel, paste0("Xgb.",t.vars[j]))
      } else {
        mgbX <- readRDS.gz(paste0("mgb.",t.vars[j],".rds"))
      }
      importance_matrix <- xgb.importance(mgbX$coefnames, model = mgbX$finalModel)
      cat("\n", file=out.file, append=TRUE)
      print(mgbX)
      cat("\n XGBoost variable importance:\n", file=out.file, append=TRUE)
      print(importance_matrix[1:25,])
      cat("--------------------------------------\n", file=out.file, append=TRUE)
      sink()
    }    
  }
}
rm(mrfX); rm(mgbX)
stopCluster(cl); closeAllConnections()
save.image()

## ------------- PREDICTIONS -----------

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)

#save(pr.dirs, file="pr.dirs.rda")
## 16,561 dirs
type.lst <- c("Int16","Byte","Byte","Byte","Byte","Byte","Byte","Int16","Int16","Int16")
names(type.lst) = t.vars
mvFlag.lst <- c(-32768, 255, 255, 255, 255, 255, 255, -32768, -32768, -32768)
names(mvFlag.lst) = t.vars

## Run per property (TAKES ABOUT 26 HOURS OF COMPUTING PER PROPERTY)
library(ranger)
library(xgboost)
library(tools)
library(parallel)
library(doParallel)
library(rgdal)
library(plyr)
## Run per property (RF = 20 tiles per minute)
for(j in c("PHIHOX","ORCDRC","OCDENS")){ ## t.vars
#for(j in t.vars){
  try( detach("package:snowfall", unload=TRUE), silent=TRUE)
  try( detach("package:snow", unload=TRUE), silent=TRUE)
  if(j %in% c("PHIHOX","PHIKCL","OCDENS")){ multiplier = 10 }
  if(j %in% c("ORCDRC","CRFVOL","SNDPPT","SLTPPT","CLYPPT","BLD.f","CECSUM")){ multiplier = 1 }
  ## Random forest predictions:
  gm = readRDS.gz(paste0("mrf.", j,".rds"))
  gm1.w = 1/gm$prediction.error
  ## Estimate amount of RAM needed per core
  cpus = unclass(round((500-50)/(3.5*(object.size(gm)/1e9))))
  cl <- parallel::makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
  x = parallel::parLapply(cl, pr.dirs, fun=function(x){ if(any(!file.exists(paste0("/data/tt/SoilGrids250m/predicted250m/", x, "/", j, "_M_sl", 1:7, "_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", split_no=NULL, varn=j, method="ranger", DEPTH.col="DEPTH.f", multiplier=multiplier, rds.file=paste0("/data/tt/SoilGrids250m/predicted250m/", x, "/", x,".rds")) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## XGBoost:
  gm = readRDS.gz(paste0("mgb.", j,".rds"))
  gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
  cpus = unclass(round((500-30)/(3.5*(object.size(gm)/1e9))))
  cl <- parallel::makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
  x = parallel::parLapply(cl, pr.dirs, fun=function(x){ if(any(!file.exists(paste0("/data/tt/SoilGrids250m/predicted250m/", x, "/", j, "_M_sl", 1:7, "_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", split_no=NULL, varn=j, method="xgboost", DEPTH.col="DEPTH.f", multiplier=multiplier, rds.file=paste0("/data/tt/SoilGrids250m/predicted250m/", x, "/", x,".rds")) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## sum up predictions:
  if(is.nan(gm1.w)|is.nan(gm2.w)){ gm1.w = 0.55; gm2.w = 0.45 } 
  ## TH: in some extreme situations it can happen that ranger results in NaN
  library(snowfall)
  sfInit(parallel=TRUE, cpus=48)
  sfExport("pr.dirs", "sum_predict_ensemble", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( sum_predict_ensemble(x, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn=j, num_splits=NULL, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]], rds.file=paste0("/data/tt/SoilGrids250m/predicted250m/", x, "/", x,".rds")) ) } )
  sfStop()
  gc(); gc()
}

## some geotifs got corrupt for unknown reason:
tif.lst <- list.files("/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("OCDENS*.tif"), recursive = TRUE, full.names = TRUE)

corrupt.tif <- function(i, x=NULL, try.raster=TRUE){
  if(try.raster==TRUE){
    try( x <- raster::raster(i) , silent=TRUE)
    if(!class(x)[1]=="RasterLayer"){ 
      return(i)
    }
  } else {
    try( x <- rgdal::GDALinfo(i, silent=TRUE) , silent=TRUE)    
    if(!class(x)[1]=="GDALobj"){ 
      return(i)
    }
  }
}
## takes 30 mins to inspect all tifs
library(snowfall)
sfInit(parallel=TRUE, cpus=48)
sfLibrary(rgdal)
sfLibrary(raster)
sfExport("corrupt.tif", "tif.lst")
del.tif <- unlist(sfLapply(tif.lst, corrupt.tif))
sfStop()
write.csv(del.tif, "corrupt_tif.csv")
unlink(del.tif)

#i = "/data/tt/SoilGrids250m/predicted250m/T46814/OCDENS_M_sl2_T46814.tif"
#GDALinfo(i)
#plot(raster(i))
#del.tif = tif.lst[sapply(c.lst, function(x){grep(x, tif.lst)})]
#unlink(del.tif)

## corrupt or missing tiles:
missing.tiles <- function(varn, pr.dirs){
  dir.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", varn, "*.tif")), full.names=TRUE, recursive=TRUE)
  out.lst <- pr.dirs[which(!pr.dirs %in% basename(dirname(dir.lst)))]
  return(out.lst)
}
sfInit(parallel=TRUE, cpus=length(t.vars))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("missing.tiles", "t.vars", "pr.dirs")
missing.lst <- sfLapply(t.vars, missing.tiles, pr.dirs=pr.dirs)
sfStop()
names(missing.lst) = t.vars
## "NA_101_074" "NA_106_069"

## ORCDRC map can lead to artifacts e.g. 8% ORC for "Shurblands" class in Turmenistan (obvious overprediction)
## this happens mainly because there are not enough training points in drylands hence the model extrapolates;
## to avoid this, we run a post-processing filter that uses precipitation and EVI images to cut out all values that are obvious artifacts (should be 0--5 ORC);
#source("HISTPR_SoilGrids.R")

## clean-up:
#for(i in t.vars){ ## c("ORCDRC", "OCSTHA")){ ##  c("BLD", "ORCDRC", "PHIHOX", "PHIKCL", "SNDPPT", "SLTPPT", "CLYPPT")  
# del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)
#}
# for(i in c("BLD.f", "ORCDRC")){ ## c("BLD", "ORCDRC", "PHIHOX") 
#   del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", i, "_*_*_*_rf?.rds")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
#del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("^CECSUM"), full.names=TRUE, recursive=TRUE)

## delete using date:
orc.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern="ORCDRC", full.names=TRUE, recursive=TRUE)
d.lst <- sapply(orc.lst, function(x){file.mtime(x)})
# #hist(d.lst, breaks=45)
# unlink(orc.lst[d.lst>unclass(as.POSIXct("2017-04-11", format="%Y-%m-%d"))])
# ocd.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern="OCDENS", full.names=TRUE, recursive=TRUE)
# d.lst <- sapply(ocd.lst, function(x){file.mtime(x)})
# unlink(ocd.lst[d.lst>unclass(as.POSIXct("2017-04-11", format="%Y-%m-%d"))])
#ocs.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern="OCSTHA", full.names=TRUE, recursive=TRUE)
#d.lst <- sapply(ocs.lst, function(x){file.mtime(x)})
# unlink(ocs.lst[d.lst>unclass(as.POSIXct("2017-04-11", format="%Y-%m-%d"))])
#del.lst <- unique(dirname((ocs.lst[d.lst>unclass(as.POSIXct("2017-04-17 00:00:00 CEST"))])))
#str(del.lst)
# for(i in c("OCSTHA", "AWCh1", "AWCh2", "AWCh3", "WWP", "AWCtS")){
#   if(i == "OCSTHA"){ 
#     unlink(as.vector(sapply(paste0("sd", 1:6), function(x){ paste0(del.lst, "/", i, "_M_", x, "_", sapply(del.lst, function(k){strsplit(k, "/")[[1]][6]}), ".tif") } )))
#   } else {
#     unlink(as.vector(sapply(1:7, function(x){ paste0(del.lst, "/", i, "_M_sl", x, "_", sapply(del.lst, function(k){strsplit(k, "/")[[1]][6]}), ".tif") } )))
#   }
# }

## Mosaick ----
l.vars = c(as.vector(sapply(c("ORCDRC","PHIHOX","PHIKCL","CRFVOL","SNDPPT","SLTPPT","CLYPPT","BLD.f","CECSUM","OCDENS","TEXMHT","AWCh1", "AWCh2", "AWCh3", "WWP", "AWCtS"), rep, 7)), rep("OCSTHA", 9), "HISTPR", "TAXNWRB", "TAXOUSDA")
#l.vars = c(as.vector(sapply(c("ORCDRC","PHIHOX","PHIKCL","CRFVOL","SNDPPT","SLTPPT","CLYPPT","BLD.f","CECSUM","OCDENS","TEXMHT"), rep, 7)), "HISTPR", "TAXNWRB", "TAXOUSDA")
l.vars.f <- ifelse(l.vars=="BLD.f", "BLDFIE", ifelse(l.vars=="CECSUM", "CECSOL", l.vars))
d.lst = c(rep(paste0("M_sl", 1:7), 16), paste0("M_sd", 1:6), paste0("M_", c(30,100,200), "cm"), "NULL", "NULL", "NULL")
#d.lst = c(rep(paste0("M_sl", 1:7), 11), "NULL", "NULL", "NULL")
filename = paste0(l.vars.f, "_", d.lst, "_250m_ll.tif")
filename = gsub("_NULL", "", filename)
source('/data/models/mosaick_functions_ll.R')
#View(data.frame(l.vars,l.vars.f,d.lst,filename))
library(snowfall)
sfInit(parallel=TRUE, cpus=ifelse(length(l.vars)>28, 28, length(l.vars)))
sfExport("l.vars", "l.vars.f", "d.lst", "make_mosaick_ll", "metasd", "sel.metasd", "filename")
out <- sfClusterApplyLB(1:length(l.vars), function(x){ try( make_mosaick_ll(varn=l.vars[x], i=d.lst[x], in.path="/data/tt/SoilGrids250m/predicted250m", ot=metasd[which(metasd$FileName == filename[x]), "DATA_FORMAT"], dstnodata=metasd[which(metasd$FileName == filename[x]), "NO_DATA"], metadata=metasd[which(metasd$FileName == filename[x]), sel.metasd]) )})
sfStop()
save.image()


## ------------- VISUALIZATIONS -----------

## world plot - overlay and plot points and maps:
xy.pnts <- ovA[!duplicated(ovA$SOURCEID),c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(all.pnts)
#plotKML(xy.pnts["SOURCEDB"], folder.name="Soil properties", file.name="SPROPS_observed.kml")
#zip(zipfile="SPROPS_observed.kmz", files="SPROPS_observed.kml", zip="zip")

require(maptools)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
no.plt <- xy.pnts@coords[,2]>-65 & xy.pnts@coords[,2]<85
png(file = "Fig_global_distribution_SPROPS.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(xy.pnts[xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
points(xy.pnts[!xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
dev.off()

## Cross-validation 10-fold (TH: this does not takes into account high spatial clustering):

source("../../cv/cv_functions.R")
cat("Results of Cross-validation:\n\n", file="resultsCV.txt")
cv_lst <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(!file.exists(paste0("CV_", t.vars[j], ".rda"))){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsCV.txt", append=TRUE)
    cat("\n", file="resultsCV.txt", append=TRUE)
    cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ovA[complete.cases(ovA[,all.vars(formulaString.lst[[j]])]),], nfold=10, idcol="SOURCEID", Log=TRUE)
    sink(file="resultsCV.txt", append=TRUE, type="output")
    print(cv_lst[[j]]$Summary)
    cat("\n", file="resultsCV.txt", append=TRUE)
    sink()
    assign(paste0("CV_", t.vars[j]), cv_lst[[j]])
    save(list=paste0("CV_", t.vars[j]), file=paste0("CV_", t.vars[j], ".rda"))
  }
}

## correlation plots:
source("../plot_hexbin.R")
plt.names <- c("SOC in g/kg", "Soil pH x 10 in H2O", "Soil pH x 10 in KCl", "Coarse fragments in %vol", "Sand fraction in %", "Silt fraction in %", "Clay fraction in %", "Bulk density (FE) in kg / m3", "CEC soil in cmolc/kg") 
names(plt.names) = t.vars
breaks.lst <- list(c(0,5,10,seq(20,1000,length=47)), seq(2.5,9.5,length=50), seq(2.5,9.5,length=50), c(0,1,2,5,seq(8,100,length=46)), seq(0,100,length=50), seq(0,100,length=50), seq(0,100,length=50), seq(450,2200,length=50), c(0,1,2,5,seq(8,450,length=26)))
names(breaks.lst) = t.vars
plt.log <- c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
names(plt.log) = t.vars

for(j in 1:length(t.vars)){
  plot_hexbin(j, breaks.lst[[t.vars[j]]], main=plt.names[t.vars[j]], in.file=paste0("CV_", t.vars[j], ".rda"), log.plot=plt.log[t.vars[j]])
}

## Comparison with indepdent point data not used for model building:
load("/data/profs/SPROPS/tmp/SPROPS.NPDB_V2.rda")
SPROPS.NPDB_V2$LOC_ID = as.factor(paste("ID", SPROPS.NPDB_V2$LONWGS84, SPROPS.NPDB_V2$LATWGS84, sep="_"))
## Correct depth (Canadian data can have negative depths):
z.min.NPDB_V2 <- ddply(SPROPS.NPDB_V2, .(SOURCEID), summarize, aggregated = min(UHDICM, na.rm=TRUE))
z.shift.NPDB_V2 <- join(SPROPS.NPDB_V2[,c("SOURCEID","SOURCEDB")], z.min.NPDB_V2, type="left")$aggregated
z.shift.NPDB_V2 <- ifelse(z.shift.NPDB_V2>0, 0, z.shift.NPDB_V2)
SPROPS.NPDB_V2$UHDICM.f <- SPROPS.NPDB_V2$UHDICM - z.shift.NPDB_V2
SPROPS.NPDB_V2$LHDICM.f <- SPROPS.NPDB_V2$LHDICM - z.shift.NPDB_V2
SPROPS.NPDB_V2$DEPTH.f <- SPROPS.NPDB_V2$UHDICM.f + (SPROPS.NPDB_V2$LHDICM.f - SPROPS.NPDB_V2$UHDICM.f)/2

## overlay points and predict:
xy.NPDB_V2 = SPROPS.NPDB_V2[!duplicated(SPROPS.NPDB_V2$LOC_ID),c("SOURCEID","LOC_ID","SOURCEDB","LONWGS84","LATWGS84")]
coordinates(xy.NPDB_V2) = ~ LONWGS84 + LATWGS84
proj4string(xy.NPDB_V2) = proj4string(tile.pol)
NPDB_V2.ov <- extract.tiled(x=xy.NPDB_V2, tile.pol=tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=54)
NPDB_V2.ovA <- join(SPROPS.NPDB_V2, NPDB_V2.ov, type="left", by="LOC_ID")
gm = readRDS.gz(paste0("mrf.ORCDRC.rds"))
sel.AV = complete.cases(NPDB_V2.ovA[,gm$forest$independent.variable.names])
NPDB_V2.ovA$SG.ORCDRC = NA
NPDB_V2.ovA$SG.ORCDRC[sel.AV] = predict(gm, NPDB_V2.ovA[sel.AV,])$predictions
CV_NPDB_V2_ORCDRC = list(list(Observed=NPDB_V2.ovA$ORCDRC, Predicted=NPDB_V2.ovA$SG.ORCDRC))
save(CV_NPDB_V2_ORCDRC, file="/data/profs/SPROPS/tmp/CV_NPDB_V2_ORCDRC.rda")
## Plot errors:
plot_hexbin("NPDB_V2_ORCDRC", c(0,800), main="NPDB_V2_ORCDRC", in.file="/data/profs/SPROPS/tmp/CV_NPDB_V2_ORCDRC.rda", log.plot=TRUE, out.file = "plot_CV_NPDB_V2_ORCDRC.png")

