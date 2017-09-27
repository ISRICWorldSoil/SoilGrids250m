## Fit models for depth to bedrock and generate predictions - SoilGrids250m
## Wei.Shangguan (shgwei@mail.sysu.edu.cn) & Tom.Hengl@isric.org 

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "hexbin", "gridExtra", "ranger", "xgboost", "stringr", "caret", "plotKML", "maptools", "maps")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/models/BDR")
load(".RData")
library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
library(devtools)
#devtools::install_github('dmlc/xgboost')
library(xgboost)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger", ref="forest_memory") ## version to deal with Memory problems
library(ranger)
library(caret)
library(hexbin)
library(gridExtra)
library(snowfall)
library(utils)
library(plotKML)
library(GSIF)
library(R.utils)
#library(doParallel)

plotKML.env(convert="convert", show.env=FALSE)
options(bitmapType='cairo')
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

## load points:
load("/data/profs/BDR/BDR.pnts.rda")
load("/data/profs/BDR/BDR_all.pnts.rda")
str(BDR.pnts@data)

## OVERLAY (takes ca 1 hr):
tile.pol = rgdal::readOGR("/data/models/tiles_ll_100km.shp", "tiles_ll_100km")
#tile.pol = readRDS("stacked250m_tiles_pol.rds")
ov <- extract.tiled(x=BDR.pnts, tile.pol=tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=48)

#sfStop()
ov$LOC_ID = BDR.pnts$LOC_ID
#str(ov)
ovA <- join(BDR_all.pnts, ov, type="left", by="LOC_ID")
write.csv(ovA, file="ov.BDR_SoilGrids250m.csv")
unlink("ov.BDR_SoilGrids250m.csv.gz")
gzip("ov.BDR_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
#load("ovA.rda")
## 1.3GB object

## BDRICM = depth to bedrock until 250 cm (censored data)
## BDRLOG = occurrence of R horizon 0 / 1
## BDTICM = absolute depth to bedrock
hist(ovA$BDRICM)
## filter out few very high values?
range(ovA$BDTICM, na.rm=TRUE)
## Remove any negative values:
summary(ovA$BDTICM)
ovA$BDTICM <- ifelse(ovA$BDTICM<0, NA, ovA$BDTICM)
#ovA$logBDTICM <- log1p(ovA$BDTICM)
hist(log1p(ovA$BDTICM), breaks=40, col="grey", xlab="log-BDTICM", main="Histogram")
## mask out water bodies and permanent ice:
#ovA$OCCGSW7.tif = ifelse(ovA$OCCGSW7.tif>100, 0, ovA$OCCGSW7.tif)
#ovA$LCEE10.tif = ifelse(ovA$LCEE10.tif==220, NA, ovA$LCEE10.tif)
#summary(is.na(ovA$LCEE10.tif))
## 44,489 fall in permanent ice or water
save.image()

## ------------- MODEL FITTING -----------
## TAKES ca 1 hr

t.vars <- c("BDRICM", "BDRLOG", "BDTICM")
z.min <- as.list(c(0,0,0))
names(z.min) = t.vars
z.max <- as.list(c(200,100,Inf))
names(z.max) = t.vars

pr.lst <- basename(list.files(path="/data/stacked250m", ".tif"))
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-unlist(sapply(c("QUAUEA3","LCEE10","N11MSD3","CSCMCF5","B02CHE3","B08CHE3","B09CHE3","S01ESA4","S02ESA4","S11ESA4","S12ESA4","BICUSG"), function(x){grep(x, pr.lst)}))]
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ ', paste(pr.lst, collapse="+")))}) ## LATWGS84 +
## For BDTICM we can not use latitude, quakes and volcano density as predictor because points are somewhat clustered in LAT space --> predictions lead to artifacts 

## sub-sample to speed up model fitting:
Nsub <- 1.2e4 
## Caret training settings (reduce number of combinations to speed up):
ctrl <- trainControl(method="repeatedcv", number=4, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:4, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
rf.tuneGrid <- expand.grid(mtry = seq(5,50,by=5))
max.points = 8e5

## clean-up:
#unlink(list.files(pattern="^mrf."))
#unlink(list.files(pattern="^t.mrf."))
#unlink(list.files(pattern="^mgb."))

## Initiate cluster
cl <- parallel::makeCluster(48)
doParallel::registerDoParallel(cl)
## Write results of model fitting into a text file:
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="BDR_resultsFit.txt")
for(j in 1:length(t.vars)){
  cat("\n", file="BDR_resultsFit.txt", append=TRUE)
  cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="BDR_resultsFit.txt", append=TRUE)
  cat("\n", file="BDR_resultsFit.txt", append=TRUE)
  LOC_ID <- ovA$LOC_ID
  out.rf <- paste0("mrf.",t.vars[j],".rds")
  if(!file.exists(out.rf) | !file.exists(paste0("mgb.",t.vars[j],".rds"))){
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
    if(!file.exists(out.rf)){
      ## fit RF model using 'ranger' (fully parallelized)
      mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)
      saveRDS.gz(mrfX, file=paste0("mrf.",t.vars[j],".rds"))
      #mrfX <- readRDS.gz(file=paste0("mrf.",t.vars[j],".rds"))
      ## Top 25 covariates:
    } else {
      mrfX = readRDS.gz(paste0("mrf.",t.vars[j],".rds"))
    }
    sink(file="BDR_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="BDR_resultsFit.txt", append=TRUE)
    xl <- as.list(ranger::importance(mrfX))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
    ## save fitting success vectors:
    if(!file.exists(paste0("RF_fit_", t.vars[j], ".csv.gz"))){
      fit.df <- data.frame(LOC_ID=LOC_ID[sel], observed=dfs[,1], predicted=predictions(mrfX))
      unlink(paste0("RF_fit_", t.vars[j], ".csv.gz"))
      write.csv(fit.df, paste0("RF_fit_", t.vars[j], ".csv"))
      gzip(paste0("RF_fit_", t.vars[j], ".csv"))
    }
    gc(); gc()
    if(!file.exists(paste0("mgb.",t.vars[j],".rds"))){
      ## fit XGBoost model (uses all points):
      #mgbX <- caret::train(formulaString.lst[[j]], data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid)
      mgbX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), size=max.points),], method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid)
      saveRDS.gz(mgbX, file=paste0("mgb.",t.vars[j],".rds"))
      ## save also binary model for prediction purposes:
      xgb.save(mgbX$finalModel, paste0("Xgb.",t.vars[j]))
    } else {
      mgbX <- readRDS.gz(paste0("mgb.",t.vars[j],".rds"))
    }
    importance_matrix <- xgb.importance(mgbX$coefnames, model = mgbX$finalModel)
    cat("\n", file="BDR_resultsFit.txt", append=TRUE)
    #mgbX <- readRDS.gz(paste0("mgb.",t.vars[j],".rds"))
    print(mgbX)
    cat("\n XGBoost variable importance:\n", file="BDR_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:25,])
    cat("--------------------------------------\n", file="BDR_resultsFit.txt", append=TRUE)
    sink()
  }
  gc(); gc()
}
rm(mrfX); rm(mgbX)
stopCluster(cl); closeAllConnections()
save.image()

## ------------- PREDICTIONS -----------

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("^T*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
#save(pr.dirs, file="pr.dirs.rda")
## 18,653 dirs

## Predictions:
type.lst <- c("Byte", "Byte", "Int32")
names(type.lst) = t.vars
mvFlag.lst <- c(255, 255, -32768)
names(mvFlag.lst) = t.vars
source("predict.BDR.R")
save.image()

## rename files:
#lst.tif = list.files("/data/tt/SoilGrids250m/predicted250m/", "BDRICM_M_", full.names=TRUE, recursive = TRUE)
#file.rename(lst.tif, gsub("_M_sl1_", "_M_", lst.tif))

## Mosaick:
r <- raster("/data/stacked250m/LCEE10.tif")
cellsize = res(r)[1]
te = paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")
filename = paste0(t.vars, "_M_250m_ll.tif")

sfInit(parallel=TRUE, cpus=ifelse(length(t.vars)>45, 45, length(t.vars)))
sfExport("t.vars", "mvFlag.lst", "make_mosaick_ll", "type.lst", "metasd", "sel.metasd", "filename", "te", "cellsize")
out <- sfClusterApplyLB(1:length(t.vars), function(x){ try( make_mosaick_ll(varn=t.vars[x], i="M", in.path="/data/tt/SoilGrids250m/predicted250m", tr=cellsize, te=te, ot=metasd[which(metasd$FileName == filename[x]), "DATA_FORMAT"], dstnodata=metasd[which(metasd$FileName == filename[x]), "NO_DATA"], metadata=metasd[which(metasd$FileName == filename[x]), sel.metasd]) )})
sfStop()

## ------------- VISUALIZATIONS -----------

x <- readGDAL("/data/tt/SoilGrids250m/predicted250m/T34122/BDRLOG_M_T34122.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDRLOG_M_T34122.kml", folder.name="R horizon", colour=band1, z.lim=c(0,100), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="BDRLOG_M_T34122.png")
x <- readGDAL("/data/tt/SoilGrids250m/predicted250m/T34122/BDTICM_M_T34122.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDTICM_M_T34122.kml", folder.name="Absolute depth in cm", colour=band1, colour_scale=SAGA_pal[[1]], raster_name="BDTICM_M_T34122.png", z.lim=c(0,7000))

## world plot - overlay and plot points and maps:
xy.pnts <- join(ovA[complete.cases(ovA[,c("SOURCEID","SOURCEDB")]),], as.data.frame(BDR.pnts[c("SOURCEID")]), type="left", match="first")
xy.pnts <- xy.pnts[!is.na(xy.pnts$LONWGS84),]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(BDR.pnts)
require(maptools)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
no.plt <- xy.pnts@coords[,2]>-65 & xy.pnts@coords[,2]<85
png(file = "Fig_global_distribution_BDR.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
## profile data:
points(xy.pnts[-which(xy.pnts$SOURCEDB=="Wells"|xy.pnts$SOURCEDB=="Simulated"&!no.plt),], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
## Wells data
points(xy.pnts[-which(!xy.pnts$SOURCEDB=="Wells"|xy.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("blue", 0.6), cex=.8, col="black")
points(xy.pnts[which(xy.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()

## world plot - overlay and plot points and maps:
shp.pnts <- ovA[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84",t.vars)]
shp.pnts <- xy.pnts[complete.cases(xy.pnts[,c("LONWGS84","LATWGS84")]),]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(BDR.pnts)
unlink("Soil_depth_training_points.shp")
writeOGR(xy.pnts, "Soil_depth_training_points.shp", "Soil_depth_training_points", "ESRI Shapefile")
unlink("Soil_depth_training_points.zip")
system("7za a Soil_depth_training_points.zip Soil_depth_training_points.*")

## ------------- CROSS-VALIDATION -----------
source("../../cv/cv_functions.R")

cat("Results of Cross-validation:\n\n", file="resultsCV_BDR.txt")
cv_lst <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(!file.exists(paste0("CV_", t.vars[j], ".rda"))){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsCV_BDR.txt", append=TRUE)
    cat("\n", file="resultsCV_BDR.txt", append=TRUE)
    ## exclude simulated points from Cross-validation:
    cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ovA[complete.cases(ovA[,all.vars(formulaString.lst[[j]])]),], nfold=3, idcol="SOURCEID", Log=TRUE)
    sink(file="resultsCV_BDR.txt", append=TRUE, type="output")
    print(cv_lst[[j]]$Summary)
    cat("\n", file="resultsCV_BDR.txt", append=TRUE)
    sink()
    assign(paste0("CV_", t.vars[j]), cv_lst[[j]])
    save(list=paste0("CV_", t.vars[j]), file=paste0("CV_", t.vars[j], ".rda"))
  }
}

source("../plot_hexbin.R")
plt.names <- c("Depth to bedrock (up to 250 cm)", "Occurrence of the R horizon", "Absolute depth to bedrock (in cm)")
names(plt.names) = t.vars
breaks.lst <- list(c(seq(0,280,length=20)), seq(0,1,length=30), seq(0,50000,length=50))
names(breaks.lst) = t.vars
plt.log <- c(FALSE, FALSE, TRUE)
names(plt.log) = t.vars
colorcut.lst = list(c(0,0.005,0.022,0.045,0.086,0.11,0.5,0.75,1), c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1), c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))
names(colorcut.lst) = t.vars

for(j in 1:length(t.vars)){
  plot_hexbin(t.vars[j], breaks.lst[[t.vars[j]]], main=plt.names[t.vars[j]], in.file=paste0("CV_", t.vars[j], ".rda"), log.plot=plt.log[t.vars[j]], colorcut=colorcut.lst[[t.vars[j]]])
}

## without the map from Pelletier et al.
formulaString.Pell = as.formula(paste('BDTICM ~ ', paste(pr.lst[-which(pr.lst %in% "ASSDAC3")], collapse="+")))
CV_Pell <- cv_numeric(formulaString.Pell, rmatrix=ovA[complete.cases(ovA[,all.vars(formulaString.Pell)]),], nfold=3, idcol="SOURCEID", Log=TRUE)
CV_Pell$Summary
save(CV_Pell, file="CV_Pell.rda")
plot_hexbin("Pell", breaks.lst[["BDTICM"]], main="Without Pelletier et al.", in.file="CV_Pell.rda", log.plot=plt.log[["BDTICM"]])
cat(paste("Variable:", all.vars(formulaString.Pell)[1]), file="resultsCV_BDR.txt", append=TRUE)
cat("\n", file="resultsCV_BDR.txt", append=TRUE)
sink(file="resultsCV_BDR.txt", append=TRUE, type="output")
print(CV_Pell$Summary)
cat("\n", file="resultsCV_BDR.txt", append=TRUE)
sink()

# ## clean-up:
# for(i in c("BDRICM", "BDRLOG", "BDTICM")){
#   del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
