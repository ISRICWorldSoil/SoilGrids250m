## Fit models for depth to bedrock and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org & shanggv@hotmail.com

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

plotKML.env(convert="convert", show.env=FALSE)
options(bitmapType='cairo')
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_cs.R")
source("../saveRDS_functions.R")
load("../equi7t3.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")
mask_value <- as.list(des$MASK_VALUE)
names(mask_value) = des$WORLDGRIDS_CODE

## points:
load("../../profs/BDR/BDR.pnts.rda")
load("../../profs/BDR/BDR_all.pnts.rda")
ov <- extract.equi7t3(x=BDR.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=48) 
## TAKES CA 30 mins!
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
## TH: filter out few very high values?
range(ovA$BDTICM, na.rm=TRUE)
## Remove any negative values:
summary(ovA$BDTICM)
ovA$BDTICM <- ifelse(ovA$BDTICM<0, NA, ovA$BDTICM)
## use log to put more emphasis on lower values?
#ovA$logBDTICM <- log1p(ovA$BDTICM)
#hist(log1p(ovA$BDTICM), breaks=40, col="grey", xlab="log-BDTICM", main="Histogram")

## ------------- MODEL FITTING -----------

t.vars <- c("BDRICM", "BDRLOG", "BDTICM")
z.min <- as.list(c(0,0,0))
names(z.min) = t.vars
z.max <- as.list(c(200,100,Inf))
names(z.max) = t.vars

pr.lst <- as.character(des$WORLDGRIDS_CODE)
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-unlist(sapply(c("QUAUEA3","VOLNOA3","C??GLC5"), function(x){grep(glob2rx(x), pr.lst)}))]
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ ', paste(pr.lst, collapse="+")))}) ## LATWGS84 +
## For BDTICM we can not use latitude, quakes and volcano density as predictor because points are somewhat clustered in LAT space --> predictions lead to artifacts 
#formulaString.lst[[3]] = as.formula(paste(t.vars[3], ' ~ ', paste(pr.lst, collapse="+")))

## sub-sample to speed up model fitting:
Nsub <- 1.5e4 
## Initiate cluster
cl <- makeCluster(48)
registerDoParallel(cl)
## Write results of model fitting into a text file:
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="BDR_resultsFit.txt")
for(j in 1:length(t.vars)){
  cat("\n", file="BDR_resultsFit.txt", append=TRUE)
  cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="BDR_resultsFit.txt", append=TRUE)
  cat("\n", file="BDR_resultsFit.txt", append=TRUE)
  LOC_ID <- ovA$LOC_ID
  ## Caret training settings (reduce number of combinations to speed up):
  ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
  gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
  rf.tuneGrid <- expand.grid(mtry = seq(10,60,by=5))
  out.rf <- paste0("mrf.",t.vars[j],".rds")
  if(!file.exists(out.rf)){
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    sel <- complete.cases(dfs)
    dfs <- dfs[sel,]
    ## optimize mtry parameter:
    if(!file.exists(gsub("mrf","t.mrf",out.rf))){
      t.mrfX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid)
      saveRDS.gz(t.mrfX, file=gsub("mrf","t.mrf",out.rf))
    } else {
      t.mrfX <- readRDS.gz(gsub("mrf","t.mrf",out.rf))
    }
    ## fit RF model using 'ranger' (fully parallelized)
    mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=300)
    saveRDS.gz(mrfX, file=paste0("mrf.",t.vars[j],".rds"))
    ## Top 15 covariates:
    sink(file="BDR_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="BDR_resultsFit.txt", append=TRUE)
    xl <- as.list(ranger::importance(mrfX))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
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
    cat("\n", file="BDR_resultsFit.txt", append=TRUE)
    print(mgbX)
    cat("\n XGBoost variable importance:\n", file="BDR_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:15,])
    cat("--------------------------------------\n", file="BDR_resultsFit.txt", append=TRUE)
    sink()
  }
}
rm(mrfX); rm(mgbX)
stopCluster(cl); closeAllConnections()

mrfX_lst <- list.files(pattern="^mrf.")
mgbX_lst <- list.files(pattern="^mgb.")
names(mrfX_lst) <- paste(sapply(mrfX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))
names(mgbX_lst) <- paste(sapply(mgbX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))

## ------------- PREDICTIONS -----------

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs1t", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
#save(pr.dirs, file="pr.dirs.rda")
## 16,561 dirs

## Split RF models (otherwise memory limit problems):
num_splits = 3
for(j in t.vars){
  mrfX = readRDS.gz(mrfX_lst[[j]])
  mrfX_final <- split_rf(mrfX, num_splits)
  for(k in 1:length(mrfX_final)){
    gm = mrfX_final[[k]]
    saveRDS.gz(gm, file=paste0("mrfX_", k, "_", j,".rds"))
  }
}
rm(mrfX); rm(mrfX_final)
gc(); gc()
save.image()

type.lst <- c("Byte", "Byte", "Int32")
names(type.lst) = t.vars
mvFlag.lst <- c(255, 255, -32768)
names(mvFlag.lst) = t.vars

## Run per property (TAKES ABOUT 20-30 HOURS OF COMPUTING PER PROPERTY)
## To avoid memory problems with RF objects, we split predictions into e.g. 3 parts
for(j in t.vars[3]){ # t.vars
  if(j=="BDRLOG"){ multiplier = 100 }
  if(j %in% c("BDRICM","BDTICM")){ multiplier = 1 }
  ## Random forest predictions (splits):
  gm = readRDS.gz(mrfX_lst[[j]])
  gm1.w = 1/gm$prediction.error
  rm(gm)
  for(k in 1:num_splits){
    gm = readRDS.gz(paste0("mrfX_", k, "_", j,".rds"))
    gc()
    cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
    sfInit(parallel=TRUE, cpus=ifelse(cpus>46, 46, cpus))
    sfExport("gm", "pr.dirs", "split_predict_n", "j", "k", "multiplier")
    sfLibrary(ranger)
    x <- sfLapply(pr.dirs, fun=function(x){ if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ try( split_predict_n(x, gm, in.path="/data/covs1t", out.path="/data/predicted", split_no=k, varn=j, method="ranger", multiplier=multiplier, depths=FALSE) ) } } )
    sfStop()
    rm(gm)
  }
  ## XGBoost:
  gm = readRDS.gz(paste0("mgb.", j,".rds"))
  gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
  cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
  gc()
  sfInit(parallel=TRUE, cpus=ifelse(cpus>46, 46, cpus))
  sfExport("gm", "pr.dirs", "split_predict_n", "j", "multiplier")
  sfLibrary(xgboost)
  x <- sfLapply(pr.dirs, fun=function(x){ try( if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ split_predict_n(x, gm, in.path="/data/covs1t", out.path="/data/predicted", varn=j, method="xgboost", multiplier=multiplier, depths=FALSE) } ) } )
  sfStop()
  rm(gm)
  ## sum up predictions:
  sfInit(parallel=TRUE, cpus=45)
  sfExport("pr.dirs", "sum_predict_ensemble", "num_splits", "j", "z.min", "z.max", "gm1.w", "gm2.w", "mvFlag.lst", "type.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ sum_predict_ensemble(x, in.path="/data/covs1t", out.path="/data/predicted", varn=j, num_splits, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]], depths=FALSE) } )  } )
  sfStop()
}

## ------------- VISUALIZATIONS -----------

x <- readGDAL("/data/predicted/NA_060_036/BDRLOG_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDRLOG_M_NA_060_036.kml", folder.name="R horizon", colour=band1, z.lim=c(0,100), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="BDRLOG_M_NA_060_036.png")
x <- readGDAL("/data/predicted/NA_060_036/BDTICM_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDTICM_M_NA_060_036.kml", folder.name="Absolute depth in cm", colour=band1, colour_scale=SAGA_pal[[1]], raster_name="BDTICM_M_NA_060_036.png", z.lim=c(0,7000))

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
    cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ovA, nfold=10, idcol="SOURCEID", h2o=TRUE, Log=TRUE)
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
breaks.lst <- list(c(seq(0,250,length=50)), seq(0,1,length=50), seq(0,50000))
names(breaks.lst) = t.vars
plt.log <- c(FALSE, FALSE, TRUE)
names(plt.log) = t.vars

for(j in 1:length(t.vars)){
  plot_hexbin(j, breaks.lst[[t.vars[j]]], main=plt.names[t.vars[j]], in.file=paste0("CV_", t.vars[j], ".rda"), log.plot=plt.log[t.vars[j]])
}

## clean-up:
# for(i in c("BDRICM", "BDRLOG", "BDTICM")){ 
#   del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
# for(i in c("BDRICM", "BDRLOG", "BDTICM")){ 
#   del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
