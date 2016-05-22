## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

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
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source('/data/models/objects_in_memory.R')
#source("../wrapper.predict_np.R")
source("../wrapper.predict_cs.R")
source("../saveRDS_functions.R")

load("../equi7t3.rda")
load("../equi7t1.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## points:
load("../../profs/SPROPS/SPROPS.pnts.rda")
load("../../profs/SPROPS/all.pnts.rda")
ov <- extract.equi7(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7=equi7t3, path="/data/covs", cpus=48) 
#str(ov)
ovA <- join(all.pnts, ov, type="left", by="LOC_ID")
## 752,161 obs
for(i in des$WORLDGRIDS_CODE){ ovA[,i] <- ifelse(ovA[,i]<= -10000, NA, ovA[,i])  }
## Check values:
hist(log1p(ovA$CECSUM))
hist(ovA$BLD)
summary(ovA$BLD)
hist(log1p(ovA$ORCDRC))
hist(ovA$PHIHOX)
summary(ovA$PHIKCL)
summary(ovA$DEPTH.f)
## mask out datasets that still need to be checked:
ovA <- ovA[!ovA$SOURCEDB %in% c("Artic"),]
## 765,539

write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
unlink("ov.SPROPS_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
#load("ovA.rda")
## 1.3GB

## ------------- MODEL FITTING -----------

t.vars <- c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM")
#lapply(ovA[,t.vars], quantile, probs=c(0.01,0.5,0.99), na.rm=TRUE)

z.min <- as.list(c(0,20,20,0,0,0,0,50,0))
names(z.min) = t.vars
z.max <- as.list(c(800,110,110,100,100,100,100,3500,2200))
names(z.max) = t.vars
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ DEPTH.f +', paste(pr.lst, collapse="+")))}) ## LATWGS84 +
#all.vars(formulaString.lst[[1]])

## Median value per cov:
#cov.m <- lapply(ovA[,all.vars(formulaString.lst[[1]])[-1]], function(x){quantile(x, probs=c(0.01,0.5,0.99), na.rm=TRUE)})
#col.m <- as.data.frame(cov.m)

## sub-sample to speed up model fitting:
Nsub <- 1.5e4 
## Initiate cluster
cl <- makeCluster(48)
registerDoParallel(cl)
## Takes 1 hour to fit all models:
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="SPROPS_resultsFit.txt")
for(j in 1:length(t.vars)){
  cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
  cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="SPROPS_resultsFit.txt", append=TRUE)
  cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
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
    ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
    mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=300)  
    saveRDS.gz(mrfX, file=paste0("mrf.",t.vars[j],".rds"))
    ## Top 15 covariates:
    sink(file="SPROPS_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
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
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(mgbX)
    cat("\n XGBoost variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:15,])
    cat("--------------------------------------\n", file="SPROPS_resultsFit.txt", append=TRUE)
    sink()
  }
}
rm(mrfX); rm(mgbX)
stopCluster(cl); closeAllConnections()

mrfX_lst <- list.files(pattern="^mrf.")
mgbX_lst <- list.files(pattern="^mgb.")
names(mrfX_lst) <- paste(sapply(mrfX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))
names(mgbX_lst) <- paste(sapply(mgbX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))
save.image()

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

type.lst <- c("Int16","Byte","Byte","Byte","Byte","Byte","Byte","Int16","Int16")
names(type.lst) = t.vars
mvFlag.lst <- c(-32768, 255, 255, 255, 255, 255, 255, -32768, -32768)
names(mvFlag.lst) = t.vars

## Run per property (TAKES ABOUT 20-30 HOURS OF COMPUTING PER PROPERTY)
## To avoid memory problems with RF objects, we split into e.g. 3 parts
for(j in t.vars[c(5,6,7,8)]){
  if(j=="PHIHOX"|j=="PHIKCL"){ multiplier = 10 }
  if(j %in% c("P", "S", "B", "Cu", "Zn")){ multiplier = 100 }
  if(j %in% c("ORCDRC","CRFVOL","SNDPPT","SLTPPT","CLYPPT","BLD","CECSUM")){ multiplier = 1 }
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
    x <- sfClusterApplyLB(pr.dirs, fun=function(x){ if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ try( split_predict_n(x, gm, in.path="/data/covs1t", out.path="/data/predicted", split_no=k, varn=j, method="ranger", multiplier=multiplier) ) } } )
    sfStop()
    closeAllConnections()
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
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ split_predict_n(x, gm, in.path="/data/covs1t", out.path="/data/predicted", varn=j, method="xgboost", multiplier=multiplier) } ) } )
  sfStop()
  rm(gm)
  closeAllConnections()
  ## sum up predictions:
  sfInit(parallel=TRUE, cpus=45)
  sfExport("pr.dirs", "sum_predict_ensemble", "num_splits", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( if(length(list.files(path = paste0("/data/predicted/", x, "/"), glob2rx(paste0("^",j,"*.tif$"))))==0){ sum_predict_ensemble(x, in.path="/data/covs1t", out.path="/data/predicted", varn=j, num_splits, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]]) } )  } )
  sfStop()
}

## corrupt or missing tiles:
missing.tiles <- function(varn, pr.dirs){
  dir.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", varn, "*.tif")), full.names=TRUE, recursive=TRUE)
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

## clean-up:
# for(i in c("BLD", "ORCDRC", "PHIHOX", "PHIKCL", "SNDPPT", "SLTPPT", "CLYPPT")){ ## c("CECSUM","CRFVOL")  
#   del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
# for(i in c("CECSUM", "CRFVOL")){ ## c("BLD", "ORCDRC", "PHIHOX") 
#   del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "_*_*_*_rf?.rds")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }

# for(i in t.vars){
#  del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

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

## Cross-validation 10-fold (TH: this does not account for high spatial clustering):

