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
source("../wrapper.predict_np.R")
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
lapply(ovA[,t.vars], quantile, probs=c(0.01,0.5,0.99), na.rm=TRUE)

z.min <- as.list(c(0,20,20,0,0,0,0,50,0))
names(z.min) = t.vars
z.max <- as.list(c(800,110,110,100,100,100,100,3500,2200))
names(z.max) = t.vars
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + DEPTH.f +', paste(pr.lst, collapse="+")))})
all.vars(formulaString.lst[[1]])

## Median value per cov:
cov.m <- lapply(ovA[,all.vars(formulaString.lst[[1]])[-1]], function(x){quantile(x, probs=c(0.01,0.5,0.99), na.rm=TRUE)})
col.m <- as.data.frame(cov.m)
mask_value <- as.list(des$MASK_VALUE)
names(mask_value) = des$WORLDGRIDS_CODE

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
  rf.tuneGrid <- expand.grid(mtry = seq(4,22,by=2))
  out.rf <- paste0("mrf.",t.vars[j],".rds")
  if(!file.exists(out.rf)){
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    sel <- complete.cases(dfs)
    dfs <- dfs[sel,]
    ## optimize mtry parameter:
    t.mrfX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid) 
    ## fit RF model using 'ranger' (fully parallelized)
    ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
    saveRDS(t.mrfX, file=gsub("mrf","t.mrf",out.rf))
    mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=291) ## 
    saveRDS(mrfX, file=paste0("mrf.",t.vars[j],".rds"))
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
    ## fit XGBoost model (uses all points):
    mgbX <- caret::train(formulaString.lst[[j]], data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid) 
    saveRDS(mgbX, file=paste0("mgb.",t.vars[j],".rds"))
    ## save also binary model for prediction purposes:
    xgb.save(mgbX$finalModel, paste0("Xgb.",t.vars[j]))
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

mrfX_lst <- list.files(pattern="mrf.")
mgbX_lst <- list.files(pattern="mgb.")
names(mrfX_lst) <- paste(sapply(mrfX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))
names(mgbX_lst) <- paste(sapply(mgbX_lst, function(x){strsplit(x, "\\.")[[1]][2]}))

## ------------- PREDICTIONS -----------

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs1t", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
#save(pr.dirs, file="pr.dirs.rda")
## 16,561 dirs

## Test it:
system.time( wrapper.predict_np(i="NA_060_036", varn="ORCDRC", gm1=mrfX_lst[[j]], gm2=mgbX_lst[[j]], in.path="/data/covs1t", out.path="/data/predicted", zmin=50, zmax=3500, mask_value=mask_value) )

## Run per property:
library(snowfall)
for(j in t.vars){
  gc(); gc()
  #gm1 <- mrfX_lst[[j]]
  gm1 <- readRDS(mrfX_lst[[j]])
  gm2 <- readRDS(mgbX_lst[[j]])
  ## divide RAM by size of model:
  cpus = unclass(round(256/(3*(object.size(gm1)/1e9+object.size(gm2)/1e9))))
  cpus <- ifelse(cpus>48, 48, cpus)
  #cl <- makeCluster(s48-cpus)
  #registerDoParallel(cl)
  sfInit(parallel=TRUE, cpus=cpus) # cpus=48, type="MPI"
  ## Error in mpi.send(x = serialize(obj, NULL), type = 4, dest = dest, tag = tag, : long vectors not supported yet: memory.c:3361
  sfExport("wrapper.predict_np", "pr.dirs", "mask_value", "z.min", "z.max", "gm1", "gm2", "j")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  sfLibrary(ranger)
  sfLibrary(xgboost)
  sfLibrary(caret)
  ## exach process up to 5GB (gm1 largest object)
  ## Can take up to an hour to export all objects to all CPUs
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( wrapper.predict_np(i=x, varn=j, gm1=gm1, gm2=gm2, in.path="/data/covs1t", out.path="/data/predicted", zmin=z.min[[j]], zmax=z.max[[j]], mask_value=mask_value) )  } )
  sfStop()
  #stopCluster(cl)
}
stopCluster(cl)

## clean-up:
# for(i in t.vars){ ## c("BLD", "ORCDRC", "PHIHOX") 
#  del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

# for(i in t.vars){
#  del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

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

