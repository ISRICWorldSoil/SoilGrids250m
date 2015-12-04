## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
#library(randomForestSRC)
library(snowfall)
library(utils)
library(plotKML)
library(GSIF)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
#source("../wrapper.predict_nCSV.R")
source("../wrapper.predict_n.R")
load("../equi7t3.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## points:
load("../../profs/SPROPS/SPROPS.pnts.rda")
load("../../profs/SPROPS/all.pnts.rda")
ov <- extract.equi7t3(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
#str(ov)
ovA <- join(all.pnts, ov, type="left", by="LOC_ID")
for(i in des$WORLDGRIDS_CODE){ ovA[,i] <- ifelse(ovA[,i]<= -10000, NA, ovA[,i])  }
write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
unlink("ov.SPROPS_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
#load("ovA.rda")
## 1.3GB

t.vars <- c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM")
summary(ovA$CECSUM)
summary(ovA$BLD)
summary(ovA$ORCDRC)
summary(ovA$PHIHOX)
summary(ovA$PHIKCL)

z.min <-c(0,20,20,0,0,0,0,0,0)
z.max <-c(800,110,110,100,100,100,100,3500,2200)
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + DEPTH +', paste(pr.lst, collapse="+")))})
#all.vars(formulaString.lst[[1]])

## H2O package more suited for large data (http://www.r-bloggers.com/benchmarking-random-forest-implementations/)
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)

## We fit two alternative models - RF and Deeplearning (http://www.rdocumentation.org/packages/h2o/functions/h2o.deeplearning)
cat("Results of 'h2o.randomForest / Deeplearning':\n\n", file="resultsFit.txt")
mrfX_path <- rep(list(NULL), length(t.vars))
mdLX_path <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(is.null(mrfX_path[[j]])){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsFit.txt", append=TRUE)
    cat("\n", file="resultsFit.txt", append=TRUE)
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    dfs.hex <- as.h2o(dfs[complete.cases(dfs),], conn = h2o.getConnection(), destination_frame = "dfs.hex")
    #str(dfs.hex@mutable$col_names)
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex, importance=TRUE) ## TAKES ONLY 2-3 MINS even with 1M points - very fast!!
    print(mrfX)
    sink(file="resultsFit.txt", append=TRUE, type="output")
    cat("\n", file="resultsFit.txt", append=TRUE)
    print(mrfX@model$variable_importances)
    cat("\n", file="resultsFit.txt", append=TRUE)
    mdLX <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex)
    print(mdLX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    sink()
    mrfX_path[[j]] = h2o.saveModel(mrfX, path="./", force=TRUE)
    mdLX_path[[j]] = h2o.saveModel(mdLX, path="./", force=TRUE)
  }
}
names(mrfX_path) = t.vars
names(mdLX_path) = t.vars

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
## 1km resolution (10-15 times faster):
system.time( wrapper.predict_n(i="NA_075_066", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )
system.time( wrapper.predict_n(i="NA_063_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )
## Run all tiles one by one:
x <- lapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )})

## 250m resolution:
system.time( wrapper.predict_n(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )
x <- sapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )})

## clean-up:
#for(i in t.vars){
#  del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
#}
