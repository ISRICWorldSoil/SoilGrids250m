## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
library(snowfall)
library(utils)
library(plotKML)
library(GSIF)
#library(randomForestSRC)
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_n.R")
load("../equi7t3.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## points:
#load("../../profs/SPROPS/SPROPS.pnts.rda")
#load("../../profs/SPROPS/all.pnts.rda")
#ov <- extract.equi7t3(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
#str(ov)
#ovA <- join(all.pnts, ov[,c("LOC_ID",des$WORLDGRIDS_CODE)], type="left", by="LOC_ID")
#write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
#unlink("ov.SPROPS_SoilGrids250m.csv.gz")
#gzip("ov.SPROPS_SoilGrids250m.csv")
#save(ovA, file="ovA.rda", compression_level="xz")
load("ovA.rda")
## 1.3GB

t.vars <- c("ORCDRC", "PHIHOX", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM")
summary(ovA$CECSUM)
summary(ovA$BLD)
z.min <-c(0,20,0,0,0,0,0,0)
z.max <-c(800,110,100,100,100,100,3500,2200)
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + DEPTH +', paste(pr.lst, collapse="+")))})
#all.vars(formulaString.lst[[1]])

#s <- sample.int(nrow(ovA), 25000)
#dfs <- ovA[s,all.vars(formulaString.ORCDRC)]
## linear model:
#formulaString.ORCDRC.log = as.formula(paste('log1p(ORCDRC) ~ LATWGS84 + DEPTH +', paste(pr.lst, collapse="+")))
#m_ORCDRC <- lm(formulaString.ORCDRC.log, data=dfs)
#summary(m_ORCDRC)$adj.r.squared
## 30% of var explained

## H2O package more suited for large data (http://www.r-bloggers.com/benchmarking-random-forest-implementations/)
cat("Results of 'h2o.randomForest':\n\n", file="resultsFit.txt")
mrfX_path <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(is.null(mrfX_path[[j]])){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsFit.txt", append=TRUE)
    cat("\n", file="resultsFit.txt", append=TRUE)
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    dfs.hex <- as.h2o(dfs[complete.cases(dfs),], conn = h2o.getConnection(), destination_frame = "dfs.hex")
    #str(dfs.hex@mutable$col_names)
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex, importance=TRUE) ## TAKES ONLY 2-3 MINS even with 1M points - very fast!!
    sink(file="resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    print(mrfX@model$variable_importances)
    cat("\n", file="resultsFit.txt", append=TRUE)
    sink()
    mrfX_path[[j]] = h2o.saveModel(mrfX, path="./", force=TRUE)
  }
}
names(mrfX_path) = t.vars

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
## test predictions:
wrapper.predict_n(i="NA_060_036", varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max)
#x <- sapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )})

## run on 20 threads with 2 threads per process
sfInit(parallel=TRUE, cpus=20)
sfExport("wrapper.predict_n", "mrfX_path", "pr.dirs", "t.vars", "z.min", "z.max")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(h2o)
x <- sfLapply(pr.dirs, fun=function(x){ try( wrapper.predict_n(i=x, varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )  } )
sfStop()
