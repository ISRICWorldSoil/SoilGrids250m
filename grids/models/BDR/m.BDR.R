## Fit models for depth to bedrock and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org & shanggv@hotmail.com

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

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_2D.R")
load("../equi7t3.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## points:
load("../../profs/BDR/BDR.pnts.rda")
load("../../profs/BDR/BDR_all.pnts.rda")
#ov <- extract.equi7t3(x=BDR.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
## TAKES CA 30 mins!
#str(ov)
#ovA <- join(BDR_all.pnts, ov, type="left", by="LOC_ID")
#write.csv(ovA, file="ov.BDR_SoilGrids250m.csv")
#unlink("ov.BDR_SoilGrids250m.csv.gz")
#gzip("ov.BDR_SoilGrids250m.csv")
#save(ovA, file="ovA.rda", compression_level="xz")
load("ovA.rda")
## 1.3GB

## BDRICM = depth to bedrock until 250 cm (sensored data)
## BDRLOG = occurrence of R horizon 0 / 1
## BDTICM = absolute depth to bedrock
t.vars <- c("BDRICM", "BDRLOG", "BDTICM")
summary(ovA$BDRICM)
summary(ovA$BDTICM)
z.min <-c(0,0,0)
z.max <-c(180,100,NA)
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + ', paste(pr.lst, collapse="+")))})

## H2O package more suited for large data (http://www.r-bloggers.com/benchmarking-random-forest-implementations/)
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)

cat("Results of 'h2o.randomForest':\n\n", file="resultsFit.txt")
mrfX_path <- rep(list(NULL), length(t.vars))
mdLX_path <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(is.null(mrfX_path[[j]])){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsFit.txt", append=TRUE)
    cat("\n", file="resultsFit.txt", append=TRUE)
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    dfs.hex <- as.h2o(dfs[complete.cases(dfs),], conn = h2o.getConnection(), destination_frame = "dfs.hex")
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex, importance=TRUE) ## TAKES ONLY 2-3 MINS 
    mdLX <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex)
    sink(file="resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    print(mrfX@model$variable_importances)
    cat("\n", file="resultsFit.txt", append=TRUE)
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
## test predictions:
<<<<<<< HEAD
wrapper.predict_2D(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max)
## Run in loop (whole world):
x <- lapply(pr.dirs, function(i){try( wrapper.predict_2D(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )})
#wrapper.predict_2D(i="NA_096_036", varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max)

x <- readGDAL("/data/predicted1km/NA_060_036/BDRLOG_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDRLOG_M_NA_060_036.kml", folder.name="R horizon", colour=band1, z.lim=c(0,100), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="BDRLOG_M_NA_060_036.png")
x <- readGDAL("/data/predicted1km/NA_060_036/BDTICM_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDTICM_M_NA_060_036.kml", folder.name="Absolute depth in cm", colour=band1, colour_scale=SAGA_pal[[1]], raster_name="BDTICM_M_NA_060_036.png") ## z.lim=c(0,5000), 
=======
wrapper.predict_2D(i="NA_060_036", varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max)
wrapper.predict_2D(i="NA_096_036", varn=t.vars, gm_path=mrfX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max)

x <- readGDAL("/data/predicted/NA_096_036/BDRLOG_M_NA_096_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDRLOG_M_NA_096_036.kml", folder.name="R horizon", colour=band1, z.lim=c(0,100), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="BDRLOG_M_NA_096_036.png")
x <- readGDAL("/data/predicted/NA_096_036/BDTICM_M_NA_096_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDTICM_M_NA_096_036.kml", folder.name="Absolute depth in cm", colour=band1, z.lim=c(0,5000), colour_scale=SAGA_pal[[1]], raster_name="BDTICM_M_NA_096_036.png")
>>>>>>> origin/master

h2o.shutdown(localH2O)
