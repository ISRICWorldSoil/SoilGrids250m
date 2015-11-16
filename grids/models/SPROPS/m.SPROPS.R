## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
library(randomForest)
library(nnet)
library(dplyr)
library(snowfall)
library(rgdal)
library(utils)
library(plotKML)
library(GSIF)
library(randomForestSRC)
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_n.R")
load("../equi7t3.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## points:
load("../../profs/SPROPS/SPROPS.pnts.rda")
load("../../profs/SPROPS/all.pnts.rda")
ov <- extract.equi7t3(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
#str(ov)
ovA <- join(all.pnts, ov[,c("LOC_ID",des$WORLDGRIDS_CODE)], type="left", by="LOC_ID")
write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
unlink("ov.SPROPS_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
## 1.3GB

## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.ORCDRC = as.formula(paste('ORCDRC ~ LATWGS84 + DEPTH +', paste(pr.lst, collapse="+")))
all.vars(formulaString.ORCDRC)

options(rf.cores=40) ## it seems that fitting of the rfsrc can not be parallelized
## subset:
s <- sample.int(nrow(ovA), 50000)
dfs <- ovA[s,]
mrf_ORCDR <- rfsrc(formulaString.ORCDRC, data=dfs[complete.cases(dfs[,all.vars(formulaString.ORCDRC)]),])
plot.rfsrc(mrf_ORCDR)
#mrf_ORCDR <- randomForest(formulaString.ORCDRC, data=ovA[!is.na(ovA$ORCDRC),])
#varImpPlot(mrf_ORCDR)
saveRDS(mrf_ORCDR, file="mrf_ORCDR.rds")

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
## test it:
wrapper.predict_n(i=varn="ORCDRC", gm=mrf_SPROPS, in.path="/data/covs", out.path="/data/predicted", z.min=0)

sfInit(parallel=TRUE, cpus=40)
sfExport("wrapper.predict_c", "mrf_ORCDR", "pr.dirs")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(randomForest)
sfLibrary(randomForestSRC)
x <- sfLapply(pr.dirs, fun=function(x){ try( wrapper.predict_n(i=x, varn="ORCDRC", gm=mrf_SPROPS, in.path="/data/covs", out.path="/data/predicted", z.min=0) )  } )
sfStop()


