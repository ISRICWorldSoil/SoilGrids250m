## Fit models for TAXOUSDA and generate predictions - SoilGrids250m
## By Tom.Hengl@isric.org with help from Marvin N. Wright <wright at imbs.uni-luebeck.de>

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "ranger", "mda", "psych", "stringr", "caret", "plotKML", "maptools", "maps")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/models/TAXOUSDA")
load(".RData")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
#devtools::install_github('topepo/caret/pkg/caret')
library(caret)
#install_bitbucket("mkuhn/parallelRandomForest", ref="parallelRandomForest")
#library(parallelRandomForest)
#library(e1071)
#library(randomForest)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger", ref="forest_memory") ## version to deal with Memory problems
library(ranger)
library(nnet)
library(dplyr)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(plotKML)
library(GSIF)

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

## class definitions:
col.legend <- read.csv("TAXOUSDA_legend.csv")
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXOUSDA.txt")
makeSAGAlegend(x=as.factor(paste(col.legend$Group)), MINIMUM=col.legend$Number, MAXIMUM=col.legend$Number+1, col_pal=col.legend$COLOR, filename="TAXOUSDA.txt")

load("/data/profs/TAXOUSDA/TAXOUSDA.pnts.rda")
str(TAXOUSDA.pnts@data)
#summary(TAXOUSDA.pnts$TAXOUSDA.f)
## 342,745 points!
xs = summary(TAXOUSDA.pnts$TAXOUSDA.f, maxsum=length(levels(TAXOUSDA.pnts$TAXOUSDA.f)))
sel.levs = attr(xs, "names")[xs > 10]
TAXOUSDA.pnts$soiltype <- TAXOUSDA.pnts$TAXOUSDA.f
TAXOUSDA.pnts$soiltype[which(!TAXOUSDA.pnts$TAXOUSDA.f %in% sel.levs)] <- NA
TAXOUSDA.pnts$soiltype <- droplevels(TAXOUSDA.pnts$soiltype)
TAXOUSDA.pnts <- TAXOUSDA.pnts[!is.na(TAXOUSDA.pnts$soiltype),]
summary(TAXOUSDA.pnts$soiltype)

## Post-processing filter:
soil.clim <- read.csv("/data/models/SOIL_CLIMATE_MATRIX.csv")
soil.clim <- soil.clim[soil.clim$Classification_system=="TAXOUSDA",-which(names(soil.clim) %in% c("Classification_system","COUNT_training","Min_lat","Max_lat","Max_elevation"))]
soil.fix <- data.frame(t(soil.clim[,-1]))
names(soil.fix) = gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", soil.clim$Name)))
soil.fix <- lapply(soil.fix, function(i){grep("x",i)})
soil.fix <- soil.fix[sapply(soil.fix, function(i){length(i)>0})]
## subset to existing classes:
soil.fix <- soil.fix[names(soil.fix) %in% levels(TAXOUSDA.pnts$soiltype)]

## OVERLAY AND FIT MODELS:
tile.pol = rgdal::readOGR("/data/models/tiles_ll_100km.shp", "tiles_ll_100km")
#tile.pol = readRDS("/data/models/stacked250m_tiles_pol.rds")
ov <- extract.tiled(x=TAXOUSDA.pnts, tile.pol=tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=56)

## TAKES ca 10 MINS FOR 300k points
#str(ov)
write.csv(ov, file="ov.TAXOUSDA_SoilGrids250m.csv")
unlink("ov.TAXOUSDA_SoilGrids250m.csv.gz")
R.utils::gzip("ov.TAXOUSDA_SoilGrids250m.csv")
save(ov, file="ov.TAXOUSDA.rda")
#summary(ov$soiltype)
#load("ov.TAXOUSDA.rda")
save.image()

## ------------- MODEL FITTING -----------

pr.lst <- basename(list.files(path="/data/stacked250m", ".tif"))
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-unlist(sapply(c("QUAUEA3","LCEE10"), function(x){grep(x, pr.lst)}))]
formulaString.USDA = as.formula(paste('soiltype ~ ', paste(pr.lst, collapse="+")))
#formulaString.USDA

## Use Caret package to optimize model fitting
## http://stackoverflow.com/questions/18705159/r-caret-nnet-package-in-multicore
## fitting takes ca 30-60 mins
Nsub <- 8e3
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
max.Mtry = round((length(all.vars(formulaString.USDA)[-1]))/3)
rf.tuneGrid <- expand.grid(mtry = seq(10,max.Mtry,by=5))
#mnetX_TAXOUSDA <- caret::train(formulaString.USDA, data=ov, method="multinom", trControl=ctrl, MaxNWts = 19000, na.action=na.omit) ## ?? minutes
## Optimize fitting of random forest:
dsf = ov[sample.int(nrow(ov), Nsub),]
dsf = dsf[complete.cases(dsf[,all.vars(formulaString.USDA)]),]
dsf$soiltype = droplevels(dsf$soiltype)
t.mrfX <- caret::train(formulaString.USDA, data=dsf, method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
## Ranger package:
mrfX_TAXOUSDA <- ranger::ranger(formulaString.USDA, ov[complete.cases(ov[,all.vars(formulaString.USDA)]),], importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, probability=TRUE, num.trees=85) ## TAKES 10 minutes
## 'num.trees' - to reduce the size of output objects 
stopCluster(cl)

cat("Results of model fitting 'nnet / randomForest':\n", file="TAXOUSDA_resultsFit.txt")
cat("\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
cat(paste("Variable:", all.vars(formulaString.USDA)[1]), file="TAXOUSDA_resultsFit.txt", append=TRUE)
cat("\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
sink(file="TAXOUSDA_resultsFit.txt", append=TRUE, type="output")
#print(mnetX_TAXOUSDA)
cat("\n Random forest model:", file="TAXOUSDA_resultsFit.txt", append=TRUE)
print(mrfX_TAXOUSDA)
cat("\n Variable importance:\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_TAXOUSDA))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
sink()

## save objects in parallel:
#saveRDS.gz(mnetX_TAXOUSDA, file="mnetX_TAXOUSDA.rds")
saveRDS.gz(mrfX_TAXOUSDA, file="mrfX_TAXOUSDA.rds")
save.image()

## ------------- PREDICTIONS -----------

## for weigths use prediction accuracy (assessed using OOB/boosting):
#gm1.w <- max(mnetX_TAXOUSDA$results$Accuracy, na.rm=TRUE)
gm2.w <- 1-mrfX_TAXOUSDA$prediction.error
lev <- mrfX_TAXOUSDA$forest$levels


## clean-up:
#del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("^TAXOUSDA_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## run all predictions in parallel
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 18,653

## test prediction:
#factor_predict_ranger(i="T38275", gm=mrfX_TAXOUSDA, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXOUSDA", col.legend=col.legend, soil.fix=soil.fix)
#factor_predict_ranger(i="T40502", gm=mrfX_TAXOUSDA, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXOUSDA", col.legend=col.legend, soil.fix=soil.fix)

## ranger model only:
## ca 14 hrs of computing
try( detach("package:snowfall", unload=TRUE), silent=TRUE)
try( detach("package:snow", unload=TRUE), silent=TRUE)
library(parallel)
library(ranger)
library(rgdal)
library(plyr)
cpus = unclass(round((500-35)/(3.5*(object.size(mrfX_TAXOUSDA)/1e9))))
cl <- parallel::makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
x = parLapply(cl, pr.dirs, fun=function(x){ try( factor_predict_ranger(i=x, gm=mrfX_TAXOUSDA, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXOUSDA", col.legend=col.legend, soil.fix=soil.fix)  ) } )
stopCluster(cl)
gc(); gc()

## Mosaick:
t.vars = paste0("TAXOUSDA_", lev)
library(snowfall)
sfInit(parallel=TRUE, cpus=ifelse(length(t.vars)>45, 45, length(t.vars)))
sfExport("t.vars", "make_mosaick_ll", "metasd", "sel.metasd")
out <- sfClusterApplyLB(1:length(t.vars), function(x){ try( make_mosaick_ll(varn=t.vars[x], i=NULL, in.path="/data/tt/SoilGrids250m/predicted250m", ot="Byte", dstnodata=255, metadata=metasd[which(metasd$FileName == t.vars[x]), sel.metasd]) )})
sfStop()


## ------------- VISUALIZATION -----------

## world plot - overlay and plot points and maps:
xy.pnts <- join(ov[,c("SOURCEID","SOURCEDB","soiltype")], as.data.frame(TAXOUSDA.pnts[c("SOURCEID")]), type="left", match="first")
xy.pnts <- xy.pnts[!duplicated(xy.pnts$SOURCEID),]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(TAXOUSDA.pnts)
lev.leg <- join(data.frame(Group=levels(xy.pnts$soiltype)), col.legend[,c("Group","COLOR")], type="left")
plotKML(xy.pnts["soiltype"], folder.name="USDA classifications", file.name="TAXOUSDA_observed.kml", colour_scale=lev.leg$COLOR)
zip(zipfile="TAXOUSDA_observed.kmz", files="TAXOUSDA_observed.kml", zip="zip")

require(maptools)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
no.plt <- xy.pnts@coords[,2]>-65 & xy.pnts@coords[,2]<85
png(file = "Fig_global_distribution_TAXOUSDA.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(xy.pnts[xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
points(xy.pnts[!xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
dev.off()

## world plot - organic vs histosols:
hist.sel <- grep("Hist", xy.pnts$soiltype)
length(hist.sel)
png(file = "Fig_global_distribution_Histosols_USDA.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(xy.pnts[-c(hist.sel, which(xy.pnts$SOURCEDB=="Simulated"), which(!no.plt)),], pch=21, bg=alpha("yellow", 0.6), cex=.8, col="grey")
points(xy.pnts[c(hist.sel),], pch=21, bg=alpha("red", 0.6), cex=1.5, col="black")
dev.off()

## plot in Google Earth:
# tmp.dir <- c("NA_060_036","NA_063_036","AF_072_048","OC_087_063","AS_072_087","AS_048_003","EU_051_012","SA_072_066","SA_090_048")
# for(i in tmp.dir){
#   r <- readGDAL(paste0("../predicted/", i, "/", ".tif"))
#   file.name <- paste0("TAXOUSDA_", i, ".kml")
#   kml(r, folder.name=i, file.name=file.name, colour=band1, layer.name="TAXOUSDA", subfolder.name="SoilGrids: USDA Soil suborder", colour_scale=col.legend$COLOR, raster_name=paste0("./TAXOUSDA_", i, ".png"))
# }

## Cross-validation 10-fold:
source("/data/cv/cv_functions.R")

## TAKES CA 1hr
formulaString2.USDA = as.formula(paste('soiltype ~ ', paste(pr.lst, collapse="+")))
test.USDA <- cv_factor(formulaString2.USDA, ov[complete.cases(ov[,all.vars(formulaString2.USDA)]),], nfold=10, idcol="SOURCEID")
str(test.USDA)
test.USDA[["Cohen.Kappa"]]
test.USDA[["Classes"]]
save(test.USDA, file="test.USDA.rda")
unlink("cv_TAXOUSDA_classes.csv.gz")
unlink("cv_TAXOUSDA_observed.csv.gz")
unlink("cv_TAXOUSDA_predicted.csv.gz")
write.csv(test.USDA[["Classes"]], "cv_TAXOUSDA_classes.csv")
write.csv(test.USDA[["Observed"]], "cv_TAXOUSDA_observed.csv")
gzip("cv_TAXOUSDA_observed.csv")
write.csv(test.USDA[["Predicted"]], "cv_TAXOUSDA_predicted.csv")
gzip("cv_TAXOUSDA_predicted.csv")
