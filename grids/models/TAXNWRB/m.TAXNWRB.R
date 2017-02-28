## Fit models for TAXNWRB and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "ranger", "mda", "psych", "stringr", "caret", "plotKML", "maptools", "maps")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/models/TAXNWRB")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
#devtools::install_github('topepo/caret/pkg/caret')
library(caret)
library(ranger)
library(nnet)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(plotKML)
library(GSIF)

#library(mnlogit)
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
col.legend <- read.csv("TAXNWRB_legend.csv")
## 118 classes
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXNWRB.txt")
makeSAGAlegend(x=as.factor(as.character(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
## training points:
load("/data/profs/TAXNWRB/TAXNWRB.pnts.rda")
#str(TAXNWRB.pnts)
## 70,560 points
## Post-processing filter:
soil.clim <- read.csv("/data/models/SOIL_CLIMATE_MATRIX.csv")
soil.clim <- soil.clim[soil.clim$Classification_system=="TAXNWRB",-which(names(soil.clim) %in% c("Classification_system","COUNT_training","Min_lat","Max_lat","Max_elevation"))]
soil.fix <- data.frame(t(soil.clim[,-1]))
names(soil.fix) = gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", soil.clim$Name)))
soil.fix <- lapply(soil.fix, function(i){grep("x",i)})
soil.fix <- soil.fix[sapply(soil.fix, function(i){length(i)>0})]
## subset to existing classes:
soil.fix <- soil.fix[names(soil.fix) %in% gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.",levels(TAXNWRB.pnts$TAXNWRB.f))))]
str(soil.fix)
## 24 classes need fixing

## spatial overlay (takes ca 20+ mins):
tile.pol = rgdal::readOGR("/data/models/tiles_ll_100km.shp", "tiles_ll_100km")
ov <- extract.tiled(x=TAXNWRB.pnts, tile.pol=tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=56)
str(ov)
## 70,560 obs. of  198 variables
write.csv(ov, file="ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB_SoilGrids250m.csv.gz")
R.utils::gzip("ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB.rda")
save(ov, file="ov.TAXNWRB.rda")
#load("ov.TAXNWRB.rda")
summary(ov$TAXNWRB.f)

## ------------- MODEL FITTING -----------

pr.lst <- basename(list.files(path="/data/stacked250m", ".tif"))
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-unlist(sapply(c("QUAUEA3","LCEE10"), function(x){grep(x, pr.lst)}))]
formulaString.WRB = as.formula(paste('TAXNWRB.f ~ ', paste(pr.lst, collapse="+")))
#formulaString.WRB

## Use Caret package to optimize model fitting
## http://stackoverflow.com/questions/18705159/r-caret-nnet-package-in-multicore
## fitting takes ca 30-60 mins
Nsub <- 8e3
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ctrl <- trainControl("boot",number=5)
max.Mtry = round((length(all.vars(formulaString.WRB)[-1]))/3)
rf.tuneGrid <- expand.grid(mtry = seq(10,max.Mtry,by=5))
#mnetX_TAXNWRB <- caret::train(formulaString.WRB, data=ov, method="multinom", trControl=ctrl, MaxNWts = 19000, na.action=na.omit)
## Optimize fitting of random forest:
dsf = ov[sample.int(nrow(ov), Nsub),]
dsf = dsf[complete.cases(dsf[,all.vars(formulaString.WRB)]),]
dsf$TAXNWRB.f = droplevels(dsf$TAXNWRB.f)
t.mrfX <- caret::train(formulaString.WRB, data=dsf, method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
## Ranger package:
mrfX_TAXNWRB <- ranger::ranger(formulaString.WRB, ov[complete.cases(ov[,all.vars(formulaString.WRB)]),], importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, probability=TRUE, num.trees=85) ## TAKES 10 minutes
## 'num.trees' - to reduce the size of output objects 
stopCluster(cl)

cat("Results of model fitting 'nnet / randomForest':\n", file="TAXNWRB_resultsFit.txt")
cat("\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
cat(paste("Variable:", all.vars(formulaString.WRB)[1]), file="TAXNWRB_resultsFit.txt", append=TRUE)
cat("\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
sink(file="TAXNWRB_resultsFit.txt", append=TRUE, type="output")
#print(mnetX_TAXNWRB)
cat("\n Random forest model:", file="TAXNWRB_resultsFit.txt", append=TRUE)
print(mrfX_TAXNWRB)
cat("\n Variable importance:\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_TAXNWRB))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
sink()

## save objects in parallel:
#saveRDS.gz(mnetX_TAXNWRB, file="mnetX_TAXNWRB.rds")
saveRDS.gz(mrfX_TAXNWRB, file="mrfX_TAXNWRB.rds")
save.image()

## Alternative: Multinomial Logit Models R Package mnlogit (https://cran.r-project.org/web/packages/mnlogit/vignettes/mnlogit.pdf) unfortunatelly not easy to install and use
#m_TAXNWRB <- mnlogit::mnlogit(formulaString.FAO, ov, ncores=48, shape="wide")


## ------------- PREDICTIONS -----------

#gm1.w <- max(mnetX_TAXNWRB$results$Accuracy, na.rm=TRUE)
gm2.w <- 1-mrfX_TAXNWRB$prediction.error
lev <- mrfX_TAXNWRB$forest$levels

## clean-up:
#del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("^TAXNWRB_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## run all predictions in parallel
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 18,653

## test prediction:
#factor_predict_ranger(i="T38275", gm=mrfX_TAXNWRB, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXNWRB", col.legend=col.legend, soil.fix=soil.fix, check.names=TRUE)
#factor_predict_ranger(i="T40505", gm=mrfX_TAXNWRB, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXNWRB", col.legend=col.legend, soil.fix=soil.fix, check.names=TRUE)

## ranger model only:
## ca 10-14 hrs of computing
try( detach("package:snowfall", unload=TRUE), silent=TRUE)
try( detach("package:snow", unload=TRUE), silent=TRUE)
library(parallel)
library(ranger)
library(rgdal)
library(plyr)
cpus = unclass(round((500-35)/(3.5*(object.size(mrfX_TAXNWRB)/1e9))))
cl <- parallel::makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
x = parLapply(cl, pr.dirs, fun=function(x){ try( factor_predict_ranger(i=x, gm=mrfX_TAXNWRB, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/tt/SoilGrids250m/predicted250m", varn="TAXNWRB", col.legend=col.legend, soil.fix=soil.fix, check.names=TRUE)  ) } )
stopCluster(cl)
gc(); gc()

## Mosaick:
t.vars = paste0("TAXNWRB_", gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", lev))))
library(snowfall)
sfInit(parallel=TRUE, cpus=ifelse(length(t.vars)>45, 45, length(t.vars)))
sfExport("t.vars", "make_mosaick_ll", "metasd", "sel.metasd")
out <- sfClusterApplyLB(1:length(t.vars), function(x){ try( make_mosaick_ll(varn=t.vars[x], i=NULL, in.path="/data/tt/SoilGrids250m/predicted250m", ot="Byte", dstnodata=255, metadata=metasd[grep(t.vars[x], metasd$FileName), sel.metasd]) )})
sfStop()

## ------------- VISUALIZATION -----------

## world plot - overlay and plot points and maps:
xy.pnts <- join(ov[,c("SOURCEID","SOURCEDB","TAXNWRB.f")], as.data.frame(TAXNWRB.pnts[c("SOURCEID")]), type="left", match="first")
xy.pnts <- xy.pnts[!duplicated(xy.pnts$SOURCEID),]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(TAXNWRB.pnts)
lev.leg <- join(data.frame(Group=levels(xy.pnts$TAXNWRB.f)), col.legend[,c("Group","COLOR")], type="left")
plotKML(xy.pnts["TAXNWRB.f"], folder.name="WRB classifications", file.name="TAXNWRB_observed.kml", colour_scale=lev.leg$COLOR)
zip(zipfile="TAXNWRB_observed.kmz", files="TAXNWRB_observed.kml", zip="zip")

library(maps)
require(maptools)
library(scales)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
no.plt <- xy.pnts@coords[,2]>-65 & xy.pnts@coords[,2]<85
png(file = "Fig_global_distribution_TAXNWRB.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(xy.pnts[xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
points(xy.pnts[!xy.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
dev.off()

## world plot - organic vs histosols:
hist.sel <- grep("Hist", xy.pnts$TAXNWRB.f)
length(hist.sel)
png(file = "Fig_global_distribution_Histosols.png", res = 150, width = 2000, height = 900)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(xy.pnts[-c(hist.sel, which(xy.pnts$SOURCEDB=="Simulated"), which(!no.plt)),], pch=21, bg=alpha("yellow", 0.6), cex=.8, col="grey")
points(xy.pnts[c(hist.sel),], pch=21, bg=alpha("red", 0.6), cex=1.5, col="black")
dev.off()

## ------------- CROSS-VALIDATION -----------

## subset to complete pairs:
## Remove classes with too little observations:
summary(ov$TAXNWRB.f)
xs = summary(ov$TAXNWRB.f, maxsum=length(levels(ov$TAXNWRB.f)))
sel.levs = attr(xs, "names")[xs > 5]
ov$TAXNWRB.f2 <- ov$TAXNWRB.f
ov$TAXNWRB.f2[which(!ov$TAXNWRB.f %in% sel.levs)] <- NA
ov$TAXNWRB.f2 <- droplevels(ov$TAXNWRB.f2)

## TAKES CA 1hr
formulaString2.FAO = as.formula(paste('TAXNWRB.f2 ~ ', paste(pr.lst, collapse="+")))
test.WRB <- cv_factor(formulaString2.FAO, ov[complete.cases(ov[,all.vars(formulaString2.FAO)]),], nfold=10, idcol="SOURCEID")
str(test.WRB)
test.WRB[["Cohen.Kappa"]]
test.WRB[["Classes"]]
save(test.WRB, file="test.WRB.rda")
unlink("cv_TAXNWRB_classes.csv.gz")
unlink("cv_TAXNWRB_observed.csv.gz")
unlink("cv_TAXNWRB_predicted.csv.gz")
write.csv(test.WRB[["Classes"]], "cv_TAXNWRB_classes.csv")
write.csv(test.WRB[["Observed"]], "cv_TAXNWRB_observed.csv")
gzip("cv_TAXNWRB_observed.csv")
write.csv(test.WRB[["Predicted"]], "cv_TAXNWRB_predicted.csv")
gzip("cv_TAXNWRB_predicted.csv")
## Confusion matrix:
write.csv(test.WRB[["Confusion.matrix"]], "cv_TAXNWRB_Confusion.matrix.csv")
closeAllConnections()
