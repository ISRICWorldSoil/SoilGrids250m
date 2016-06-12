## Fit models for TAXOUSDA and generate predictions - SoilGrids250m
## By Tom.Hengl@isric.org with help from Marvin N. Wright <wright at imbs.uni-luebeck.de>

setwd("/data/models/TAXOUSDA")
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
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../saveRDS_functions.R")
source("../wrapper.predict_cs.R")

load("../equi7t3.rda")
load("../equi7t1.rda")
des <- read.csv("../SoilGrids250m_COVS250m.csv")

## class definitions:
col.legend <- read.csv("TAXOUSDA_legend.csv")
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXOUSDA.txt")
makeSAGAlegend(x=as.factor(paste(col.legend$Group)), MINIMUM=col.legend$Number, MAXIMUM=col.legend$Number+1, col_pal=col.legend$COLOR, filename="TAXOUSDA.txt")

load("../../profs/TAXOUSDA/TAXOUSDA.pnts.rda")
## 58,124 points!
## Post-processing filter:
soil.clim <- read.csv("../SOIL_CLIMATE_MATRIX.csv")
soil.clim <- soil.clim[soil.clim$Classification_system=="TAXOUSDA",-which(names(soil.clim) %in% c("Classification_system","COUNT_training","Min_lat","Max_lat","Max_elevation"))]
soil.fix <- data.frame(t(soil.clim[,-1]))
names(soil.fix) = gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", soil.clim$Name)))
soil.fix <- lapply(soil.fix, function(i){grep("x",i)})
soil.fix <- soil.fix[sapply(soil.fix, function(i){length(i)>0})]

## OVERLAY AND FIT MODELS:
ov <- extract.equi7t3(x=TAXOUSDA.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=48)
## TAKES ca 10 MINS FOR 40k points
#str(ov)
## remove all NA values:
#NA.rows <- which(!complete.cases(ov))
#ov$LATWGS84 <- TAXOUSDA.pnts@coords[,2]
write.csv(ov, file="ov.TAXOUSDA_SoilGrids250m.csv")
unlink("ov.TAXOUSDA_SoilGrids250m.csv.gz")
gzip("ov.TAXOUSDA_SoilGrids250m.csv")
save(ov, file="ov.TAXOUSDA.rda")
summary(ov$TAXOUSDA.f)
#load("ov.TAXOUSDA.rda")

## ------------- MODEL FITTING -----------

pr.lst <- des$WORLDGRIDS_CODE
formulaString.USDA = as.formula(paste('TAXOUSDA.f ~ ', paste(pr.lst, collapse="+"))) 
## Exclude latitude otherwise artifacts are possible ' LATWGS84 + '
formulaString.USDA

## Use Caret package to optimize model fitting
## http://stackoverflow.com/questions/18705159/r-caret-nnet-package-in-multicore
## fitting takes ca 30-60 mins
Nsub <- 1.5e4
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ctrl <- trainControl("boot",number=5)
max.Mtry = round((length(all.vars(formulaString.USDA)[-1]))/3)
rf.tuneGrid <- expand.grid(mtry = seq(5,max.Mtry,by=5))
#m_TAXOUSDA <- nnet::multinom(formulaString.USDA, ov, MaxNWts = 19000)
mnetX_TAXOUSDA <- caret::train(formulaString.USDA, data=ov, method="multinom", trControl=ctrl, MaxNWts = 19000, na.action=na.omit) ## 20 minutes
## Optimize fitting of random forest:
t.mrfX <- caret::train(formulaString.USDA, data=ov[sample.int(nrow(ov), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid)
## Ranger package (https://github.com/imbs-hl/ranger)
mrfX_TAXOUSDA <- ranger::ranger(formulaString.USDA, ov[complete.cases(ov[,all.vars(formulaString.USDA)]),], importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, probability=TRUE) ## TAKES 10 minutes
## to reduce the size of output objects use: 'num.trees=291'
stopCluster(cl)

cat("Results of model fitting 'nnet / randomForest':\n", file="TAXOUSDA_resultsFit.txt")
cat("\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
cat(paste("Variable:", all.vars(formulaString.USDA)[1]), file="TAXOUSDA_resultsFit.txt", append=TRUE)
cat("\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
sink(file="TAXOUSDA_resultsFit.txt", append=TRUE, type="output")
print(mnetX_TAXOUSDA)
cat("\n Random forest model:", file="TAXOUSDA_resultsFit.txt", append=TRUE)
print(mrfX_TAXOUSDA)
cat("\n Variable importance:\n", file="TAXOUSDA_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_TAXOUSDA))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
sink()

## save objects in parallel:
saveRDS.gz(mnetX_TAXOUSDA, file="mnetX_TAXOUSDA.rds")
saveRDS.gz(mrfX_TAXOUSDA, file="mrfX_TAXOUSDA.rds")
save.image()

## ------------- PREDICTIONS -----------

## for weigths use prediction accuracy (assessed using OOB/boosting):
gm1.w <- max(mnetX_TAXOUSDA$results$Accuracy, na.rm=TRUE)
gm2.w <- 1-mrfX_TAXOUSDA$prediction.error
mnetX_TAXOUSDA_final <- mnetX_TAXOUSDA$finalModel
rm(mnetX_TAXOUSDA)
gc(); gc()

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXOUSDA*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx("^TAXOUSDA*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## test predictions for a sample area:
#system.time( wrapper.predict_c(i="NA_060_036", varn="TAXOUSDA", gm1=mnetX_TAXOUSDA_final, gm2=mrfX_TAXOUSDA, in.path="/data/covs1t", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix, gm1.w=gm1.w, gm2.w=gm2.w) )
## plot in GE:
#x <- readGDAL("/data/predicted/NA_060_036/TAXOUSDA_Xeralfs_NA_060_036.tif")
#kml(x, file.name="TAXOUSDA_Xeralfs_NA_060_036.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_NA_060_036.png")

## run all predictions in parallel
pr.dirs <- basename(dirname(list.files(path="/data/covs1t", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 16,561 dirs

## Split models otherwise too large in size:
mrfX_TAXOUSDA_final <- split_rf(mrfX_TAXOUSDA)
rm(mrfX_TAXOUSDA)
gc(); gc()

## First, predict randomForest in parallel
## TAKES 8 hours
## Split to minimize problems of object size
for(j in 1:length(mrfX_TAXOUSDA_final)){
  gm = mrfX_TAXOUSDA_final[[j]]
  sfInit(parallel=TRUE, cpus=48)
  sfExport("gm", "pr.dirs", "split_predict_c", "j")
  sfLibrary(ranger)
  x <- sfLapply(pr.dirs, fun=function(x){ try( split_predict_c(x, gm, in.path="/data/covs1t", out.path="/data/predicted", split_no=j, varn="TAXOUSDA") )  } )
  sfStop()
}

## Second, predict nnet model in parallel
## TAKES ca 4 hrs
#cpus = unclass(round((256-30)/(2*(object.size(mnetX_TAXOUSDA_final)/1e9))))
sfInit(parallel=TRUE, cpus=48)
sfExport("mnetX_TAXOUSDA_final", "pr.dirs", "predict_nnet")
sfLibrary(nnet)
x <- sfLapply(pr.dirs, fun=function(x){ try( predict_nnet(x, mnetX_TAXOUSDA_final, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXOUSDA") )  } )
sfStop()

## Finally, sum up all predictions and generate geotifs
## TAKES ca 2.5 hrs
## this will also remove all temporary files
lev = mnetX_TAXOUSDA_final$lev
sfInit(parallel=TRUE, cpus=48)
sfExport("pr.dirs", "sum_predictions", "gm1.w", "gm2.w", "soil.fix", "lev", "col.legend")
sfLibrary(rgdal)
sfLibrary(plyr)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( sum_predictions(x, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXOUSDA", gm1.w=gm1.w, gm2.w=gm2.w, col.legend=col.legend, soil.fix=soil.fix, lev=lev) )  } )
sfStop()

## most probable class fix:
sfInit(parallel=TRUE, cpus=46)
sfExport("pr.dirs", "most_probable_fix", "lev", "col.legend") ## pr.dirs
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(raster)
sfLibrary(sp)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( most_probable_fix(x, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXOUSDA", col.legend=col.legend, lev=lev) )  } )
sfStop()


## ------------- VISUALIZATION -----------

## world plot - overlay and plot points and maps:
xy.pnts <- join(ov[,c("SOURCEID","SOURCEDB","TAXOUSDA.f")], as.data.frame(TAXOUSDA.pnts[c("SOURCEID")]), type="left", match="first")
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(TAXOUSDA.pnts)
plotKML(xy.pnts["TAXOUSDA.f"], folder.name="USDA classifications", file.name="TAXOUSDA_observed.kml")
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
hist.sel <- grep("Hist", xy.pnts$TAXOUSDA.f)
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
source("../../cv/cv_functions.R")
## Remove classes with too little observations:
summary(ov$TAXOUSDA.f)
xs = summary(ov$TAXOUSDA.f, maxsum=length(levels(ov$TAXOUSDA.f)))
sel.levs = attr(xs, "names")[xs > 5]
ov$TAXOUSDA.f2 <- ov$TAXOUSDA.f
ov$TAXOUSDA.f2[which(!ov$TAXOUSDA.f %in% sel.levs)] <- NA
ov$TAXOUSDA.f2 <- droplevels(ov$TAXOUSDA.f2)

## TAKES CA 1hr
formulaString2.USDA = as.formula(paste('TAXOUSDA.f2 ~ ', paste(pr.lst, collapse="+")))
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
