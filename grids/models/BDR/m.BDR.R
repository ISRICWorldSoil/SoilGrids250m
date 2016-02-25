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
options(bitmapType='cairo')
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
ov <- extract.equi7t3(x=BDR.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=48) 
## TAKES CA 30 mins!
#str(ov)
ovA <- join(BDR_all.pnts, ov, type="left", by="LOC_ID")
write.csv(ovA, file="ov.BDR_SoilGrids250m.csv")
unlink("ov.BDR_SoilGrids250m.csv.gz")
gzip("ov.BDR_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
#load("ovA.rda")
## 1.3GB

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
z.min <- c(0,0,0)
z.max <- c(250,100,Inf)
pr.lst <- des$WORLDGRIDS_CODE
## remove some predictors that might lead to artifacts (buffer maps and land cover):
pr.lst <- pr.lst[-sapply(c("QUAUEA3","VOLNOA3","C??GLC5"), function(x){grep(
  glob2rx(x), pr.lst)})]
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + ', paste(pr.lst, collapse="+")))})
## For BDTICM we can not use latitude as predictor because points are somewhat clustered in LAT space 
formulaString.lst[[3]] = as.formula(paste(t.vars[3], ' ~ ', paste(pr.lst, collapse="+")))

## H2O package more suited for large data (http://www.r-bloggers.com/benchmarking-random-forest-implementations/)
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)

## Write results of model fitting into a text file:
cat("Results of 'h2o.randomForest':\n\n", file="resultsFit.txt")
mrfX_path <- rep(list(NULL), length(t.vars))
mdLX_path <- rep(list(NULL), length(t.vars))
#mrfX_path[[2]] = "/data/models/BDR/DRF_model_R_1450371688480_15"
#mdLX_path[[2]] = "/data/models/BDR/DeepLearning_model_R_1450371688480_20"
for(j in 1:length(t.vars)){
  if(is.null(mrfX_path[[j]])){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsFit.txt", append=TRUE)
    cat("\n", file="resultsFit.txt", append=TRUE)
    LOC_ID <- ovA$LOC_ID
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    sel <- complete.cases(dfs)
    dfs <- dfs[sel,]
    dfs.hex <- as.h2o(dfs, destination_frame = "dfs.hex")
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex) 
    mdLX <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex)
    sink(file="resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    ## Top 15 covariates:
    print(mrfX@model$variable_importances[1:15,])
    cat("\n", file="resultsFit.txt", append=TRUE)
    print(mdLX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    sink()
    mrfX_path[[j]] = h2o.saveModel(mrfX, path="./", force=TRUE)
    mdLX_path[[j]] = h2o.saveModel(mdLX, path="./", force=TRUE)
    ## save fitting success:
    fit.df <- data.frame(LOC_ID=LOC_ID[sel], observed=dfs[,1], predicted=as.data.frame(h2o.predict(mrfX, dfs.hex, na.action=na.pass))$predict)
    unlink(paste0("RF_fit_", t.vars[j], ".csv.gz"))
    write.csv(fit.df, paste0("RF_fit_", t.vars[j], ".csv"))
    gzip(paste0("RF_fit_", t.vars[j], ".csv"))
  }
}
names(mrfX_path) = t.vars
names(mdLX_path) = t.vars
write.table(mrfX_path, file="mrfX_path.txt")
write.table(mdLX_path, file="mdLX_path.txt")
#mrfX_path = as.list(read.table("mrfX_path.txt"))
#mdLX_path = as.list(read.table("mdLX_path.txt"))

## world plot - overlay and plot points and maps:
shp.pnts <- ovA[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84",t.vars)]
shp.pnts <- xy.pnts[complete.cases(xy.pnts[,c("LONWGS84","LATWGS84")]),]
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(BDR.pnts)
unlink("Soil_depth_training_points.shp")
writeOGR(xy.pnts, "Soil_depth_training_points.shp", "Soil_depth_training_points", "ESRI Shapefile")
unlink("Soil_depth_training_points.zip")
system("7za a Soil_depth_training_points.zip Soil_depth_training_points.*")

## ------------- PREDICTIONS -----------

## Predict per tile:
#pr.dirs <- basename(dirname(list.files(path="/data/covs1km", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2353 dirs
## test predictions:
wrapper.predict_2D(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max)
#wrapper.predict_2D(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max)
## Run in loop (whole world):
#x <- lapply(pr.dirs, function(i){try( wrapper.predict_2D(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )})
x <- lapply(pr.dirs, function(i){try( wrapper.predict_2D(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )})

x <- readGDAL("/data/predicted/NA_060_036/BDRLOG_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDRLOG_M_NA_060_036.kml", folder.name="R horizon", colour=band1, z.lim=c(0,100), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="BDRLOG_M_NA_060_036.png")
x <- readGDAL("/data/predicted/NA_060_036/BDTICM_M_NA_060_036.tif")
x.ll <- reproject(x)
kml(x.ll, file.name="BDTICM_M_NA_060_036.kml", folder.name="Absolute depth in cm", colour=band1, colour_scale=SAGA_pal[[1]], raster_name="BDTICM_M_NA_060_036.png", z.lim=c(0,5000))

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

h2o.shutdown()

## clean-up:
# for(i in c("BDRICM", "BDRLOG", "BDTICM")){ 
#   del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
# for(i in c("BDRICM", "BDRLOG", "BDTICM")){ 
#   del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
