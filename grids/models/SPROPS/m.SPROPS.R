## Fit models for soil properties and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
#library(randomForestSRC)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
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
ov <- extract.equi7t3(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=48) 
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

write.csv(ovA, file="ov.SPROPS_SoilGrids250m.csv")
unlink("ov.SPROPS_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_SoilGrids250m.csv")
save(ovA, file="ovA.rda", compression_level="xz")
#load("ovA.rda")
## 1.3GB

t.vars <- c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM")
lapply(ovA[,t.vars], quantile, probs=c(0.01,0.5,0.99), na.rm=TRUE)

z.min <-c(0,20,20,0,0,0,0,0,0)
z.max <-c(800,110,110,100,100,100,100,3500,2200)
## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + DEPTH.f +', paste(pr.lst, collapse="+")))})
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
    LOC_ID <- ovA$LOC_ID
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    sel <- complete.cases(dfs)
    dfs <- dfs[sel,]
    dfs.hex <- as.h2o(dfs, destination_frame = "dfs.hex")
    #str(dfs.hex@mutable$col_names)
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex) 
    mdLX <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex)
    sink(file="resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n", file="resultsFit.txt", append=TRUE)
    ## Top 15 covariates:
    print(mrfX@model$variable_importances[1:15,])
    print(mdLX)
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

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
## 1km resolution (10-15 times faster):
system.time( wrapper.predict_n(i="NA_075_066", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )
system.time( wrapper.predict_n(i="NA_063_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )
## Bulk density only:
#system.time( wrapper.predict_n(i="NA_063_036", varn=t.vars[8], gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min[8], z.max=z.max[8]) )
## Run all tiles one by one:
x <- lapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs1km", out.path="/data/predicted1km", z.min=z.min, z.max=z.max) )})
## TAKES 3-4 DAYS!

## 250m resolution:
system.time( wrapper.predict_n(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )
#system.time( wrapper.predict_n(i="SA_087_057", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )
## all tiles one by one:
#x <- sapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )})
## TAKES >2 WEEKS!!

## clean-up:
# for(i in t.vars){ ## c("BLD", "ORCDRC", "PHIHOX") 
#  del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

# for(i in t.vars){
#  del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }
h2o.shutdown()

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

## Cross-validation 10-fold (TH: this does not account for high spatial clustering!):
library(h2o)
h2o.init(nthreads = -1)

source("../../cv/cv_functions.R")
cat("Results of Cross-validation:\n\n", file="resultsCV.txt")
cv_lst <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(!file.exists(paste0("CV_", t.vars[j], ".rda"))){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsCV.txt", append=TRUE)
    cat("\n", file="resultsCV.txt", append=TRUE)
    cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ovA, nfold=10, idcol="SOURCEID", h2o=TRUE, Log=TRUE)
    sink(file="resultsCV.txt", append=TRUE, type="output")
    print(cv_lst[[j]]$Summary)
    cat("\n", file="resultsCV.txt", append=TRUE)
    sink()
    assign(paste0("CV_", t.vars[j]), cv_lst[[j]])
    save(list=paste0("CV_", t.vars[j]), file=paste0("CV_", t.vars[j], ".rda"))
  }
}
h2o.shutdown()

## correlation plots:
source("../plot_hexbin.R")
plt.names <- c("SOC in g/kg", "Soil pH x 10 in H2O", "Soil pH x 10 in KCl", "Coarse fragments in %vol", "Sand fraction in %", "Silt fraction in %", "Clay fraction in %", "Bulk density (FE) in kg / m3", "CEC soil in cmolc/kg") 
names(plt.names) = t.vars
breaks.lst <- list(c(0,5,10,seq(20,1000,length=47)), seq(2.5,9.5,length=50), seq(2.5,9.5,length=50), c(0,1,2,5,seq(8,100,length=46)), seq(0,100,length=50), seq(0,100,length=50), seq(0,100,length=50), seq(450,2200,length=50), c(0,1,2,5,seq(8,450,length=26)))
names(breaks.lst) = t.vars
plt.log <- c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
names(plt.log) = t.vars

for(j in 1:length(t.vars)){
  plot_hexbin(j, breaks.lst[[t.vars[j]]], main=plt.names[t.vars[j]], in.file=paste0("CV_", t.vars[j], ".rda"), log.plot=plt.log[t.vars[j]])
}

## export to shapefile:
#sel <- CV_CECSUM[[1]]$Predicted>40 & CV_CECSUM[[1]]$Predicted<90 & CV_CECSUM[[1]]$Observed>2 & CV_CECSUM[[1]]$Observed<40 
#summary(sel)
#CEC.res <- join(CV_CECSUM[[1]][sel,], all.pnts, type="left", by="SOURCEID")
#CEC.res <- CEC.res[!is.na(CEC.res$LONWGS84),]
#coordinates(CEC.res) <- ~LONWGS84+LATWGS84
#writeOGR(CEC.res, "CEC.res.shp", "CEC.res", "ESRI Shapefile")
#hist(CEC.res$DEPTH)
#zip("CEC.res.*", files="CEC.res.zip")
