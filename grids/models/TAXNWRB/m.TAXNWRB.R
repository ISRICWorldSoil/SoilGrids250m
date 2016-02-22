## Fit models for TAXNWRB and generate predictions - SoilGrids250m
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
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(plotKML)
library(GSIF)
library(maps)
library(scales)
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_c.R")
load("../equi7t3.rda")

## class definitions:
col.legend <- read.csv("TAXNWRB_legend.csv")
## 107 classes
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXNWRB.txt")
makeSAGAlegend(x=as.factor(as.character(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
des <- read.csv("../SoilGrids250m_COVS250m.csv")
#load("../cov.lst.rda") ## list of all Geotifs in 'covs' dir
load("../../profs/TAXNWRB/TAXNWRB.pnts.rda")
#str(TAXNWRB.pnts)
## 64623 points
## Post-processing filter:
soil.clim <- read.csv("../SOIL_CLIMATE_MATRIX.csv")
soil.clim <- soil.clim[soil.clim$Classification_system=="TAXNWRB",-which(names(soil.clim) %in% c("Classification_system","COUNT_training","Min_lat","Max_lat","Max_elevation"))]
soil.fix <- data.frame(t(soil.clim[,-1]))
names(soil.fix) = gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", soil.clim$Name)))
soil.fix <- lapply(soil.fix, function(i){grep("x",i)})
soil.fix <- soil.fix[sapply(soil.fix, function(i){length(i)>0})]

## OVERLAY (ca 20+ mins):
ov <- extract.equi7t3(x=TAXNWRB.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=48) 
## TAKES ca 10 MINS FOR 40k points
str(ov)
## 64623 obs. of  174 variables
#ov$LATWGS84 <- TAXNWRB.pnts@coords[,2]
## remove all NA values:
for(i in levels(des$WORLDGRIDS_CODE)){ ov[,i] <- ifelse(ov[,i]<= -10000, NA, ov[,i])  }
write.csv(ov, file="ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB_SoilGrids250m.csv.gz")
gzip("ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB.rda")
save(ov, file="ov.TAXNWRB.rda")
summary(ov$TAXNWRB.f)

## ------------- MODEL FITTING -----------

pr.lst <- des$WORLDGRIDS_CODE
formulaString.FAO = as.formula(paste('TAXNWRB.f ~ LATWGS84 + ', paste(pr.lst, collapse="+")))
formulaString.FAO

## TAKES > 20 mins to fit...
m_TAXNWRB <- nnet::multinom(formulaString.FAO, ov, MaxNWts = 19000)
## goodness of fit:
cout.m <- as.factor(paste(predict(m_TAXNWRB, newdata=ov, na.action = na.pass)))
cf <- mda::confusion(cout.m, as.character(ov[,"TAXNWRB.f"]))
## remove missing classes:
a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
c.kappa = psych::cohen.kappa(cf[a,b])
ac <- sum(diag(cf))/sum(cf)*100
message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
## 32%
message(paste("Map purity:", signif(ac, 3)))
## 23%
saveRDS(m_TAXNWRB, file="m_TAXNWRB.rds")

## subset to complete pairs:
ovA <- ov[-m_TAXNWRB$na.action,]
nrow(ovA)
ovA$TAXNWRB.f <- as.factor(paste(ovA$TAXNWRB.f))
## 2nd model:
mrf_TAXNWRB <- randomForest(formulaString.FAO, data=ovA)
#mrf_TAXNWRB 
## OOB estimate of error rate: ??%
varImpPlot(mrf_TAXNWRB)
#object.size(mrf_TAXNWRB)
saveRDS(mrf_TAXNWRB, file="mrf_TAXNWRB.rds")

## validate model:
#s <- sample.int(nrow(ovA), 5000)
#cv.mrf_TAXNWRB <- randomForest(x=ovA[-s,all.vars(formulaString.FAO)[-1]], y=ovA$TAXNWRB.f[-s], xtest=ovA[s,all.vars(formulaString.FAO)[-1]], ytest=ovA$TAXNWRB.f[s])
#cv.mrf_TAXNWRB$test$confusion[,"class.error"]
## http://stats.stackexchange.com/questions/30691/how-to-interpret-oob-and-confusion-matrix-for-random-forest

## ------------- PREDICTIONS -----------

## predict for sample locations:
wrapper.predict_c(i="NA_060_036", varn="TAXNWRB", gm1=m_TAXNWRB, gm2=mrf_TAXNWRB, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix)

## TH: TAKES ABOUT 12-14 HOURS
## PROBLEMS WITH RAM (HENCE <25 CORES) BECAUSE THE prediction location/data objects are LARGE
## 'Error in unserialize(socklist[[n]]) : error reading from connection'
## http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/#Memory_load

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXNWRB*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2353 dirs
sfInit(parallel=TRUE, cpus=15)
sfExport("wrapper.predict_c", "predict_df", "m_TAXNWRB", "mrf_TAXNWRB", "pr.dirs", "col.legend", "soil.fix")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(plyr)
sfLibrary(nnet)
sfLibrary(randomForest)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( wrapper.predict_c(i=x, varn="TAXNWRB", gm1=m_TAXNWRB, gm2=mrf_TAXNWRB, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix) )  } )
sfStop()

## world plot - overlay and plot points and maps:
xy.pnts <- join(ovA[,c("SOURCEID","SOURCEDB","TAXNWRB.f")], as.data.frame(TAXNWRB.pnts[c("SOURCEID")]), type="left", match="first")
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(TAXNWRB.pnts)
plotKML(xy.pnts["TAXNWRB.f"], folder.name="WRB classifications", file.name="TAXNWRB_observed.kml")
zip(zipfile="TAXNWRB_observed.kmz", files="TAXNWRB_observed.kml", zip="zip")

require(maptools)
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

## Cross-validation 10-fold:
source("../../cv/cv_functions.R")
## TAKES CA 1hr
test.WRB <- cv_factor(formulaString.FAO, ovA, nfold=10, idcol="SOURCEID")
str(test.WRB)
test.WRB[["Cohen.Kappa"]]
test.WRB[["Classes"]]
save(test.WRB, file="test.WRB.rda")
write.csv(test.WRB[["Classes"]], "cv_TAXNWRB_classes.csv")
write.csv(test.WRB[["Observed"]], "cv_TAXNWRB_observed.csv")
gzip("cv_TAXNWRB_observed.csv")
write.csv(test.WRB[["Predicted"]], "cv_TAXNWRB_predicted.csv")
gzip("cv_TAXNWRB_predicted.csv")
