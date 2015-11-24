## Fit models for TAXOUSDA and generate predictions - SoilGrids250m
## By Tom.Hengl@isric.org

setwd("/data/models/TAXOUSDA")
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
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("../extract.equi7t3.R")
source("../wrapper.predict_c.R")
load("../equi7t3.rda")

## class definitions:
col.legend <- read.csv("TAXOUSDA_legend.csv")
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
#makeSAGAlegend(x=as.factor(paste(col.legend$Group)), MINIMUM=col.legend$Number, MAXIMUM=col.legend$Number+1, col_pal=col.legend$COLOR, filename="TAXOUSDA.txt")
des <- read.csv("../SoilGrids250m_COVS250m.csv")
load("../../profs/TAXOUSDA/TAXOUSDA.pnts.rda")
## 54,311 points!

## OVERLAY AND FIT MODELS:
ov <- extract.equi7t3(x=TAXOUSDA.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40)
## TAKES ca 10 MINS FOR 40k points
#str(ov)
## remove all NA values:
for(i in des$WORLDGRIDS_CODE){ ov[,i] <- ifelse(ov[,i]<= -32767, NA, ov[,i])  }
#NA.rows <- which(!complete.cases(ov))
ov$ES2MOD5 <- ifelse(is.na(ov$ES2MOD5)|ov$ES2MOD5<0, ov$ES3MOD5, ov$ES2MOD5)
ov$ES1MOD5 <- ifelse(is.na(ov$ES1MOD5)|ov$ES1MOD5<0, ov$ES2MOD5, ov$ES1MOD5) 
ov$ES6MOD5 <- ifelse(is.na(ov$ES6MOD5)|ov$ES6MOD5<0, ov$ES5MOD5, ov$ES6MOD5)
ov$EX1MOD5 <- ifelse(is.na(ov$EX1MOD5), ov$EX2MOD5, ov$EX1MOD5) 
ov$EX6MOD5 <- ifelse(is.na(ov$EX6MOD5), ov$EX5MOD5, ov$EX6MOD5)
#ov$LATWGS84 <- TAXOUSDA.pnts@coords[,2]

write.csv(ov, file="ov.TAXOUSDA_SoilGrids250m.csv")
unlink("ov.TAXOUSDA_SoilGrids250m.csv.gz")
gzip("ov.TAXOUSDA_SoilGrids250m.csv")
save(ov, file="ov.TAXOUSDA.rda")
#save(ov, file="../ov.rda")
summary(ov$TAXOUSDA.f)

## FIT MODELS:
pr.lst <- des$WORLDGRIDS_CODE
formulaString.USDA = as.formula(paste('TAXOUSDA.f ~ LATWGS84 + ', paste(pr.lst, collapse="+")))
formulaString.USDA

## TAKES > 20 mins to fit...
m_TAXOUSDA <- nnet::multinom(formulaString.USDA, ov, MaxNWts = 11000)
# group ‘Gelands’ is empty
#str(fitted(m_TAXOUSDA))
#head(signif(fitted(m_TAXOUSDA),3))
## goodness of fit:
cout.m <- as.factor(paste(predict(m_TAXOUSDA, newdata=ov, na.action = na.pass)))
cf <- mda::confusion(cout.m, as.character(ov[,"TAXOUSDA.f"]))
## remove missing classes:
a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
c.kappa = psych::cohen.kappa(cf[a,b])
ac <- sum(diag(cf))/sum(cf)*100
message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
## 35%
message(paste("Map purity:", signif(ac, 3)))
## 34%
saveRDS(m_TAXOUSDA, file="m_TAXOUSDA.rds")

## subset to complete pairs:
ovA <- ov[-m_TAXOUSDA$na.action,]
## 37973 points
ovA$TAXOUSDA.f <- as.factor(paste(ovA$TAXOUSDA.f))

## run in parallel?
## 2nd model:
mrf_TAXOUSDA <- randomForest(formulaString.USDA, data=ovA)
## validate the model:
#s <- sample.int(nrow(ovA), 5000)
#cv.mrf_TAXOUSDA <- randomForest(x=ovA[-s,all.vars(formulaString.USDA)[-1]], y=ovA$TAXOUSDA.f[-s], xtest=ovA[s,all.vars(formulaString.USDA)[-1]], ytest=ovA$TAXOUSDA.f[s])
#cv.mrf_TAXOUSDA$test$confusion[,"class.error"]
## http://stats.stackexchange.com/questions/30691/how-to-interpret-oob-and-confusion-matrix-for-random-forest

## 3rd model (Support Vector Machines):
#svm_TAXOUSDA <- svm(formulaString.USDA, data=ovA, probability=TRUE, cross=5)
#svm_TAXOUSDA$tot.accuracy
## TAKES >2 hrs!!

#mrf_TAXOUSDA 
## OOB estimate of error rate: 55.33%
varImpPlot(mrf_TAXOUSDA)
#object.size(mrf_TAXOUSDA)
saveRDS(mrf_TAXOUSDA, file="mrf_TAXOUSDA.rds")

#ss <- sample(1:nrow(ov), size=1000)
#mkk_TAXOUSDA <- kknn(formulaString.USDA, train=ov, test=ov, distance=1, kernel="triangular")

## predict for sample locations:
#wrapper.predict_c(i="NA_060_036", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend) ## gm3=svm_TAXOUSDA,
#wrapper.predict_c(i="EU_036_015", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend)
#wrapper.predict_c(i="EU_030_015", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend)
#wrapper.predict_c(i="OC_087_063", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend)
#wrapper.predict_c(i="EU_051_012", varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend)
## plot in GE:
#x <- readGDAL("/data/predicted/NA_060_036/TAXOUSDA_Xeralfs_NA_060_036.tif")
#kml(x, file.name="TAXOUSDA_Xeralfs_NA_060_036.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_NA_060_036.png")
#x <- readGDAL("/data/predicted/OC_087_063/TAXOUSDA_Xeralfs_OC_087_063.tif")
#kml(x, file.name="TAXOUSDA_Xeralfs_OC_087_063.kml", folder.name="Xeralfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Xeralfs_OC_087_063.png")
#x <- readGDAL("/data/predicted/EU_051_012/TAXOUSDA_Aqualfs_EU_051_012.tif")
#kml(x, file.name="TAXOUSDA_Aqualfs_EU_051_012.kml", folder.name="Aqualfs", colour=band1, z.lim=c(0,60), colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name="TAXOUSDA_Aqualfs_EU_051_012.png")

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXOUSDA*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## run all predictions in parallel
## TH: TAKES ABOUT 8-10 HOURS
## THERE IS A PROBLEM WITH RAM (HENCE <25 CORES) -> THE mrf_TAXOUSDA object is LARGE
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
sfInit(parallel=TRUE, cpus=15)
sfExport("wrapper.predict_c", "mrf_TAXOUSDA", "m_TAXOUSDA", "pr.dirs", "col.legend")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(plyr)
sfLibrary(nnet)
sfLibrary(randomForest)
## must be 'sfClusterApplyLB' otherwise 2-3 times slower
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( wrapper.predict_c(i=x, varn="TAXOUSDA", gm1=m_TAXOUSDA, gm2=mrf_TAXOUSDA, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend) )  } )
sfStop()

## plot in Google Earth:
# tmp.dir <- c("NA_060_036","NA_063_036","AF_072_048","OC_087_063","AS_072_087","AS_048_003","EU_051_012","SA_072_066","SA_090_048")
# for(i in tmp.dir){
#   r <- readGDAL(paste0("../predicted/", i, "/", ".tif"))
#   file.name <- paste0("TAXOUSDA_", i, ".kml")
#   kml(r, folder.name=i, file.name=file.name, colour=band1, layer.name="TAXOUSDA", subfolder.name="SoilGrids: USDA Soil suborder", colour_scale=col.legend$COLOR, raster_name=paste0("./TAXOUSDA_", i, ".png"))
# }

