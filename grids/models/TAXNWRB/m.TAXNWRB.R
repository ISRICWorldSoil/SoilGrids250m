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
library(snowfall)
library(rgdal)
library(utils)
library(R.utils)
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
col.legend <- read.csv("TAXNWRB_legend.csv")
## 107 classes
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXNWRB.txt")
<<<<<<< HEAD
makeSAGAlegend(x=as.factor(as.character(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
=======
makeSAGAlegend(x=as.factor(paste(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
>>>>>>> origin/master
des <- read.csv("../SoilGrids250m_COVS250m.csv")
#load("../cov.lst.rda") ## list of all Geotifs in 'covs' dir
load("../../profs/TAXNWRB/TAXNWRB.pnts.rda")
#str(TAXNWRB.pnts)
## 34915 points!
#load("../ov.rda") ## check if the values exist already

## OVERLAY (ca 20+ mins):
ov <- extract.equi7t3(x=TAXNWRB.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
## TAKES ca 10 MINS FOR 40k points
str(ov)
## 42121 x 162 vars
#ov$LATWGS84 <- TAXNWRB.pnts@coords[,2]
## remove all NA values:
for(i in des$WORLDGRIDS_CODE){ ov[,i] <- ifelse(ov[,i]<= -10000, NA, ov[,i])  }
write.csv(ov, file="ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB_SoilGrids250m.csv.gz")
gzip("ov.TAXNWRB_SoilGrids250m.csv")
save(ov, file="ov.TAXNWRB.rda")
summary(ov$TAXNWRB.f)

## FIT MODELS:
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
## 29%
message(paste("Map purity:", signif(ac, 3)))
## 22%
saveRDS(m_TAXNWRB, file="m_TAXNWRB.rds")

## subset to complete pairs:
ovA <- ov[-m_TAXNWRB$na.action,]
nrow(ovA)
## 34491 points
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

## predict for sample locations:
#wrapper.predict_c(i="NA_063_036", varn="TAXNWRB", gm1=m_TAXNWRB, gm2=mrf_TAXNWRB, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend)

## run all predictions in parallel
## TH: TAKES ABOUT 12-14 HOURS
## THERE IS A PROBLEM WITH RAM (HENCE <25 CORES) BECAUSE THE prediction location/data objects are LARGE
## 'Error in unserialize(socklist[[n]]) : error reading from connection'
## http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/#Memory_load
## 'Error in unserialize(socklist[[n]]) : error reading from connection'
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2356 dirs
sfInit(parallel=TRUE, cpus=15)
sfExport("wrapper.predict_c", "m_TAXNWRB", "mrf_TAXNWRB", "pr.dirs", "col.legend")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(plyr)
sfLibrary(nnet)
sfLibrary(randomForest)
## 'sfClusterApplyLB' would have been better but it breaks as soon as the available RAM runs out
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( wrapper.predict_c(i=x, varn="TAXNWRB", gm1=m_TAXNWRB, gm2=mrf_TAXNWRB, in.path="/data/covs", out.path="/data/predicted", col.legend=col.legend) )  } )
sfStop()

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXNWRB*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

