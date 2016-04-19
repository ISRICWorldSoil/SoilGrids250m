## Fit models for TAXNWRB and generate predictions - SoilGrids250m
## Tom.Hengl@isric.org

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
col.legend <- read.csv("TAXNWRB_legend.csv")
## 107 classes
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXNWRB.txt")
makeSAGAlegend(x=as.factor(as.character(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
## training points:
load("../../profs/TAXNWRB/TAXNWRB.pnts.rda")
#str(TAXNWRB.pnts)
## 64,623 points
## Post-processing filter:
soil.clim <- read.csv("../SOIL_CLIMATE_MATRIX.csv")
soil.clim <- soil.clim[soil.clim$Classification_system=="TAXNWRB",-which(names(soil.clim) %in% c("Classification_system","COUNT_training","Min_lat","Max_lat","Max_elevation"))]
soil.fix <- data.frame(t(soil.clim[,-1]))
names(soil.fix) = gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", soil.clim$Name)))
soil.fix <- lapply(soil.fix, function(i){grep("x",i)})
soil.fix <- soil.fix[sapply(soil.fix, function(i){length(i)>0})]

## spatial overlay (takes ca 20+ mins):
ov <- extract.equi7(x=TAXNWRB.pnts, y=des$WORLDGRIDS_CODE, equi7=equi7t3, path="/data/covs", cpus=48) 
str(ov)
## 64,623 obs. of  174 variables
write.csv(ov, file="ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB_SoilGrids250m.csv.gz")
gzip("ov.TAXNWRB_SoilGrids250m.csv")
unlink("ov.TAXNWRB.rda")
save(ov, file="ov.TAXNWRB.rda")
#load("ov.TAXNWRB.rda")
summary(ov$TAXNWRB.f)

## ------------- MODEL FITTING -----------

pr.lst <- des$WORLDGRIDS_CODE
formulaString.WRB = as.formula(paste('TAXNWRB.f ~ ', paste(pr.lst, collapse="+"))) 
## Exclude latitude otherwise artifacts are possible ' LATWGS84 + '
formulaString.WRB

## Use Caret package to optimize model fitting
## http://stackoverflow.com/questions/18705159/r-caret-nnet-package-in-multicore
## fitting takes ca 30-60 mins
Nsub <- 1.5e4
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ctrl <- trainControl("boot",number=5)
max.Mtry = round((length(all.vars(formulaString.WRB)[-1]))/3)
rf.tuneGrid <- expand.grid(mtry = seq(5,max.Mtry,by=5))
#m_TAXNWRB <- nnet::multinom(formulaString.WRB, ov, MaxNWts = 19000)
mnetX_TAXNWRB <- caret::train(formulaString.WRB, data=ov, method="multinom", trControl=ctrl, MaxNWts = 19000, na.action=na.omit) ## 10 minutes
## Optimize fitting of random forest:
t.mrfX <- caret::train(formulaString.WRB, data=ov[sample.int(nrow(ov), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid)
## Ranger package (https://github.com/imbs-hl/ranger)
mrfX_TAXNWRB <- ranger::ranger(formulaString.WRB, ov[complete.cases(ov[,all.vars(formulaString.WRB)]),], importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, probability=TRUE) ## TAKES 10 minutes
## to reduce the size of output objects use: 'num.trees=291'
stopCluster(cl)

cat("Results of model fitting 'nnet / randomForest':\n", file="TAXNWRB_resultsFit.txt")
cat("\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
cat(paste("Variable:", all.vars(formulaString.WRB)[1]), file="TAXNWRB_resultsFit.txt", append=TRUE)
cat("\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
sink(file="TAXNWRB_resultsFit.txt", append=TRUE, type="output")
print(mnetX_TAXNWRB)
cat("\n Random forest model:", file="TAXNWRB_resultsFit.txt", append=TRUE)
print(mrfX_TAXNWRB)
cat("\n Variable importance:\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
xl <- as.list(ranger::importance(mrfX_TAXNWRB))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:15]])))
sink()

## save objects in parallel:
saveRDS.gz(mnetX_TAXNWRB, file="mnetX_TAXNWRB.rds")
saveRDS.gz(mrfX_TAXNWRB, file="mrfX_TAXNWRB.rds")

## Alternative: Multinomial Logit Models R Package mnlogit (https://cran.r-project.org/web/packages/mnlogit/vignettes/mnlogit.pdf) unfortunatelly not easy to install and use
#m_TAXNWRB <- mnlogit::mnlogit(formulaString.FAO, ov, ncores=48, shape="wide")


## ------------- PREDICTIONS -----------

## for weigths use prediction accuracy (assessed using OOB/boosting):
gm1.w <- max(mnetX_TAXNWRB$results$Accuracy, na.rm=TRUE) 
gm2.w <- 1-mrfX_TAXNWRB$prediction.error
## We only need to export the final model from caret to make predictions
## (http://stackoverflow.com/questions/21096909/difference-betweeen-predictmodel-and-predictmodelfinalmodel-using-caret-for)
mnetX_TAXNWRB_final <- mnetX_TAXNWRB$finalModel
rm(mnetX_TAXNWRB)
gc(); gc()

## test predictions for a sample area:
#system.time( wrapper.predict_c(i="NA_061_036", varn="TAXNWRB", gm1=mnetX_TAXNWRB_final, gm2=mrfX_TAXNWRB, in.path="/data/covs1t", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix, mask_value=mask_value, gm1.w=gm1.w, gm2.w=gm2.w) )
#rm(gm1); rm(gm2)
#gc(); gc()

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXNWRB*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## run all predictions in parallel
pr.dirs <- basename(dirname(list.files(path="/data/covs1t", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 16,561 dirs

## Split models otherwise too large in size:
mrfX_TAXNWRB_final <- split_rf(mrfX_TAXNWRB)
rm(mrfX_TAXNWRB)
gc(); gc()

## First, predict randomForest in parallel
## TAKES 8 hours
## Split to minimize problems of object size
for(j in 1:length(mrfX_TAXNWRB_final)){
  gm = mrfX_TAXNWRB_final[[j]]
  cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
  sfInit(parallel=TRUE, cpus=cpus)
  sfExport("gm", "pr.dirs", "split_predict_c", "j")
  sfLibrary(ranger)
  x <- sfLapply(pr.dirs, fun=function(x){ try( split_predict_c(x, gm, in.path="/data/covs1t", out.path="/data/predicted", split_no=j, varn="TAXNWRB") )  } ) 
  sfStop()
}

## Second, predict nnet model in parallel
## TAKES ca 4 hrs
cpus = unclass(round((256-30)/(3.5*(object.size(mnetX_TAXNWRB_final)/1e9))))
sfInit(parallel=TRUE, cpus=cpus)
sfExport("mnetX_TAXNWRB_final", "pr.dirs", "predict_nnet")
sfLibrary(nnet)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( predict_nnet(x, mnetX_TAXNWRB_final, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXNWRB") )  } )
sfStop()

## Finally, sum up all predictions and generate geotifs
## TAKES ca 2.5 hrs
## this will also remove all temporary files
lev = mnetX_TAXNWRB_final$lev
sfInit(parallel=TRUE, cpus=35)
sfExport("pr.dirs", "sum_predictions", "gm1.w", "gm2.w", "soil.fix", "lev", "col.legend")
sfLibrary(rgdal)
sfLibrary(plyr)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( sum_predictions(x, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXNWRB", gm1.w=gm1.w, gm2.w=gm2.w, col.legend=col.legend, soil.fix=soil.fix, lev=lev) )  } )
sfStop()

## ------------- VISUALIZATION -----------

## world plot - overlay and plot points and maps:
xy.pnts <- join(ovA[,c("SOURCEID","SOURCEDB","TAXNWRB.f")], as.data.frame(TAXNWRB.pnts[c("SOURCEID")]), type="left", match="first")
coordinates(xy.pnts) = ~ LONWGS84+LATWGS84
proj4string(xy.pnts) = proj4string(TAXNWRB.pnts)
plotKML(xy.pnts["TAXNWRB.f"], folder.name="WRB classifications", file.name="TAXNWRB_observed.kml")
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
ovA <- ov[-m_TAXNWRB$na.action,]
nrow(ovA)
ovA$TAXNWRB.f <- as.factor(paste(ovA$TAXNWRB.f))
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
## Confusion matrix:
write.csv(test.WRB[["Confusion.matrix"]], "cv_TAXNWRB_Confusion.matrix.csv")
