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
library(randomForest)
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
source("../wrapper.predict_c.R")
source("../saveRDS_functions.R")
load("../equi7t3.rda")
load("../equi7t1.rda")

## class definitions:
col.legend <- read.csv("TAXNWRB_legend.csv")
## 107 classes
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)
unlink("TAXNWRB.txt")
makeSAGAlegend(x=as.factor(as.character(col.legend$Group)), MINIMUM=1:nrow(col.legend), MAXIMUM=1:nrow(col.legend)+1, col_pal=col.legend$COLOR, filename="TAXNWRB.txt")
des <- read.csv("../SoilGrids250m_COVS250m.csv")
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
## correction values:
mask_value <- as.list(des$MASK_VALUE)
names(mask_value) = des$WORLDGRIDS_CODE

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
formulaString.FAO = as.formula(paste('TAXNWRB.f ~ ', paste(pr.lst, collapse="+")))  ## LATWGS84 +
formulaString.FAO

## Use Caret package to optimize model fitting
## http://stackoverflow.com/questions/18705159/r-caret-nnet-package-in-multicore
Nsub <- 1.5e4 
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ctrl <- trainControl("boot",number=5)
rf.tuneGrid <- expand.grid(mtry = seq(4,22,by=2))
#m_TAXNWRB <- nnet::multinom(formulaString.FAO, ov, MaxNWts = 19000)
mnetX_TAXNWRB <- caret::train(formulaString.FAO, data=ov, method="multinom", trControl=ctrl, MaxNWts = 19000, na.action=na.omit) ## >10 minutes
t.mrfX <- caret::train(formulaString.FAO, data=ov[sample.int(nrow(ov), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid)
## Ranger package (https://github.com/imbs-hl/ranger)
mrfX_TAXNWRB <- ranger::ranger(formulaString.FAO, ov[complete.cases(ov[,all.vars(formulaString.FAO)]),], importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, probability=TRUE, num.trees=291) ## 2-5 minutes
#mrfX_TAXNWRB <- randomForest::randomForest(importance=TRUE, na.action=na.omit, mtry=t.mrfX$bestTune$mtry)
stopCluster(cl)

cat("Results of model fitting 'nnet / randomForest':\n", file="TAXNWRB_resultsFit.txt")
cat("\n", file="TAXNWRB_resultsFit.txt", append=TRUE)
cat(paste("Variable:", all.vars(formulaString.FAO)[1]), file="TAXNWRB_resultsFit.txt", append=TRUE)
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
system.time( wrapper.predict_c(i="NA_061_036", varn="TAXNWRB", gm1=mnetX_TAXNWRB_final, gm2=mrfX_TAXNWRB, in.path="/data/covs1t", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix, mask_value=mask_value, gm1.w=gm1.w, gm2.w=gm2.w) )
rm(gm1); rm(gm2)
gc(); gc()

## clean-up:
#del.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXNWRB*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

pr.dirs <- basename(dirname(list.files(path="/data/covs1t", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 16,561 dirs

## TAKES about 18 hrs
## There is some incompatibility / name overlap between parallel and snowfall:
detach("package:snowfall", unload=TRUE)
library(snowfall)
cpus = unclass(round((256-60)/(3.5*(object.size(mnetX_TAXNWRB_final)/1e9+object.size(mrfX_TAXNWRB)/1e9))))
## leave 2 CPU's free for background processes:
cpus <- ifelse(cpus>46, 46, cpus)
sfInit(parallel=TRUE, cpus=cpus)
sfExport("wrapper.predict_c", "predict_df", "pr.dirs", "col.legend", "soil.fix", "mnetX_TAXNWRB_final", "mrfX_TAXNWRB", "mask_value", "gm1.w", "gm2.w")
## export takes >5 mins
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(nnet)
#sfLibrary(caret)
sfLibrary(ranger)
x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( wrapper.predict_c(i=x, varn="TAXNWRB", gm1=mnetX_TAXNWRB_final, gm2=mrfX_TAXNWRB, in.path="/data/covs1t", out.path="/data/predicted", col.legend=col.legend, soil.fix=soil.fix, mask_value=mask_value, gm1.w=gm1.w, gm2.w=gm2.w) )  } )
sfStop()

## corrupt or missing tiles:
val.lst <- list.files(path="/data/predicted", pattern=glob2rx("^TAXNWRB_*_*_*.tif$"), full.names=TRUE, recursive=TRUE)
missing.lst <- pr.dirs[which(!pr.dirs %in% basename(dirname(val.lst)))]
str(missing.lst)
m <- readRDS(paste0("/data/covs1t/", missing.lst[1], "/", missing.lst[1], ".rds"))
str(m@data)

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
