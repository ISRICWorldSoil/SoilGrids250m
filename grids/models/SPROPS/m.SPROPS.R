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
ov <- extract.equi7t3(x=SPROPS.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="/data/covs", cpus=40) 
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
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + DEPTH +', paste(pr.lst, collapse="+")))})
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
    dfs <- ovA[,all.vars(formulaString.lst[[j]])]
    dfs.hex <- as.h2o(dfs[complete.cases(dfs),], conn = h2o.getConnection(), destination_frame = "dfs.hex")
    #str(dfs.hex@mutable$col_names)
    mrfX <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex, importance=TRUE) 
    mdLX <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString.lst[[j]])), training_frame=dfs.hex)
    sink(file="resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    print(mrfX@model$variable_importances)
    print(mdLX)
    sink()
    mrfX_path[[j]] = h2o.saveModel(mrfX, path="./", force=TRUE)
    mdLX_path[[j]] = h2o.saveModel(mdLX, path="./", force=TRUE)
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
## TAKES 3 DAYS!

## 250m resolution:
system.time( wrapper.predict_n(i="NA_060_036", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )
#system.time( wrapper.predict_n(i="SA_087_057", varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )
## all tiles one by one:
x <- sapply(pr.dirs, function(i){try( wrapper.predict_n(i, varn=t.vars, gm_path1=mrfX_path, gm_path2=mdLX_path, in.path="/data/covs", out.path="/data/predicted", z.min=z.min, z.max=z.max) )})
## TAKES >2 WEEKS!!

## clean-up:
# for(i in c("BLD", "ORCDRC", "PHIHOX")){
#  del.lst <- list.files(path="/data/predicted1km", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

# for(i in t.vars){
#  del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#  unlink(del.lst)
# }

h2o.shutdown()

## Cross-validation 10-fold:
source("../../cv/cv_functions.R")
library(h2o)
h2o.init(nthreads = -1)

cat("Results of Cross-validation:\n\n", file="resultsCV.txt")
cv_lst <- rep(list(NULL), length(t.vars))
for(j in 1:length(t.vars)){
  if(!file.exists(paste0("CV_", t.vars[j], ".rda"))){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsCV.txt", append=TRUE)
    cat("\n", file="resultsCV.txt", append=TRUE)
    cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ovA, nfold=10, idcol="SOURCEID", h2o=TRUE)
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
plt.names <- c("Soil organic carbon (g/kg)", "Soil pH x 10 in H2O ", "Soil pH x 10 in KCl", "Coarse fragments volumetric in percent", "Soil texture fraction sand in percent", "Soil texture fraction silt in percent", "Soil texture fraction clay in percent", "Bulk density (fine earth) in kg / cubic-meter", "Cation exchange capacity in cmolc/kg") 
names(plt.names) = t.vars
breaks.lst <- list(c(0,5,10,seq(20,1000,length=47)), seq(3.5,9.5,length=50), seq(2.5,9.5,length=50), c(0,1,2,5,seq(8,100,length=46)), seq(0,100,length=50), seq(0,100,length=50), seq(0,100,length=50), seq(150,2500,length=50), c(0,1,2,5,seq(8,800,length=26)))
names(breaks.lst) = t.vars

plot_CV <- function(j, breaks, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)){
  out.file = paste0("plot_CV_", t.vars[j], ".png")
  if(!file.exists(out.file)){
    load(paste0("CV_", t.vars[j], ".rda"))
    assign("m", get(paste0("CV_", t.vars[j]))[[1]])
    png(file = out.file, res = 150, width=850, height=850, type="cairo")
    if(any(t.vars[j] %in% c("ORCDRC","CECSUM","CRFVOL"))){
      d.meas <- min(m$Observed, na.rm=TRUE)
      out <- m$Predicted+ifelse(d.meas==0, 1, d.meas)
      meas <- m$Observed+ifelse(d.meas==0, 1, d.meas)
      lim <- range(breaks)+ifelse(d.meas==0, 1, d.meas)
      plt <- hexbinplot(out~meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=plt.names[t.vars[j]], xlab="measured", ylab="predicted (SoilGrids250m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=25, density=40, xlim=lim, ylim=lim, panel=pfun, colorcut=colorcut)
    } else {
      lim <- range(breaks)
      plt <- hexbinplot(m$Predicted~m$Observed, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=plt.names[t.vars[j]], xlab="measured", ylab="predicted (SoilGrids250m)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, xlim=lim, ylim=lim, asp=1, xbins=25, density=40, panel=pfun, colorcut=colorcut)
    }
    print(plt)
    dev.off()
  }
}

for(j in 1:length(t.vars)){
  plot_CV(j, breaks.lst[[t.vars[j]]])
}

