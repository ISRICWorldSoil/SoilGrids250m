## Cross-validation example
## by: Tom.Hengl@isric.org, Maria.RuiperezGonzales@wur.nl and Gerard.Heuvelink@wur.nl

library(sp)
library(randomForest)
library(nnet)
library(plotKML)
library(GSIF)
library(plyr)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(h2o)
library(scales)
source("cv_functions.R")

## load the data
set.seed(42)
data(eberg)
data(eberg_grid)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
eberg_grid@data <- cbind(eberg_grid@data, eberg_spc@predicted@data)
## overlay points and grids:
ov <- over(eberg, eberg_grid)
m <- cbind(ov, eberg@data)

## clean-up target variable:
summary(m$TAXGRSC)
m$soiltype <- NA
for(i in levels(m$TAXGRSC)){
 sel <- grep(pattern=i, m$TAXGRSC)
 if(length(sel)>5){
  m$soiltype[sel] <- i
 }
}
## regression matrix:
m <- m[complete.cases(m[,1:(ncol(eberg_grid)+2)]),]
m$soiltype <- as.factor(m$soiltype)
cov.lst <- paste0("PC", 1:11)

## Cross-validation factor-type variable:
formulaString = as.formula(paste('soiltype ~ ', paste(cov.lst, collapse="+")))
test.CLASS <- cv_factor(formulaString, rmatrix=m, nfold=5, idcol="ID")
str(test.CLASS)
test.CLASS[["Classes"]]

## Cross-validation numeric soil property:
ovP <- over(eberg, eberg_grid)
mP <- cbind(ovP, eberg@data)
## clean-up target variable:
summary(m$SNDMHT_A)
mP <- na.omit(mP)
formulaStringP = as.formula(paste('SNDMHT_A ~ ', paste(cov.lst, collapse="+")))
test.prop <- cv_numeric(formulaStringP, rmatrix=mP, nfold=5, idcol="ID")
str(test.prop)

## with h2o software:
h2o.init(nthreads = -1)
test.prop <- cv_numeric(formulaStringP, rmatrix=mP, nfold=5, idcol="ID", h2o=TRUE)
str(test.prop)

## Edgeroi data set:
data(edgeroi)
edgeroi.spc <- join(edgeroi$sites, edgeroi$horizons, type='inner')
coordinates(edgeroi.spc) <- ~ LONGDA94 + LATGDA94
proj4string(edgeroi.spc) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
edgeroi.spc <- spTransform(edgeroi.spc, CRS("+init=epsg:28355"))
## load the 250 m grids:
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
## overlay points and grids:
ov2 <- over(edgeroi.spc, edgeroi.grids)
m2 <- cbind(ov2, edgeroi.spc@data)
m2$DEPTH <- m2$UHDICM + (m2$LHDICM - m2$UHDICM)/2
formulaStringP2 = ORCDRC ~ DEMSRT5+TWISRT5+PMTGEO5+EV1MOD5+EV2MOD5+EV3MOD5+DEPTH
mP2 <- m2[complete.cases(m2[,all.vars(formulaStringP2)]),]

h2o.init(nthreads = -1)
test.ORC <- cv_numeric(formulaStringP2, rmatrix=mP2, nfold=5, idcol="SOURCEID", h2o=TRUE, Log=TRUE)
str(test.ORC)
## Plot CV results (use log-scale):
plt0 <- xyplot(test.ORC[[1]]$Predicted~test.ORC[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (machine learning)")
plt0 <- plt0 + layer(panel.abline(0,1,lty=1,lw=2,col="black"))
plt0

## Hexbin plot
d.meas <- min(test.ORC[[1]]$Observed, na.rm=TRUE)
pred <- test.ORC[[1]]$Predicted+ifelse(d.meas==0, 1, d.meas)
meas <- test.ORC[[1]]$Observed+ifelse(d.meas==0, 1, d.meas)
lim <- range(test.ORC[[1]]$Observed, na.rm=TRUE)
## Correlation plot:
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.abline(0,1,lty=1,lw=2,col="black")
  ## To plot RMSE around the 1:1 line:
  #panel.abline(0+test.ORC$Summary$logRMSE,1,lty=2,lw=2,col="black")
  #panel.abline(0-test.ORC$Summary$logRMSE,1,lty=2,lw=2,col="black")
}
plt <- hexbinplot(pred~meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main="Organic carbon in g/kg (accuracy assessment)", xlab="measured", ylab="predicted (ensemble)", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=25, density=40, xlim=lim, ylim=lim, panel=pfun)
plt

h2o.shutdown()