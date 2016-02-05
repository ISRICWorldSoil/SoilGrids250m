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
library(h2o)
h2o.init(nthreads = -1)
test.prop <- cv_numeric(formulaStringP, rmatrix=mP, nfold=5, idcol="ID", h2o=TRUE)
str(test.prop)
plot(test.prop[[1]][,1:2], xlim=c(0,100), ylim=c(0,100), asp=1)
h2o.shutdown()
