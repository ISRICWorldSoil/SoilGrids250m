## Fit model for TAXOUSDA

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(dplyr)
library(snowfall)
library(rgdal)
library(nnet)
library(psych)

source("extract.equi7t3.R")
load("../covs/equi7t3.rda")
des <- read.csv("../covs/SoilGrids250m_COVS250m.csv")
load("TAXOUSDA.pnts.rda")
ov <- extract.equi7t3(x=TAXOUSDA.pnts, y=des$WORLDGRIDS_CODE, equi7t3=equi7t3, path="../covs", cpus=40)
str(ov)
## 37,785 profiles
write.csv(ov, file="ov.TAXOUSDA_SoilGrids250m.csv")
save(ov, file="ov.TAXOUSDA.rda")
pr.lst <- des$WORLDGRIDS_CODE
formulaString.USDA = as.formula(paste('TAXOUSDA.f ~ ', paste(pr.lst, collapse="+")))
formulaString.USDA
## TAKES > 20 mins to fit / can not be parallelized?
m_TAXOUSDA <- nnet::multinom(formulaString.USDA, ov, MaxNWts = 9000)
# groups Anthrepts Ustox are empty
str(fitted(m_TAXOUSDA))
## 32,201 points
head(signif(fitted(m_TAXOUSDA),3))
## goodness of fit:
cout.m <- as.factor(paste(predict(m_TAXOUSDA, newdata=ov, na.action = na.pass)))
cf <- mda::confusion(cout.m, as.character(ov[,"TAXOUSDA.f"]))
## remove missing classes:
a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
c.kappa = psych::cohen.kappa(cf[a,b])
ac <- sum(diag(cf))/sum(cf)*100
message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
message(paste("Map purity:", signif(ac, 3)))

## Alternative models??
#mrf_TAXOUSDA <- randomForestSRC::rfsrc(formulaString.USDA, ov)

## predict for sample locations:
load("m_TAXOUSDA.rda")
str(m_TAXOUSDA)
## create dirs:
dir.lst <- list.dirs("D:/SoilGrids250m/covs")[-1]
Sys.chmod(gsub("covs", "predicted", dir.lst))

wrapper.predict_c(i="NA_060_036", mg=m_TAXOUSDA, in.path="D:/SoilGrids250m/covs", out.path="D:/SoilGrids250m/predicted")
