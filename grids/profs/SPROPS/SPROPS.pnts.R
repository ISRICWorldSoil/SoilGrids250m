## Prepare callibration points for soil organic carbon, pH, BLD, texture fractions and similar (SoilGrids250m)
## Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(GSIF)
library(rgdal)
library(lattice)
library(scales)
library(plotKML)
library(dplyr)

## list of input data sets:
tax.lst <- list.files(path="/data/soilstorage/SoilData", pattern=glob2rx("SPROPS.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})
## 29 data sets
names(in.lst) = basename(tax.lst)

## Simulated sand dunes points:
load("../TAXOUSDA/deserts.pnt.rda")
SPROPS.sim <- as.data.frame(spTransform(deserts.pnt, CRS("+proj=longlat +datum=WGS84")))
SPROPS.sim[,1] <- NULL
SPROPS.sim <- plyr::rename(SPROPS.sim, c("x"="LONWGS84", "y"="LATWGS84"))
SPROPS.sim$SOURCEID <- paste("SIM", 1:nrow(SPROPS.sim), sep="_")
SPROPS.sim$SOURCEDB = "Simulated"
## we assume - top-soil sub-soil are exactly the same
SPROPS.sim <- rbind(SPROPS.sim, SPROPS.sim)
##  http://www.jstor.org/stable/30055331
SPROPS.sim$ORCDRC = 0
SPROPS.sim$SNDPPT = 98
SPROPS.sim$SLTPPT = 2
SPROPS.sim$CLYPPT = 0
SPROPS.sim$CRFVOL = 2
SPROPS.sim$PHIHOX = 8.1
SPROPS.sim$UHDICM = c(rep(0, nrow(deserts.pnt)), rep(5, nrow(deserts.pnt))) 
SPROPS.sim$LHDICM = c(rep(5, nrow(deserts.pnt)), rep(200, nrow(deserts.pnt)))
SPROPS.sim$DEPTH = c(rep(2.5, nrow(deserts.pnt)), rep(97.5, nrow(deserts.pnt)))
str(SPROPS.sim)

## add to the list:
in.lst[[length(tax.lst)+1]] <- SPROPS.sim
## Bind everything together:
all.pnts <- dplyr::bind_rows(in.lst)
str(all.pnts)
## Warning message:
## In bind_rows_(x, .id) : Unequal factor levels: coercing to character
## 809818 obs. of  20 variables
all.pnts$LOC_ID <- as.factor(paste("ID", all.pnts$LONWGS84, all.pnts$LATWGS84, sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## 154,412 unique locations!
## Filter out remaining values outside of physical range:
all.pnts$SNDPPT <- ifelse(all.pnts$SNDPPT<0|all.pnts$SNDPPT>100, NA, all.pnts$SNDPPT)
hist(all.pnts$SNDPPT, col="blue", xlab="Sand in %", breaks=40, main=sum(!is.na(all.pnts$SNDPPT)))
all.pnts$SLTPPT <- ifelse(all.pnts$SLTPPT<0|all.pnts$SLTPPT>100, NA, all.pnts$SLTPPT)
hist(all.pnts$SLTPPT, col="blue", xlab="Silt in %", breaks=40, main=sum(!is.na(all.pnts$SLTPPT)))
all.pnts$CLYPPT <- ifelse(all.pnts$CLYPPT<0|all.pnts$CLYPPT>100, NA, all.pnts$CLYPPT)
hist(all.pnts$CLYPPT, col="blue", xlab="Clay in %", breaks=40, main=sum(!is.na(all.pnts$CLYPPT)))
all.pnts$PHIHOX <- ifelse(all.pnts$PHIHOX<2|all.pnts$PHIHOX>12, NA, all.pnts$PHIHOX)
hist(all.pnts$PHIHOX, col="blue", xlab="pH in H2O", breaks=40, main=sum(!is.na(all.pnts$PHIHOX)))
all.pnts$PHIKCL <- ifelse(all.pnts$PHIKCL<2|all.pnts$PHIKCL>12, NA, all.pnts$PHIKCL)
hist(all.pnts$PHIKCL, col="blue", xlab="pH in KCl", breaks=40, main=sum(!is.na(all.pnts$PHIKCL)))
all.pnts$ORCDRC <- ifelse(all.pnts$ORCDRC>800, NA, all.pnts$ORCDRC)
hist(log1p(all.pnts$ORCDRC), col="blue", xlab="log-Organic carbon", breaks=40, main=sum(!is.na(all.pnts$ORCDRC)))
all.pnts$CRFVOL <- ifelse(all.pnts$CRFVOL<0|all.pnts$CRFVOL>100, NA, all.pnts$CRFVOL)
## some groupings around very high values for BLD - an error?
all.pnts$BLD <- ifelse(all.pnts$BLD<50|all.pnts$BLD>2500, NA, all.pnts$BLD)
## BLD -> 630,372 missing values!
hist(all.pnts$BLD, col="blue", xlab="Bulk density", breaks=40, main=sum(!is.na(all.pnts$BLD)))
all.pnts$CECSUM <- ifelse(all.pnts$CECSUM<0, 0, ifelse(all.pnts$CECSUM>600, NA, all.pnts$CECSUM))
hist(log1p(all.pnts$CECSUM), col="blue", xlab="log-CEC", breaks=40, main=sum(!is.na(all.pnts$CECSUM)))
summary(all.pnts)
summary(as.factor(all.pnts$SOURCEDB))
## Scale the depths so they all start at 0 cm (soil surface) -> otherwise we might miss some important "O", "H" or similar topsoil / organic horizons;
## See also: https://groups.google.com/forum/#!topic/soilgrids-dev/zT0dy_XTuHQ
hist(all.pnts$UHDICM, breaks=30)
## There are some profiles that indicate depths beyond 10 m - probably not a good idea to include into modelling?
#View(all.pnts[which(all.pnts$UHDICM>1000),])
all.pnts <- all.pnts[all.pnts$UHDICM<600,]
hist(log1p(all.pnts$LHDICM), breaks=30)
## the following takes >15 mins:
z.min <- ddply(all.pnts, .(SOURCEID), summarize, aggregated = min(UHDICM, na.rm=TRUE))
z.shift <- join(all.pnts[,c("SOURCEID","SOURCEDB")], z.min, type="left")$aggregated
z.shift <- ifelse(z.shift>0, 0, z.shift)
summary(z.shift)
all.pnts$UHDICM.f <- all.pnts$UHDICM - z.shift
all.pnts$LHDICM.f <- all.pnts$LHDICM - z.shift
str(all.pnts[which(z.shift<0)[720],])
str(all.pnts[which(z.shift<0)[721],])
all.pnts$DEPTH.f <- all.pnts$UHDICM.f + (all.pnts$LHDICM.f - all.pnts$UHDICM.f)/2
summary(all.pnts$UHDICM.f, breaks=30)
## check distributions:
par(mar=c(11,5,3,1))
boxplot(log1p(ORCDRC)~SOURCEDB, all.pnts[all.pnts$LHDICM<30,], main="log1p(ORCDRC)", las=2)
boxplot(PHIHOX~SOURCEDB, all.pnts, main="PHIHOX", las=2)
## Organic carbon as function of depth:
OC.m <- lm(log1p(ORCDRC)~log1p(DEPTH.f), all.pnts)
summary(OC.m)
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.abline(coef(OC.m)[1],coef(OC.m)[2],lty=1,lw=2,col="black")
}
pal = R_pal[["bpy_colors"]][1:18]
library(hexbin)
hexbinplot(OC.m$model$`log1p(ORCDRC)`~OC.m$model$`log1p(DEPTH.f)`, colramp=colorRampPalette(pal), ylab="ORCDRC", xlab="depth", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, asp=1, xbins=30, ybins=30, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) ## scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE))
## save data file:
save(all.pnts, file="all.pnts.rda")

SPROPS.pnts <- as.data.frame(all.pnts[!duplicated(all.pnts$LOC_ID),c("LOC_ID","SOURCEID","SOURCEDB","LONWGS84","LATWGS84")])
SPROPS.pnts <- SPROPS.pnts[!is.na(SPROPS.pnts$LONWGS84),]
coordinates(SPROPS.pnts) <- ~ LONWGS84+LATWGS84
proj4string(SPROPS.pnts) <- "+proj=longlat +datum=WGS84"
length(SPROPS.pnts)
#154,408
summary(as.factor(SPROPS.pnts$SOURCEDB))
save(SPROPS.pnts, file="SPROPS.pnts.rda")
save.image()
