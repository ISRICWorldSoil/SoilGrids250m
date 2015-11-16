# title         : rw_Russia.R
# purpose       : Reading and writing of Russian profiles (234 profiles);
# reference     : Russian SOIL REFERENCE PROFILES from LAND RESOURCES OF RUSSIA CD-ROM [http://webarchive.iiasa.ac.at/Research/FOR/russia_cd/guide.htm]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Sept 2013.
# inputs        : DBF files with some serious formatting problems
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : For more info see: Stolbovoi, V., McCallum, I., 2002. Land Resources of Russia. CD-ROM.

library(aqp)
library(sp)
library(rgdal)
#library(foreign)
library(plyr)

profs <- read.csv("profs.csv", sep=";") 
str(profs) 
## copy rows:
for(i in 2:nrow(profs)){
  if(is.na(profs[i,"LAT"])){ profs[i,"LAT"] <- profs[i-1,"LAT"] }
  if(is.na(profs[i,"LONG"])){ profs[i,"LONG"] <- profs[i-1,"LONG"] }
  if(is.na(profs[i,"SOURCE"])){ profs[i,"SOURCE"] <- profs[i-1,"SOURCE"] }
  if(is.na(profs[i,"PM"])|profs[i,"PM"]==""){ profs[i,"PM"] <- profs[i-1,"PM"] }
}
x <- strsplit(paste(profs$DEPTH), "-")
profs$UHDICM <- as.numeric(sapply(x, function(y){y[1]}))
profs$LHDICM <- as.numeric(sapply(x, function(y){y[2]}))
profs$SOURCEID <- paste("RUS", profs$RUS, profs$LONG, profs$LAT, sep="_")
profs.f <- rename(profs, c("LAT"="LATWGS84", "LONG"="LONWGS84", "CLAY"="CLYPPT", "SILT"="SLTPPT", "CARB"="ORCDRC", "PH.H2O."="PHIHO5", "BULK_DENSI"="BLD", "D_ROCK"="BDRICM"))
profs.f$SNDPPT <- profs.f$SAND1+profs.f$SAND2
#str(profs.f)
sel <- c(which(rowSums(profs.f[,c("SLTPPT","CLYPPT","SNDPPT")]) < 93), which(rowSums(profs.f[,c("SLTPPT","CLYPPT","SNDPPT")]) > 107))
profs.f[sel,c("SLTPPT","CLYPPT","SNDPPT")] <- NA
summary(rowSums(profs.f[,c("SLTPPT","CLYPPT","SNDPPT")]))
## 490 horizons were removed!
summary(profs.f$ORCDRC)
profs.f$ORCDRC <- profs.f$ORCDRC * 10
summary(profs.f$BLD)
plyr:::nunique(profs.f$SOURCEID)
## find non-unique BDRICM values:

## convert to SPC class:
profs.spc <- profs.f[,c("SOURCEID","LATWGS84","LONWGS84","PM","UHDICM","LHDICM","PHIHO5","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC","BLD","BDRICM")]
#View(profs.spc)
depths(profs.spc) <- SOURCEID ~ UHDICM + LHDICM
site(profs.spc) <- ~ LONWGS84 + LATWGS84 # + PM
profs.spc@site$SOURCEDB = "RusSPDB"
sp <- profs.spc@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
## 4 spatial duplicates!

## Copy information about the bedrock:
bdr.d <- aggregate(profs.spc@horizons$BDRICM, list(profs.spc@horizons$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
m <- merge(profs.spc@site, bdr.d, all.y=FALSE)
m$BDRICM <- ifelse(m$BDRICM<2, NA, m$BDRICM)

## export:
wsp13 <- list(sites=m, horizons=profs.spc@horizons[,c("SOURCEID","UHDICM","LHDICM","PHIHO5","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC","BLD")])
str(wsp13)
lapply(as.list(wsp13$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp13$horizons), function(x){sum(!is.na(x))})
save(wsp13, file="../wsp13.rda", compress="xz")


## end of script;