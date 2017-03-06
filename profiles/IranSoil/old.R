## Data from (http://www.iransoil.com/webgis.html):
library(gdata)
#perl <- gdata:::findPerl("perl")
perl = "C:/Perl64/bin/perl.exe"
tmp <- read.xls("Iransoil.xlsx", perl=perl, sheet=1)
## takes 2-3 mins to read!
profs2 <- tmp
str(profs2)
#library(sp)
#xy <- profs2[,c("longitude","latitude","classification.soil.taxonomy.")]
#xy <- xy[!is.na(xy$longitude)&!is.na(xy$latitude),]
#coordinates(xy) <- ~longitude+latitude
#proj4string(xy) <- "+proj=longlat +datum=WGS84"
#plotKML(xy)
profs2.f <- rename(profs2, c("classification.soil.taxonomy."="TAXNUSDA", "longitude"="LONWGS84", "latitude"="LATWGS84", "X..sand."="SNDPPT", "X..silt."="SLTPPT", "X..clay."="CLYPPT", "X.pH"="PHIHO5", "Bulk.Density.Ï.b..g.cm3."="BLD", "X..O.C"="ORCDRC"))
profs2.f$UHDICM = 0
profs2.f$LHDICM = 30
profs2.f$ORCDRC <- as.numeric(paste(profs2.f$ORCDRC)) * 10
summary(profs2.f$ORCDRC)
profs2.f$PHIHO5 <- as.numeric(paste(profs2.f$PHIHO5))
summary(profs2.f$PHIHO5)
profs2.f$SOURCEID <- as.factor(paste("IranSoil", profs2.f$id, sep="_"))
profs2.f <- profs2.f[!is.na(profs2.f$LONWGS84)&!is.na(profs2.f$LATWGS84),]

depths(profs2.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs2.f) <- ~ LONWGS84 + LATWGS84 + TAXNUSDA
## check if there are still some duplicates
profs2.f@site$SOURCEID[duplicated(profs2.f@site$SOURCEID)]
## spatial duplicates:
sp2 <- profs2.f@site
coordinates(sp2) <- ~ LONWGS84 + LATWGS84
sp::zerodist(sp2)
## 7 duplicates!

## prepare a SITE and hor table table:
wsp24 <- list(sites=profs2.f@site, horizons=profs2.f@horizons[,c("SNDPPT","SLTPPT","CLYPPT","PHIHO5","BLD","ORCDRC")])
str(wsp24)
wsp24$sites$SOURCEDB = "IranSoil"
lapply(as.list(wsp24$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp24$horizons), function(x){sum(!is.na(x))})
wsp24$sites$TAXNUSDA <- as.character(wsp24$sites$TAXNUSDA)
save(wsp24, file="../wsp24.rda", compress="xz")