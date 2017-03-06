
PHIHO5 <- merge(AnalyticalSamples[,c("Id","SiteId")], pHHO5.tbl, by.x="Id", by.y="SampleId")
names(PHIHO5)[3] <- "PHIHO5"
horizons <- merge(horizons, PHIHO5[,c("Id","PHIHO5")], all.x=TRUE, all.y=FALSE, sort=FALSE)
pHKCL.tbl <- subset(AnalyticalResults, ValueId==2)[,c("SampleId","Value")]
names(pHKCL.tbl)[2] <- "PHIKCL"
PHIKCL <- merge(AnalyticalSamples[,c("Id","SiteId")], pHKCL.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, PHIKCL[,c("Id","PHIKCL")], all.x=TRUE, all.y=FALSE, sort=FALSE)
CRF.tbl <- subset(AnalyticalResults, ValueId==22)[,c("SampleId","Value")]
names(CRF.tbl)[2] <- "CRFVOL"
CRFVOL <- merge(AnalyticalSamples[,c("Id","SiteId")], CRF.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, CRFVOL[,c("Id","CRFVOL")], all.x=TRUE, all.y=FALSE, sort=FALSE)
ORC.tbl <- subset(AnalyticalResults, ValueId==4)[,c("SampleId","Value")]
names(ORC.tbl)[2] <- "ORCDRC"
ORCDRC <- merge(AnalyticalSamples[,c("Id","SiteId")], ORC.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, ORCDRC[,c("Id","ORCDRC")], all.x=TRUE, all.y=FALSE, sort=FALSE)
SND.tbl <- subset(AnalyticalResults, ValueId==28)[,c("SampleId","Value")]
names(SND.tbl)[2] <- "SNDPPT"
SNDPPT <- merge(AnalyticalSamples[,c("Id","SiteId")], SND.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, SNDPPT[,c("Id","SNDPPT")], all.x=TRUE, all.y=FALSE, sort=FALSE)
SLT.tbl <- subset(AnalyticalResults, ValueId==31)[,c("SampleId","Value")]
names(SLT.tbl)[2] <- "SLTPPT"
SLTPPT <- merge(AnalyticalSamples[,c("Id","SiteId")], SLT.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, SLTPPT[,c("Id","SLTPPT")], all.x=TRUE, all.y=FALSE, sort=FALSE)
CLY.tbl <- subset(AnalyticalResults, ValueId==32)[,c("SampleId","Value")]
names(CLY.tbl)[2] <- "CLYPPT"
CLYPPT <- merge(AnalyticalSamples[,c("Id","SiteId")], CLY.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, CLYPPT[,c("Id","CLYPPT")], all.x=TRUE, all.y=FALSE, sort=FALSE)
CEC.tbl <- subset(AnalyticalResults, ValueId==14)[,c("SampleId","Value")]
names(CEC.tbl)[2] <- "CEC"
CEC <- merge(AnalyticalSamples[,c("Id","SiteId")], CEC.tbl, by.x="Id", by.y="SampleId")
horizons <- merge(horizons, CEC[,c("Id","CEC")], all.x=TRUE, all.y=FALSE, sort=FALSE)
BLD.tbl <- subset(AnalyticalResults, ValueId==34)[,c("SampleId","Value")]
names(BLD.tbl) <- c("Id","BLD")
#BLD <- merge(AnalyticalSamples[,c("Id","SiteId")], BLD.tbl, by.x="Id", by.y="SampleId")  ## problem with double Id's for BLD!!
## reformat columns:
horizons$PHIHOX <- as.numeric(horizons$PHIHO5)
horizons$PHIKCL <- as.numeric(horizons$PHIKCL)
horizons$CRFVOL <- as.numeric(horizons$CRFVOL)
horizons$ORCDRC <- as.numeric(horizons$ORCDRC)*10
horizons$SNDPPT <- as.numeric(horizons$SNDPPT)
horizons$SLTPPT <- as.numeric(horizons$SLTPPT)
horizons$CLYPPT <- as.numeric(horizons$CLYPPT)
horizons$CEC <- as.numeric(horizons$CEC)
## clean up:
sel.c <- complete.cases(horizons[,c("UHDICM","LHDICM")]) ## ,"PHIKCL","CRFVOL","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC"
horizons <- join(horizons[sel.c,], BLD.tbl, type="left")
horizons$BLD <- as.numeric(horizons$BLD)


## convert to SPC class:
profs <- join(sites.f[,c("SOURCEID","TAXNWRB","TAXGWRB","TAXNUSDA","LATWGS84","LONWGS84","TIMESTRR","BDRICM")], horizons[sel.c,-which(names(horizons)=="SiteId")])
depths(profs) <- SOURCEID ~ UHDICM + LHDICM
site(profs) <- ~ LONWGS84 + LATWGS84 + TIMESTRR + TAXNWRB + TAXGWRB + TAXNUSDA + BDRICM
#coordinates(profs) <- ~ LONWGS84 + LATWGS84
sp <- profs@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
str(sel)
## 86 spatial duplicates!
plyr:::nunique(profs@site$SOURCEID)
plyr:::nunique(profs@horizons$SOURCEID)
#str(profs)

## export:
wsp12 <- list(sites=profs@site, horizons=profs@horizons[,c("SOURCEID","UHDICM","LHDICM","CRFVOL","PHIHOX","PHIKCL","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC","BLD")])
str(wsp12)
wsp12$sites$TAXGWRB <- as.character(wsp12$sites$TAXGWRB)
wsp12$sites$TAXNWRB <- NULL
lapply(as.list(wsp12$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp12$horizons), function(x){sum(!is.na(x))})
save(wsp12, file="../wsp12.rda", compress="xz")
write.csv(wsp12$horizons, file="ISIS_horizons.csv")
write.csv(wsp12$sites, file="ISIS_sites.csv")

isis <- wsp12
save(isis, file="D:/Rdev/GSIF/pkg/data/isis.rda", compress="xz")

library(GSIF)
library(plotKML)
profs.geo <- as.geosamples(profs)
profs.PH <- subset(profs.geo, "PHIHO5")
profs.PH$observedValue <- as.numeric(profs.PH$observedValue)
summary(profs.PH$observedValue)
coordinates(profs.PH) <- ~ longitude + latitude
proj4string(profs.PH) <- "+proj=longlat +datum=WGS84"
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
data(R_pal)
kml(profs.PH, shape=shape, colour=observedValue, file="ISIS_PHIHO5.kml", altitude=200+profs.PH$altitude*100, extrude=TRUE, labels=observedValue, colour_scale=R_pal[["pH_pal"]], altitudeMode="relativeToGround", kmz=TRUE)
kml(profs, file="ISIS_profiles.kml", var.name="PHIHO5", balloon = TRUE, kmz=TRUE)

library(plotGoogleMaps)
pal <- R_pal[["pH_pal"]]
mp <- plotGoogleMaps(profs.PH, zcol='observedValue', colPalette=pal)
## end of the script;

isis.pnt <- readOGR("Profiles_real.shp", "Profiles_real")
str(isis.pnt@data)
names(isis.pnt@data)[1] <- "SampleId"
isis.tbl <- data.frame(isis.pnt)
isis.fao1988.tbl <- merge(isis.tbl, tax.fao1988, by="SampleId")
coordinates(isis.fao1988.tbl) <-~coords.x1+coords.x2
proj4string(isis.fao1988.tbl) <- CRS("+proj=longlat +ellps=WGS84")
isis.fao1988.tbl$Value <- as.factor(isis.fao1988.tbl$Value)
# spplot(isis.fao1988.tbl, "Value", col.regions=bpy.colors(length(levels(isis.fao1988.tbl$Value))))


