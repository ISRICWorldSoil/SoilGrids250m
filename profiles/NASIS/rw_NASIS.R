# title         : rw_NASIS.R
# purpose       : Reading and writing of NASIS pedon data;
# reference     : NASIS [http://soils.usda.gov/technical/nasis/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : "nasis_desc.mdb" MS Access dbs
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : LARGE dataset!

## Download the database:
# download.file("http://globalsoilmap.net/data/nasis_desc.7z", destfile=paste(getwd(),"nasis_desc.7z",sep="/"))

library(RODBC)
library(aqp)
library(sp)
library(plyr)
library(rgdal)
load("nasis_desc.RData")
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}

# ------------------------------------------------------------
# Fetch tables - NASIS
# ------------------------------------------------------------

cNASIS <- odbcConnect(dsn="NASIS")
sqlTables(cNASIS)$TABLE_NAME
## get tables:
site <- sqlFetch(cNASIS, "nasis_site", stringsAsFactors=FALSE)
str(site)  # 23,121 profiles!
## ID column: "siteiid"
pedon <- sqlFetch(cNASIS, "nasis_pedon", stringsAsFactors=FALSE)
str(pedon) # "siteiidref"
siteobs <- sqlFetch(cNASIS, "nasis_siteobs", stringsAsFactors=FALSE)
str(siteobs) # "siteiidref"
## horizon data:
phorizon <- sqlFetch(cNASIS, "nasis_phorizon", stringsAsFactors=FALSE)
str(phorizon) # "peiidref"
summary(phorizon$ec)

# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

# add missing columns:
phorizon$DEPTH <- phorizon$hzdept + (phorizon$hzdepb - phorizon$hzdept)/2
# mask out profiles with missing coordinates:
site <- site[!is.na(site$longdegrees)&!is.na(site$latdegrees),]
site$longitude_it <- ifelse(site$longdir=="west", "W", "E")
site$latitude_it <- ifelse(site$latdir=="north"|site$latdir=="North", "N", "S")
site$LAT <- cols2dms(site$latdegrees, site$latminutes, site$latseconds, site$latitude_it)
site$LON <- cols2dms(site$longdegrees, site$longminutes, site$longseconds, site$longitude_it)
summary(site$LAT) # 563 NA's
summary(site$LON)
# no use of profiles that have no coordinates:
site <- site[!is.na(site$LAT)&!is.na(site$LON),]
summary(as.factor(site$horizdatnm))
sp.ll <- site[,c("LAT","LON")]
coordinates(sp.ll) <- ~LON+LAT
proj4string(sp.ll) <- CRS("+proj=longlat +datum=NAD83")
ll <- data.frame(spTransform(sp.ll, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
site$LONWGS84 <- ll[,1]
site$LATWGS84 <- ll[,2]


# ------------------------------------------------------------
# create the horizon and site tables (takes time!)
# ------------------------------------------------------------

# sites table:
s0 <- merge(siteobs, pedon, by.y="siteiidref", by.x="siteiidref", all=TRUE)
str(s0)
s1 <- merge(site, s0, by.x="siteiid", by.y="siteiidref", all.x=TRUE)
## NOTE: the connector does not have the same name;
# plot(s1$LON, s1$LAT, pch="+")
s1$siteiid <- as.factor(s1$siteiid)
s1 <- s1[!duplicated(s1$siteiid),]  
# 18,034 points! about 3000 points disappeared!
## fix some columns ??

# horizons table:
hs <- merge(phorizon, s1[,c("siteiid","peiid")], by.x="peiidref", by.y="peiid", all.y=TRUE)
summary(as.factor(hs$seqnum))
# take out profiles with > 14 horizons
hs <- hs[!is.na(hs$DEPTH)&hs$seqnum<15&!is.na(hs$seqnum),]
# remove duplicate horizon keys:
hs <- hs[!duplicated(hs$phiid),]
length(levels(as.factor(hs$phiid)))
# 117218
summary(as.factor(hs$texture))
hs <- merge(hs, s1[,c("siteiid", "usiteid")], by="siteiid")


# ------------------------------------------------------------
# SoilProfileCollection
# ------------------------------------------------------------

sel1 <- c("siteiid","obsdate","usiteid","LONWGS84","LATWGS84","pedonpurpose","pedontype","geomposhill","hillslopeprof","shapeacross","shapedown","drainagecl","siteperm","runoff","soinmassamp","taxorder","taxsuborder","taxgrtgroup","taxsubgrp","taxpartsize","taxtempregime","bedrckdepth")
sel2 <- c("siteiid","usiteid","phiid","hzname","seqnum","hzdept","hzdepb","texture","rupresblkmst","phfield","phdetermeth","bounddistinct","boundtopo")
NASIS_all <- list(sites=s1[,sel1], horizons=hs[,sel2])
str(NASIS_all)

## horizons table / rename columns:
sites <- rename(NASIS_all$sites, c("siteiid"="SOURCEID", "obsdate"="TIMESTRR", "taxsubgrp"="TAXNUSDA", bedrckdepth="BDRICM"))
sites$SOURCEDB <- "NASIS" 
horizons <- rename(NASIS_all$horizons, c("siteiid"="SOURCEID","hzdept"="UHDICM", "hzdepb"="LHDICM", "phfield"="PHIHO5", "texture"="TEXMHT"))
sites$TIMESTRR <- as.Date(sites$TIMESTRR)

profs.f <- join(sites[,c("SOURCEID","LONWGS84","LATWGS84","TIMESTRR","TAXNUSDA","SOURCEDB", "BDRICM")], horizons[,c("SOURCEID","UHDICM","LHDICM","PHIHO5","TEXMHT")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
## THIS TAKES 5-10 mins and about 4GB RAM!!!
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TAXNUSDA + TIMESTRR + SOURCEDB + BDRICM
## check if there are still some duplicates
summary(duplicated(profs.f@site$SOURCEID))
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
str(sp::zerodist(sp))
## 1966 duplicate points!
# prepare a SITE and hor table table:
wsp1 <- list(sites=profs.f@site, horizons=profs.f@horizons)
str(wsp1)
save(wsp1, file="../wsp1.rda", compress="xz")
lapply(as.list(wsp1$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp1$horizons), function(x){sum(!is.na(x))})

# subset to Indiana state:
sel <- s1$LAT>36.5&s1$LAT<42.5&s1$LON< -84.5&s1$LON> -87.5
s1_in <- s1[sel,sel1]
hs_in <- hs[hs$siteiid %in% s1_in$siteiid,sel2]
NASIS <- list(sites=s1_in, horizons=hs_in)
str(NASIS)
plot(s1_in$LON, s1_in$LAT, pch="+")
# Save object:
save(NASIS, file="NASIS.rda", compress="xz")
save.image("nasis_desc.RData")

# end of script;
