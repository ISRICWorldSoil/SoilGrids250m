
##visual check ## no detected issues
View(TAXNWRB.CanSIS$TAXNWRB)
library(tools)
print(showNonASCII(TAXNWRB.CanSIS$TAXNWRB))
View(TAXOUSDA.CanSIS$TAXOUSDA)
print(showNonASCII(as.character(TAXOUSDA.CanSIS$TAXOUSDA)))

fecd.xy <- fecd$dbf[!duplicated(fecd$dbf$FECD_ID),c("FECD_ID","LAT_DD","LONG_DD","SUB_GROUP","ORDER")]
str(fecd.xy)
## 706 profiles
coordinates(fecd.xy) <- ~ LONG_DD+LAT_DD
proj4string(fecd.xy) <- "+proj=longlat +datum=WGS84"
fecd.xy$CSCS <- paste(fecd.xy$SUB_GROUP, fecd.xy$ORDER, sep="_")
plotKML(fecd.xy["CSCS"])


pedon_merge <- join(pedon_id[,c("SOURCEID","PEDON","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB","TAXOUSDA")], pedon_class[,c("PEDON","TAX_ORDER","TAX_GTGRP","TAX_SBGRP")], by="PEDON", type="left")
write.csv(pedon_merge, "pedon_merge.csv")


cGSPD <- odbcConnect(dsn="CANSIS")
sqlTables(cGSPD)$TABLE_NAME
# get horizon table:
pedon_site <- sqlFetch(cGSPD, "pedon_site", stringsAsFactors=FALSE, as.is=TRUE)
str(pedon_site)
# coordinates:
pedon_id <- sqlFetch(cGSPD, "pedon_id", stringsAsFactors=FALSE, as.is=TRUE)
str(pedon_id)



# chemical and physical properties:
pedon_chemical <- sqlFetch(cGSPD, "pedon_chemical", stringsAsFactors=FALSE, as.is=TRUE)
# filter all negative values
for(j in names(pedon_chemical)){
  pedon_chemical[,j] <- ifelse(pedon_chemical[,j]<0, NA, pedon_chemical[,j])
}
str(pedon_chemical)
summary(pedon_chemical$ORGCARB, na.rm=TRUE)
## Organic carbon in permille
pedon_chemical$ORCDRC <- signif(pedon_chemical$ORGCARB * 10/1.724, 3)
summary(pedon_chemical$ORCDRC)
## TH: This is a bit strange - the mean organic carbon is > 5 percent;
pedon_chemical$HORID <- as.factor(paste(pedon_chemical$PEDON, pedon_chemical$LAYERNO, sep="_"))

pedon_physical <- sqlFetch(cGSPD, "pedon_physical", stringsAsFactors=FALSE, as.is=TRUE)
# filter all negative values
for(j in names(pedon_physical)){
  pedon_physical[,j] <- ifelse(pedon_physical[,j]<0, NA, pedon_physical[,j])
}
str(pedon_physical)
pedon_physical$HORID <- as.factor(paste(pedon_physical$PEDON, pedon_physical$LAYERNO, sep="_"))
## morphological properties:
pedon_morphology <- sqlFetch(cGSPD, "pedon_morphology", stringsAsFactors=FALSE, as.is=TRUE)
pedon_morphology$HORID <- as.factor(paste(pedon_morphology$PEDON, pedon_morphology$LAYERNO, sep="_"))

## soil classes:
pedon_class <- sqlFetch(cGSPD, "pedon_class", stringsAsFactors=FALSE, as.is=TRUE)
summary(as.factor(pedon_class$TAX_ORDER))
summary(as.factor(pedon_class$TAX_GTGRP))

# merge tables:
horizon1 <- merge(pedon_chemical, pedon_physical, by="HORID")
horizon <- merge(horizon1, pedon_morphology[,c("HORID","HORIZON","TEXTURE","COLOR1","COLOR2")])
site1 <- merge(pedon_id, pedon_site, by="PEDON")
site <- merge(site1, pedon_class[,c("PEDON","TAX_ORDER","TAX_GTGRP","TAX_SBGRP","FAMILY_PSIZECLAS","FAMILY_DEPTH")], by="PEDON")
# str(site)
# str(horizon)

# attach unique IDs:
site$SOURCEID <- as.factor(paste("CAN", site$PEDON, sep="_"))
horizon$SOURCEID <- as.factor(paste("CAN", horizon$PEDON.x, sep="_"))
## rename some columns:
names(site)[which(names(site)=="YYYYMMDD")] <- "TIMESTRR"
site$TIMESTRR <- as.Date(site$TIMESTRR, format="%Y")
names(site)[which(names(site)=="TAX_ORDER")] <- "TAXOCAN"
names(site)[which(names(site)=="TAX_SBGRP")] <- "TAXGCAN"

names(horizon)[which(names(horizon)=="PS_TSAND")] <- "SNDPPT"
summary(horizon$SLTPPT)
names(horizon)[which(names(horizon)=="PS_50_2U")] <- "SLTPPT"
names(horizon)[which(names(horizon)=="PS_2UCLAY")] <- "CLYPPT"
summary(horizon$CLYPPT)
names(horizon)[which(names(horizon)=="PH1")] <- "PHIHO5"
summary(horizon$PHIHO5)
horizon$PHIHO5 <- ifelse(horizon$PHIHO5<2 | horizon$PHIHO5 > 10, NA, horizon$PHIHO5)
names(horizon)[which(names(horizon)=="PH2")] <- "PHIKCL"
horizon$PHIKCL <- ifelse(horizon$PHIKCL<2 | horizon$PHIKCL > 10, NA, horizon$PHIKCL)
names(horizon)[which(names(horizon)=="UDEPTH.x")] <- "UHDICM"
names(horizon)[which(names(horizon)=="LDEPTH.x")] <- "LHDICM"
names(horizon)[which(names(horizon)=="BD")] <- "BLD"
horizon$BLD <- ifelse(horizon$BLD < .5, NA, horizon$BLD)
summary(horizon$BLD)
horizon$CRFVOL <- horizon$PS_VCSAND # + horizon$PS_NO10 + horizon$PS_NO4
## depth to bedrock (R horizon):
sel.n <- which(horizon$SOURCEID %in% site$SOURCEID[grep("lit", site$FAMILY_DEPTH, ignore.case=TRUE)])
sel.r <- unique(c(sel.n, grep("R", horizon$HORIZON, ignore.case=FALSE, fixed=FALSE)))
horizon$BDRICM <- NA
horizon$BDRICM[sel.r] <- pmax(horizon$UHDICM[sel.r], horizon$LHDICM[sel.r], na.rm=TRUE)
bdr.d <- aggregate(horizon$BDRICM, list(horizon$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
site <- merge(site, bdr.d, all.y=FALSE)
site$BDRICM <- ifelse(site$BDRICM<0, NA, site$BDRICM)

# subset to complete data:
sel.s <- !is.na(site$LONWGS84)&!is.na(site$LATWGS84)
site <- site[sel.s,c("SOURCEID", "LONWGS84", "LATWGS84", "TIMESTRR", "TAXOCAN", "TAXGCAN", "BDRICM")]
sel.h <- !is.na(horizon$UHDICM)&!is.na(horizon$LHDICM)
horizon <- horizon[sel.h,c("SOURCEID", "UHDICM", "LHDICM", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHO5", "PHIKCL", "ORCDRC")]

## convert to SPC class:
profs.f <- join(site, horizon, type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TIMESTRR + TAXOCAN + TAXGCAN + BDRICM
## check if there are still some duplicates
# TH: there are about 30-40 profiles with the same ID but different coordinates
sel <- profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
str(sp::zerodist(sp))
## 208 spatial duplicates...

profs.f@site$SOURCEDB = "CanSIS"
# prepare a SITE and hor table table:
wsp3 <- list(sites=profs.f@site[,c("SOURCEID","LONWGS84","LATWGS84","TIMESTRR","SOURCEDB","BDRICM")], horizons=profs.f@horizons)
str(wsp3)
lapply(as.list(wsp3$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp3$horizons), function(x){sum(!is.na(x))})
save(wsp3, file="../wsp3.rda", compress="xz")

## convert to geosamples:
coordinates(profs.f) <- ~LONWGS84 + LATWGS84
proj4string(profs.f) <- CRS("+proj=longlat +datum=WGS84")
## convert to geosamples:
cansis.geo <- as.geosamples(profs.f)
save(cansis.geo, file="cansis.geo.rda")