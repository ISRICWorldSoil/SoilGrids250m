# title         : rw_ISIS.R
# purpose       : Reading and writing of ISRIC ISIS (785 profiles);
# reference     : ISRIC ISIS soil profile DB [http://www.isric.org/projects/isric-soil-information-system-isis]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL
# inputs        : "ISISExpFeb2012.mdb"
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : horrible DB!

library(RODBC)
library(aqp)
library(plyr)
library(sp)
library(rgdal)
library(reshape)

cGSPD <- odbcConnect(dsn="ISIS")
# odbcGetInfo(cGSPD)
sqlTables(cGSPD)$TABLE_NAME  # 245 tables
ClassificationResults <- sqlFetch(cGSPD, "ClassificationResults", as.is=TRUE)  ## classifications, we need ValueID=209  WRB Soil Group
ClassificationSamples <- sqlFetch(cGSPD, "ClassificationSamples", as.is=TRUE)
SitedescriptionResults <- sqlFetch(cGSPD, "SitedescriptionResults", as.is=TRUE)
str(SitedescriptionResults)
SitedescriptionSamples <- sqlFetch(cGSPD, "SitedescriptionSamples", as.is=TRUE)
AnalyticalSamples <- sqlFetch(cGSPD, "AnalyticalSamples", as.is=TRUE)
AnalyticalResults <- sqlFetch(cGSPD, "AnalyticalResults", as.is=TRUE)
str(AnalyticalResults)
Valuedescriptors <- sqlFetch(cGSPD, "Valuedescriptors", as.is=TRUE)
str(Valuedescriptors)
Sites <- sqlFetch(cGSPD, "Sites", as.is=TRUE)

## create sites table:
sites <- data.frame(SiteId=Sites$Id, SOURCEID=paste(Sites$CountryISO, Sites$SiteNumber, sep=""))
tax.wrb <- subset(ClassificationResults, ValueId==209)[,c("SampleId","Value")] # 262 completed only;
xx <- merge(ClassificationSamples[,c("Id","SiteId")], tax.wrb, by.x="Id", by.y="SampleId")
sites$TAXGWRB <- merge(sites[,c("SiteId","SOURCEID")], xx[,c("SiteId","Value")], all.x=TRUE)$Value
tax.fao1988 <- subset(ClassificationResults, ValueId==185)[,c("SampleId","Value")] # 575 points
xx2 <- merge(ClassificationSamples[,c("Id","SiteId")], tax.fao1988, by.x="Id", by.y="SampleId")
sites$TAXNWRB <- merge(sites[,c("SiteId","SOURCEID")], xx2[,c("SiteId","Value")], all.x=TRUE)$Value
tax.kst <- subset(ClassificationResults, ValueId==195)[,c("SampleId","Value")] # 723 points
xx3 <- merge(ClassificationSamples[,c("Id","SiteId")], tax.kst, by.x="Id", by.y="SampleId")
sites$TAXNUSDA <- merge(sites[,c("SiteId","SOURCEID")], xx3[,c("SiteId","Value")], all.x=TRUE)$Value
lat.tbl <- subset(SitedescriptionResults, ValueId==235)[,c("SampleId","Value")]
xx4 <- merge(SitedescriptionSamples[,c("Id","SiteId")], lat.tbl, by.x="Id", by.y="SampleId")
sites$LATWGS84 <- as.numeric(merge(sites[,c("SiteId","SOURCEID")], xx4[,c("SiteId","Value")], all.x=TRUE)$Value)
lon.tbl <- subset(SitedescriptionResults, ValueId==236)[,c("SampleId","Value")]
xx5 <- merge(SitedescriptionSamples[,c("Id","SiteId")], lon.tbl, by.x="Id", by.y="SampleId")
sites$LONWGS84 <- as.numeric(merge(sites[,c("SiteId","SOURCEID")], xx5[,c("SiteId","Value")], all.x=TRUE)$Value)
year.tbl <- subset(SitedescriptionResults, ValueId==224)[,c("SampleId","Value")]
xx6 <- merge(SitedescriptionSamples[,c("Id","SiteId")], year.tbl, by.x="Id", by.y="SampleId")
sites$TIMESTRR <- merge(sites[,c("SiteId","SOURCEID")], xx6[,c("SiteId","Value")], all.x=TRUE)$Value
sites$TIMESTRR <- as.Date(ifelse(sites$TIMESTRR==0, NA, sites$TIMESTRR), format="%Y") 
depth1.tbl <- subset(SitedescriptionResults, ValueId==249)[,c("SampleId","Value")]
xx7 <- merge(SitedescriptionSamples[,c("Id","SiteId")], depth1.tbl, by.x="Id", by.y="SampleId")
sites$BDRICM <- as.numeric(merge(sites[,c("SiteId","SOURCEID")], xx7[,c("SiteId","Value")], all.x=TRUE)$Value)
summary(sites$BDRICM)
sites$BDRICM <- ifelse(is.na(sites$BDRICM)|sites$BDRICM>250, 250, sites$BDRICM)
sites$SOURCEDB = "ISIS"

View(sites)
length(levels(sites$SOURCEID))
plyr:::nunique(sites$SOURCEID)
sel.c <- !is.na(sites$LONWGS84)&!is.na(sites$LATWGS84)
sites.f <- sites[sel.c,]
sites.f[sites.f$SOURCEID=="CI2","LATWGS84"] = 5.883333

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXOUSDA.ISIS <- sites.f[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNUSDA")]
## 785 profiles
TAXOUSDA.ISIS <- TAXOUSDA.ISIS[!is.na(TAXOUSDA.ISIS$TAXNUSDA)&nchar(paste(TAXOUSDA.ISIS$TAXNUSDA))>0,]
coordinates(TAXOUSDA.ISIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.ISIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.ISIS["TAXNUSDA"])
save(TAXOUSDA.ISIS, file="TAXOUSDA.ISIS.rda")

TAXNWRB.ISIS <- sites.f[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
## 785 profiles
TAXNWRB.ISIS <- TAXNWRB.ISIS[!is.na(TAXNWRB.ISIS$TAXNWRB)&nchar(paste(TAXNWRB.ISIS$TAXNWRB))>0,]
coordinates(TAXNWRB.ISIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.ISIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.ISIS["TAXNWRB"])
save(TAXNWRB.ISIS, file="TAXNWRB.ISIS.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

# prepare horizons table:
hs <- data.frame(SampleId=AnalyticalSamples$Id, UHDICM=AnalyticalSamples$Top, LHDICM=AnalyticalSamples$Bottom, SiteId=AnalyticalSamples$SiteId)
hs$LHDICM <- as.numeric(gsub(">", "", hs$LHDICM))
hs <- join(hs, sites.f[,c("SiteId","SOURCEDB","SOURCEID","LONWGS84","LATWGS84")], type="left")
pHHO5.tbl <- subset(AnalyticalResults, ValueId==1)[,c("SampleId","Value")]
names(pHHO5.tbl)[2] <- "PHIHOX"
pHKCL.tbl <- subset(AnalyticalResults, ValueId==2)[,c("SampleId","Value")]
names(pHKCL.tbl)[2] <- "PHIKCL"
CRF.tbl <- subset(AnalyticalResults, ValueId==22)[,c("SampleId","Value")]
names(CRF.tbl)[2] <- "CRFVOL"
ORC.tbl <- subset(AnalyticalResults, ValueId==4)[,c("SampleId","Value")]
names(ORC.tbl)[2] <- "ORCDRC"
SND.tbl <- subset(AnalyticalResults, ValueId==28)[,c("SampleId","Value")]
names(SND.tbl)[2] <- "SNDPPT"
SLT.tbl <- subset(AnalyticalResults, ValueId==31)[,c("SampleId","Value")]
names(SLT.tbl)[2] <- "SLTPPT"
CLY.tbl <- subset(AnalyticalResults, ValueId==32)[,c("SampleId","Value")]
names(CLY.tbl)[2] <- "CLYPPT"
CEC.tbl <- subset(AnalyticalResults, ValueId==14)[,c("SampleId","Value")]
names(CEC.tbl)[2] <- "CECSUM"
BLD.tbl <- subset(AnalyticalResults, ValueId==34)[,c("SampleId","Value")]
names(BLD.tbl)[2] <- "BLD"
Dep.tbl <- 
horizons <- plyr::join_all(list(hs, pHHO5.tbl, pHKCL.tbl, CRF.tbl, ORC.tbl, SND.tbl, SLT.tbl, CLY.tbl, CEC.tbl, BLD.tbl))
str(horizons)
## 6497
for(j in c("PHIKCL","PHIHOX","ORCDRC","CRFVOL","SNDPPT","SLTPPT","CLYPPT","CECSUM","BLD")){
  horizons[,j] <- as.numeric(horizons[,j])
}
horizons$ORCDRC <- horizons$ORCDRC*10    ## OC in permilles
summary(horizons$ORCDRC)
summary(horizons$PHIHOX)
horizons$BLD <- horizons$BLD*1000
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
horizons$SAMPLEID <- make.unique(paste(horizons$SampleId))

SPROPS.ISIS <- horizons[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIHOX","ORCDRC","BLD","CECSUM","CRFVOL")]
SPROPS.ISIS <- SPROPS.ISIS[!is.na(SPROPS.ISIS$LONWGS84) & !is.na(SPROPS.ISIS$LATWGS84) & !is.na(SPROPS.ISIS$DEPTH),]
str(SPROPS.ISIS)
## 5616
save(SPROPS.ISIS, file="SPROPS.ISIS.rda")
plot(SPROPS.ISIS$LONWGS84, SPROPS.ISIS$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

BDR.ISIS <- sites[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","BDRICM")]
str(BDR.ISIS)
summary(BDR.ISIS$BDRICM<250)
## 131 point
save(BDR.ISIS, file="BDR.ISIS.rda")

