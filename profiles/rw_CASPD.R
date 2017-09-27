## Reading and writing of the Canadian national soil profile data (3605 profiles);
## Tom.Hengl@isric.org
## The National Soil Profile Database for Canada obtained from the CanSIS;
## For more info contact: Xiaoyuan Geng, Canadian Soil Information Service (CanSIS);


## data publicly available:
#download.file("http://wms1.agr.gc.ca/xgeng/pedon.zip", "pedon.zip")
#unzip("pedon.zip")

library(aqp)
library(sp)
library(plyr)
library(plotKML)
library(foreign)
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}

## read exported tables:
pedon_physical <- read.csv("CAN_pedon_physical.csv")
pedon_physical$SAMPLEID <- make.unique(paste("CAN", pedon_physical$PEDON, pedon_physical$LAYERNO, sep="_"))
pedon_chemical <- read.csv("CAN_pedon_chemical.csv")
pedon_chemical$SAMPLEID <- make.unique(paste("CAN", pedon_chemical$PEDON, pedon_chemical$LAYERNO, sep="_"))
pedon_morphology <- read.csv("CAN_pedon_morphology.csv")
pedon_morphology$SOURCEID <- as.factor(paste("CAN", pedon_morphology$PEDON, sep="_"))

horizons <- plyr::join(pedon_physical, pedon_chemical[,c("SAMPLEID","PH1","PH2","ORGCARB","CEC_BUFF")], by="SAMPLEID")
horizons$SOURCEID <- as.factor(paste("CAN", horizons$PEDON, sep="_"))
## replace "-99.9":
for(j in c("PS_TSAND", "PS_50_2U", "PS_2UCLAY", "PH1", "PH2", "UDEPTH", "LDEPTH", "BD", "PS_VCSAND", "PS_NO10", "PS_NO4", "ORGCARB", "CEC_BUFF")){ 
  horizons[,j] <- ifelse(horizons[,j] < -9, NA, horizons[,j])
}
horizons$CRFVOL <- horizons$PS_VCSAND ## rowSums(horizons[,c("PS_VCSAND", "PS_NO10", "PS_NO4")], na.rm=TRUE)
summary(horizons$CRFVOL)
summary(horizons$ORGCARB) ## These numbers are quite high! Is this really SOC or SOM?
#horizons$ORCDRC <- signif(horizons$ORGCARB * 10/1.724, 3)
horizons$ORCDRC <- round(horizons$ORGCARB * 10)
horizons <- rename(horizons, c("PS_TSAND"="SNDPPT", "PS_50_2U"="SLTPPT", "PS_2UCLAY"="CLYPPT", "PH1"="PHIHOX", "PH2"="PHIKCL", "UDEPTH"="UHDICM", "LDEPTH"="LHDICM", "BD"="BLD", "CEC_BUFF"="CECSUM"))
horizons$BLD <- horizons$BLD*1000
## many missing lower depths:
horizons$LHDICM <- ifelse(is.na(horizons$LHDICM)&!is.na(horizons$UHDICM), horizons$UHDICM+30, horizons$LHDICM)
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2

## hits on visual check:
## 5 UHDICM values < 800cm
str(horizons[which(horizons$UHDICM>800),])
## 8 LHDICM values < 800cm
str(horizons[which(horizons$LHDICM>800),])
## 599 obs with deeper 'upper' horizons
View(horizons[which(horizons$LHDICM-horizons$UHDICM<=0),])
##987 obs with messy texture fractions
TexSum <- rowSums(horizons[,c("CLYPPT","SNDPPT","SLTPPT")], na.rm=TRUE)
hist(TexSum)
selTex <- TexSum > 110 | TexSum < 90
str(horizons[which(TexSum>0&TexSum<90),])
##4 obs above pH 10.5 - rare potential for EXTREME sodic soils (but above 10.5 is questionable calibration?? Or an arid mineral spring ;P) check spatial location for clustering.
str(horizons[which(horizons$PHIHOX>=10.5),])
## as above but 3 obs. Check spatial location and with soil chemist.
str(horizons[which(horizons$PHIKCL>=9.5),])
#104 obs with 0 OC but values for CEC
str(horizons[which(horizons$ORCDRC<=0),])
View(horizons[which(horizons$ORCDRC<=0),])
##2 obs of CEC = 0
str(horizons[which(horizons$CECSUM<=0),])
##23 obs with BLD = 0
str(horizons[which(horizons$BLD==0),])
##9 obs with < 2 pH 
str(horizons[which(horizons$PHIHOX<2),])
## 2 obs > 300 (but close...) - check
str(horizons[which(horizons$CECSUM>=300),])
 
## Classification:
pedon_cor <- read.csv("CAN_classes.csv")
pedon_id <- read.csv("CAN_pedon_id.csv")
pedon_class <- read.csv("CAN_pedon_class.csv")
summary(pedon_class$TAX_ORDER)
summary(pedon_class$TAX_GTGRP)
pedon_class$SOURCEID = paste("CAN", pedon_class$PEDON, sep="_")

pedon_id$LONWGS84 <- cols2dms(pedon_id$LOCLONGDEG, pedon_id$LOCLONGMIN, pedon_id$LOCLONGSEC, rep("W", length(pedon_id$LOCLONGDEG)))
pedon_id$LONWGS84 <- ifelse(pedon_id$LONWGS84>-1, NA, pedon_id$LONWGS84) 
pedon_id$LATWGS84 <- cols2dms(pedon_id$LOCLATDEG, pedon_id$LOCLATMIN, pedon_id$LOCLATSEC, rep("N", length(pedon_id$LOCLATDEG)))
pedon_id$LATWGS84 <- ifelse(pedon_id$LATWGS84<25, NA, pedon_id$LATWGS84)
pedon_id$TIMESTRR <- as.Date(paste(pedon_id$YEAR_), format="%Y")
pedon_id$SOURCEDB = "CanSIS"
#plot(pedon_id$LONWGS84, pedon_id$LATWGS84)
str(pedon_id)
pedon_id <- join(pedon_id, pedon_class, by="PEDON", type="left", match="first")
pedon_id$TAX_ID <- paste(pedon_id$TAX_ORDER, pedon_id$TAX_GTGRP, sep="_")
pedon_merge <- join(pedon_id, pedon_cor[,c("TAX_ID","TAXNWRB","TAXOUSDA")], by="TAX_ID", type="left", match="first")
pedon_merge$SOURCEID = paste("CAN", pedon_merge$PEDON, sep="_")
pedon_merge <- pedon_merge[!is.na(pedon_merge$LATWGS84)&!is.na(pedon_merge$LONWGS84),]
## look at some examples
pedon_merge[1:5,c("PEDON","TAX_ORDER","TAX_GTGRP","TAX_SBGRP","TAXNWRB","TAXOUSDA")]
#pedon_merge.xy <- pedon_merge[,c("LONWGS84","LATWGS84","TAX_ID")]
#coordinates(pedon_merge.xy) <- ~ LONWGS84 + LATWGS84 
#proj4string(pedon_merge.xy) <- "+proj=longlat +datum=WGS84"
#plotKML(pedon_merge.xy["TAX_ID"])

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXOUSDA.CanSIS <- pedon_merge[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXOUSDA")]
TAXOUSDA.CanSIS <- TAXOUSDA.CanSIS[!is.na(TAXOUSDA.CanSIS$TAXOUSDA)&!is.na(TAXOUSDA.CanSIS$LONWGS84)&nchar(paste(TAXOUSDA.CanSIS$TAXOUSDA))>0,]
str(TAXOUSDA.CanSIS)
## 3000 profiles
coordinates(TAXOUSDA.CanSIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.CanSIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.CanSIS["TAXOUSDA"])
## REMOVE SOME NON-SENSICAL POINTS FALLING IN THE OCEAN?
#rem <- which(TAXOUSDA.CanSIS$SOURCEID %in% c("CAN_4408","CAN_4407","CAN_5285","CAN_5288"))
save(TAXOUSDA.CanSIS, file="TAXOUSDA.CanSIS.rda")

TAXNWRB.CanSIS <- pedon_merge[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.CanSIS <- TAXNWRB.CanSIS[!is.na(TAXNWRB.CanSIS$TAXNWRB)&!is.na(TAXNWRB.CanSIS$LONWGS84)&nchar(paste(TAXNWRB.CanSIS$TAXNWRB))>0,]
TAXNWRB.CanSIS$TAXNWRB <- iconv(paste(TAXNWRB.CanSIS$TAXNWRB), "latin1", "UTF-8", "")
str(TAXNWRB.CanSIS)
## 5938 profiles
coordinates(TAXNWRB.CanSIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.CanSIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.CanSIS["TAXNWRB"])
save(TAXNWRB.CanSIS, file="TAXNWRB.CanSIS.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.CanSIS <- plyr::join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","PHIKCL","BLD","ORCDRC","CECSUM")], pedon_id[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], type="left")
SPROPS.CanSIS <- SPROPS.CanSIS[!is.na(SPROPS.CanSIS$LONWGS84) & !is.na(SPROPS.CanSIS$LATWGS84) & !is.na(SPROPS.CanSIS$DEPTH),]
SPROPS.CanSIS$PHIKCL <- ifelse(SPROPS.CanSIS$PHIKCL>11|SPROPS.CanSIS$PHIKCL<2.5, NA, SPROPS.CanSIS$PHIKCL)
SPROPS.CanSIS$PHIHOX <- ifelse(SPROPS.CanSIS$PHIHOX>11.5|SPROPS.CanSIS$PHIHOX<2.5, NA, SPROPS.CanSIS$PHIHOX)
View(SPROPS.CanSIS)
hist(SPROPS.CanSIS$CECSUM)
hist(SPROPS.CanSIS$ORCDRC, breaks=30)
## pool of high values
SPROPS.CanSIS$BLD <- ifelse(SPROPS.CanSIS$BLD<100, NA, SPROPS.CanSIS$BLD)
hist(SPROPS.CanSIS$BLD, breaks=30)
## 17,535 horizons
save(SPROPS.CanSIS, file="SPROPS.CanSIS.rda")
SPROPS.CanSIS.xy = SPROPS.CanSIS
coordinates(SPROPS.CanSIS.xy) <- ~ LONWGS84+LATWGS84
SPROPS.CanSIS.xy$ORCDRC <- round(SPROPS.CanSIS.xy$ORCDRC)
proj4string(SPROPS.CanSIS.xy) <- "+proj=longlat +datum=WGS84"
unlink("SPROPS.CanSIS.shp")
writeOGR(SPROPS.CanSIS.xy, "SPROPS.CanSIS.shp", "SPROPS.CanSIS", "ESRI Shapefile")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

levels(pedon_morphology$HORIZON)
sel.n <- which(horizons$SOURCEID %in% pedon_id$SOURCEID[grep("lit", pedon_class$FAMILY_DEPTH, ignore.case=TRUE)])
sel.r <- which(horizons$SOURCEID %in% pedon_morphology$SOURCEID[grep("^R", pedon_morphology$HORIZON, ignore.case=FALSE, fixed=FALSE)])
sel.r2 <- which(horizons$SOURCEID %in% pedon_morphology$SOURCEID[grep("2R", pedon_morphology$HORIZON, ignore.case=FALSE, fixed=FALSE)])
sel.t <- unique(c(sel.n, sel.r, sel.r2))
horizons$BDRICM <- NA
horizons[sel.t,"BDRICM"] <- pmax(horizons$UHDICM[sel.t], horizons$LHDICM[sel.t])
bdr.d <- aggregate(horizons$BDRICM, list(horizons$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.CanSIS <- join(bdr.d, pedon_id[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
BDR.CanSIS$BDRICM <- ifelse(is.infinite(BDR.CanSIS$BDRICM), 250, BDR.CanSIS$BDRICM)
summary(BDR.CanSIS$BDRICM<250)
## 209 points
str(BDR.CanSIS)
save(BDR.CanSIS, file="BDR.CanSIS.rda")

# ------------------------------------------------------------
# Forest Ecosystem Carbon Database (FECD)
# ------------------------------------------------------------

## http://www.cfs.nrcan.gc.ca/bookstore_pdfs/25626.pdf
spd <- read.csv("SPD_JOIN.csv")
str(spd)
spd$SNDPPT <- 100-(spd$silt+spd$clay)
summary(spd$cec)
spd$BLD <- spd$bulkdens * 1000
spd$ORCDRC <- spd$orgcarb * 10
summary(spd$ORCDRC)
spd$SOURCEID <- paste("FECD", spd$prov, spd$csite, sep="_")
spd$SAMPLEID <- make.unique(paste("FECD", spd$prov, spd$csite, sep="_"))
## longitudes miss "-"??
spd$long <- -spd$long
horizons2a <- rename(spd, replace=c("long"="LONWGS84", "lat"="LATWGS84", "hor"="HZDTXT", "top"="UHDICM", "bottom"="LHDICM", "silt"="SLTPPT", "clay"="CLYPPT", "cec"="CECSUM"))

fecd <- read.delim("FECD_JOIN_ALL.csv")
str(fecd)
summary(fecd$ORG_C_TOT.N.11.0)
hist(fecd$ORG_C_TOT.N.11.0)
summary(fecd$BD.N.9.0) ## many zeros / rounded numbers from conversion to DBF
fecd$SOURCEID <- paste("FECD", fecd$PROFILE_ID.N.12.0, sep="")

horizons2 <- rename(fecd, replace=c("LONG_DD.N.10.2"="LONWGS84", "LAT_DD.N.8.2"="LATWGS84", "PH_H2O.C.9"="PHIHOX", "PH_CACL2.C.9"="PHIKCL", "HORIZON.C.9"="HZDTXT", "SAND.C.9"="SNDPPT", "SILT.C.9"="SLTPPT", "CLAY.C.9"="CLYPPT", "CEC.C.9"="CECSUM", "UPPER.N.9.0"="UHDICM", "LOWER.N.9.0"="LHDICM"))
horizons2$SAMPLEID <- make.unique(paste("FECD", horizons2$PROFILE_ID.N.12.0, horizons2$HORIZON_ID.N.9.0, sep="_"))
summary(horizons2$PHIHOX)
summary(horizons2$SNDPPT)
## Soil organic carbon in g/kg (-CACO3):
horizons2$ORCDRC <- signif(10* ifelse((horizons2$PHIHOX > 7)&!is.na(horizons2$CACO3.C.9), horizons2$ORG_C_TOT.N.11.0 - .12 * horizons2$CACO3.C.9, horizons2$ORG_C_TOT.N.11.0), 3)
summary(horizons2$ORCDRC)
#horizons2$BLD <- ifelse(horizons2$BD<.15, NA, horizons2$BD*1000)
#hist(horizons2$BLD)
summary(horizons2$CF_CLASS.C.9)
CRF.tbl <- data.frame(CF_CLASS.C.9=levels(horizons2$CF_CLASS.C.9), CRFVOL=c(NA,5,70,20,23,48))
horizons2$CRFVOL <- join(horizons2, CRF.tbl, type="left")$CRFVOL
horizons2$TIMESTRR <- as.Date(paste(horizons2$DATE.N.9.0), format="%Y")

xy.FECD <- fecd[,c("SOURCEID","LAT_DD.N.8.2","LONG_DD.N.10.2","TEXTURE.C.9","ESTABLISHM.C.16","PROVINCE.C.9","SLOPE.C.9","ORDER.C.11","PARENT_MAT.C.20","DRAINAGE_C.C.18")]
xy.FECD$LOC_ID <- paste("ID", xy.FECD$LONG_DD.N.10.2, xy.FECD$LAT_DD.N.8.2, sep="_")
xy.FECD <- xy.FECD[!duplicated(xy.FECD$LOC_ID),]
coordinates(xy.FECD) <- ~ LONG_DD.N.10.2+LAT_DD.N.8.2
proj4string(xy.FECD) <- "+proj=longlat +datum=WGS84"
plotKML(xy.FECD, balloon=TRUE, kmz=TRUE)


SPROPS.FECD <- rbind.fill(list(horizons2[,c("SOURCEID","TIMESTRR","LONWGS84","LATWGS84","SAMPLEID","HZDTXT","UHDICM","LHDICM","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","PHIKCL","ORCDRC","CECSUM")], horizons2a[,c("SOURCEID","LONWGS84","LATWGS84","SAMPLEID","HZDTXT","UHDICM","LHDICM","CLYPPT","SNDPPT","SLTPPT","ORCDRC","CECSUM","BLD")]))
SPROPS.FECD$DEPTH <- SPROPS.FECD$UHDICM + (SPROPS.FECD$LHDICM - SPROPS.FECD$UHDICM)/2
SPROPS.FECD <- SPROPS.FECD[!is.na(SPROPS.FECD$LONWGS84) & !is.na(SPROPS.FECD$LATWGS84) & !is.na(SPROPS.FECD$DEPTH),]
SPROPS.FECD$PHIKCL <- ifelse(SPROPS.FECD$PHIKCL>11|SPROPS.FECD$PHIKCL<2.5, NA, SPROPS.FECD$PHIKCL)
SPROPS.FECD$PHIHOX <- ifelse(SPROPS.FECD$PHIHOX>11.5|SPROPS.FECD$PHIHOX<2.5, NA, SPROPS.FECD$PHIHOX)
SPROPS.FECD$BLD <- ifelse(SPROPS.FECD$BLD<100, NA, SPROPS.FECD$BLD)
hist(SPROPS.FECD$CECSUM, breaks=30)
hist(SPROPS.FECD$ORCDRC)
hist(SPROPS.FECD$PHIHOX, breaks=30)
hist(SPROPS.FECD$CLYPPT, breaks=30)
hist(SPROPS.FECD$BLD, breaks=30)
SPROPS.FECD$SOURCEDB = "Can_FECD"
str(SPROPS.FECD)
## 10,761 horizons
save(SPROPS.FECD, file="SPROPS.FECD.rda")
write.csv(SPROPS.FECD, file="SPROPS.FECD.csv")

horizons2a$LOC_ID <- paste("ID", horizons2a$LONWGS84, horizons2a$LATWGS84, sep="_")
sites2a <- horizons2a[!duplicated(horizons2a$LOC_ID),c("SOURCEID","LONWGS84","LATWGS84","taxa_can")]
sites2a$SOURCEDB = "Can_FECD"
summary(sites2a$taxa_can)
str(sites2a)
## merge with correlation table:
pedon_cor$taxa_can <- paste(pedon_cor$TAX_SBGRP, pedon_cor$TAX_GTGRP, sep=".")

TAXNWRB.FECD <- join(sites2a, pedon_cor, type="left", match="first")[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","TAXNWRB")]
TAXNWRB.FECD <- TAXNWRB.FECD[!is.na(TAXNWRB.FECD$TAXNWRB)&!is.na(TAXNWRB.FECD$LONWGS84)&nchar(paste(TAXNWRB.FECD$TAXNWRB))>0,]
str(TAXNWRB.FECD)
## 1214 profiles
coordinates(TAXNWRB.FECD) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.FECD) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.FECD["TAXNWRB"],balloon=TRUE)
save(TAXNWRB.FECD, file="TAXNWRB.FECD.rda")

TAXOUSDA.FECD <- join(sites2a, pedon_cor, type="left", match="first")[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","TAXOUSDA")]
TAXOUSDA.FECD <- TAXOUSDA.FECD[!is.na(TAXOUSDA.FECD$TAXOUSDA)&!is.na(TAXOUSDA.FECD$LONWGS84)&nchar(paste(TAXOUSDA.FECD$TAXOUSDA))>0,]
str(TAXOUSDA.FECD)
## 1214 profiles
coordinates(TAXOUSDA.FECD) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.FECD) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.FECD["TAXOUSDA"])
save(TAXOUSDA.FECD, file="TAXOUSDA.FECD.rda")

# ------------------------------------------------------------
# Canadian Upland Forest Carbon Database
# ------------------------------------------------------------

## read exported tables:
PROFILES_simple <- read.csv("PROFILES_simple.csv")
PROFILES_simple$SOURCEID <- paste("CUFCDB", PROFILES_simple$LOCATION_ID, sep="_")
PROFILES_simple$SAMPLEID <- make.unique(paste(PROFILES_simple$SOURCEID, PROFILES_simple$ID, sep="_"))
## locations:
SITES_simple <- read.csv("SITES_simple.csv")
SITES_simple$SOURCEID <- make.unique(paste("CUFCDB", SITES_simple$LOCATION_ID, sep="_"))
SITES_simple <- rename(SITES_simple, replace=c("LONGITUDE"="LONWGS84", "LATITUDE"="LATWGS84"))
summary(SITES_simple$Loc_accuracy)
sel.loc <- SITES_simple$Loc_accuracy=="high"
SITES_simple$SOURCEDB = "CUFCDB"

## rename columns:
PROFILES_simple$UHDICM <- PROFILES_simple$UPPER_HZN_LIMIT
PROFILES_simple$LHDICM <- PROFILES_simple$UPPER_HZN_LIMIT+PROFILES_simple$HZN_THICKNESS
PROFILES_simple$ORCDRC <- round(PROFILES_simple$ORG_CARBON_PCT*10)
PROFILES_simple$PHIHOX = NA
#PROFILES_simple$PHIKCL = NA
summary(PROFILES_simple$pH_H2O_CACL2)
sel.pH = PROFILES_simple$pH_H2O_CACL2=="H2O"|PROFILES_simple$pH_H2O_CACL2=="unknown"
#sel.CaCl2 = PROFILES_simple$pH_H2O_CACL2=="CaCl2"|PROFILES_simple$pH_H2O_CACL2=="CACL2"
PROFILES_simple$PHIHOX[sel.pH] <- PROFILES_simple$pH[sel.pH]
#PROFILES_simple$PHIKCL[sel.CaCl2] <- PROFILES_simple$pH[sel.CaCl2]
PROFILES_simple$BLD <- round(PROFILES_simple$BULK_DENSITY*1000)  
PROFILES_simple$HZDTXT <- PROFILES_simple$HORIZON
PROFILES_simple$CRFVOL <- round(PROFILES_simple$CF_VOLPCT)
#PROFILES_simple$TIMESTRR <- as.Date(paste(PROFILES_simple$YEAR_), format="%Y")
## dates missing?

SPROPS.CUFCDB <- join(PROFILES_simple[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","BLD","ORCDRC","PHIHOX","CRFVOL","HZDTXT")], SITES_simple[sel.loc,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB")]) ## ,"TIMESTRR"
#SPROPS.CUFCDB$PHIKCL <- ifelse(SPROPS.CUFCDB$PHIKCL>11|SPROPS.CUFCDB$PHIKCL<2.5, NA, SPROPS.CUFCDB$PHIKCL)
SPROPS.CUFCDB$PHIHOX <- ifelse(SPROPS.CUFCDB$PHIHOX>11.5|SPROPS.CUFCDB$PHIHOX<2.5, NA, SPROPS.CUFCDB$PHIHOX)
SPROPS.CUFCDB$BLD <- ifelse(SPROPS.CUFCDB$BLD<100, NA, SPROPS.CUFCDB$BLD)
hist(SPROPS.CUFCDB$ORCDRC)
hist(SPROPS.CUFCDB$PHIHOX, breaks=30)
hist(SPROPS.CUFCDB$BLD, breaks=30)
str(SPROPS.CUFCDB)
## 17,519 horizons
save(SPROPS.CUFCDB, file="SPROPS.CUFCDB.rda")

## merge with correlation table:
SITES_simple$CSSC_Great_Groups = SITES_simple$Great_Group
TAXNWRB.CUFCDB <- join(SITES_simple[sel.loc,], pedon_cor, type="left", match="first")[,c("SOURCEID","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.CUFCDB <- TAXNWRB.CUFCDB[!is.na(TAXNWRB.CUFCDB$TAXNWRB)&!is.na(TAXNWRB.CUFCDB$LONWGS84)&nchar(paste(TAXNWRB.CUFCDB$TAXNWRB))>0,]
str(TAXNWRB.CUFCDB)
## 2046 profiles
coordinates(TAXNWRB.CUFCDB) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.CUFCDB) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.CUFCDB["TAXNWRB"],balloon=TRUE)
save(TAXNWRB.CUFCDB, file="TAXNWRB.CUFCDB.rda")

TAXOUSDA.CUFCDB <- join(SITES_simple[sel.loc,], pedon_cor, type="left", match="first")[,c("SOURCEID","LONWGS84","LATWGS84","TAXOUSDA")]
TAXOUSDA.CUFCDB <- TAXOUSDA.CUFCDB[!is.na(TAXOUSDA.CUFCDB$TAXOUSDA)&!is.na(TAXOUSDA.CUFCDB$LONWGS84)&nchar(paste(TAXOUSDA.CUFCDB$TAXOUSDA))>0,]
str(TAXOUSDA.CUFCDB)
## 2046 profiles
coordinates(TAXOUSDA.CUFCDB) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.CUFCDB) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.CUFCDB["TAXOUSDA"],balloon=TRUE)
save(TAXOUSDA.CUFCDB, file="TAXOUSDA.CUFCDB.rda")

# ------------------------------------------------------------
# NPDB_V2
# The National Pedon database (NPDB) hosted by Agriculture and Agri-food Canada from agriculture regions in Canada
# ------------------------------------------------------------

## locations:
NPDB_loc <- read.csv("NPDB_V2_sum_source_info.csv")
NPDB_loc$SOURCEID <- paste("NPDB_V2", NPDB_loc$PEDON_ID, sep="_")
## physical and chemical props:
NPDB_chem <- read.csv("NPDB_V2_sum_chemical.csv")
NPDB_chem$SOURCEID <- paste("NPDB_V2", NPDB_chem$PEDON_ID, sep="_")
NPDB_chem$SAMPLEID <- make.unique(paste(NPDB_chem$SOURCEID, NPDB_chem$LAYER_ID, sep="_"))
NPDB_psych <- read.csv("NPDB_V2_sum_physical.csv")
NPDB_psych$SOURCEID <- paste("NPDB_V2", NPDB_psych$PEDON_ID, sep="_")
NPDB_psych$SAMPLEID <- make.unique(paste(NPDB_psych$SOURCEID, NPDB_psych$LAYER_ID, sep="_"))
NPDB_raw <- read.csv("NPDB_V2_sum_horizons_raw.csv")
NPDB_raw$SOURCEID <- paste("NPDB_V2", NPDB_raw$PEDON_ID, sep="_")
NPDB_raw$SAMPLEID <- make.unique(paste(NPDB_raw$SOURCEID, NPDB_raw$LAYER_ID, sep="_"))
NPDB_hor = plyr::join_all(list(NPDB_raw, NPDB_chem, NPDB_psych), by = c("SAMPLEID"), type = "full")
str(NPDB_hor)

NPDB_loc <- plyr::rename(NPDB_loc, replace=c("DD_LONG"="LONWGS84", "DD_LAT"="LATWGS84"))
NPDB_loc$SOURCEDB = "NPDB_V2"
NPDB_loc$TIMESTRR = as.Date(paste(NPDB_loc$CAL_YEAR), format="%Y")

## rename columns:
NPDB_hor <- plyr::rename(NPDB_hor, replace=c("U_DEPTH"="UHDICM", "L_DEPTH"="LHDICM", "PH_H2O"="PHIHOX", "PHCACL2"="PHICAL", "CEC"="CECSUM", "T_CLAY"="CLYPPT", "T_SILT"="SLTPPT", "T_SAND"="SNDPPT"))
NPDB_hor$ORCDRC <- round(ifelse(NPDB_hor$CARB_ORG==9, NA, NPDB_hor$CARB_ORG)*10)
NPDB_hor$CRFVOL <- NPDB_hor$C_SAND+NPDB_hor$VC_SAND ## TH: needs to be checked!
summary(NPDB_hor$ORCDRC)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.00    5.00   14.00   98.71   52.00  701.00    9319 
NPDB_hor$BLD = NPDB_hor$BULK_DEN * 1000
summary(NPDB_hor$BLD)
NPDB_hor$HZDTXT = paste(NPDB_hor$HZN_MAS, NPDB_hor$HZN_SUF, sep="_")

SPROPS.NPDB_V2 <- plyr::join(NPDB_hor[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","BLD","ORCDRC","PHIHOX","CRFVOL","HZDTXT")], NPDB_loc[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","TIMESTRR")]) ## 
SPROPS.NPDB_V2$PHIHOX <- ifelse(SPROPS.NPDB_V2$PHIHOX>11.5|SPROPS.NPDB_V2$PHIHOX<2.5, NA, SPROPS.NPDB_V2$PHIHOX)
SPROPS.NPDB_V2$BLD <- ifelse(SPROPS.NPDB_V2$BLD<50 | SPROPS.NPDB_V2$BLD>2400, NA, SPROPS.NPDB_V2$BLD)
hist(SPROPS.NPDB_V2$ORCDRC)
hist(SPROPS.NPDB_V2$PHIHOX, breaks=30)
hist(SPROPS.NPDB_V2$BLD, breaks=30)
str(SPROPS.NPDB_V2)
## 21,416 horizons
save(SPROPS.NPDB_V2, file="SPROPS.NPDB_V2.rda")

# ------------------------------------------------------------
# SPD from the Canadian Forest Service (CFS) for the non-agricultural areas in Canada
# ------------------------------------------------------------

## already joined:
fecd <- read.csv("fecd_join_all.csv")
fecd$SOURCEID <- paste("FECD", fecd$fecd_id, sep="_")
fecd <- plyr::rename(fecd, replace=c("long"="LONWGS84", "lat"="LATWGS84", "upper"="UHDICM", "lower"="LHDICM", "ph_h2o"="PHIHOX", "ph_cacl2"="PHICAL", "cec"="CECSUM", "clay"="CLYPPT", "silt"="SLTPPT", "sand"="SNDPPT", "horizon"="HZDTXT"))
levels(fecd[,"cf_class"])
fecd$CRFVOL <- join(fecd["cf_class"], data.frame(CRFVOL=c(3,70,20,23,46), cf_class=c("<10%",">65%","10-30%","23","31-65%")))$CRFVOL ## TH: needs to be checked!
fecd$ORCDRC <- round(ifelse(fecd$org_c_tot<0, NA, fecd$org_c_tot)*10)
summary(fecd$ORCDRC)
fecd$BLD = ifelse(fecd$bd < 50/1000 | fecd$bd > 2400/1000, NA, fecd$bd) * 1000
summary(fecd$BLD)
fecd$SOURCEDB = "FECD"
fecd$SAMPLEID = paste("FECD", fecd$profile_id, fecd$horizon_id, sep="_")
fecd$TIMESTRR = as.Date(paste(fecd$date), format="%Y")
SPROPS.FECD <- fecd[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","BLD","ORCDRC","PHIHOX","CRFVOL","HZDTXT","LONWGS84","LATWGS84","SOURCEDB","TIMESTRR")] 
str(SPROPS.FECD)
## 3547
#save(SPROPS.FECD, file="SPROPS.FECD.rda")


# ------------------------------------------------------------
# Depth to bedrock CUFCDB
# ------------------------------------------------------------

levels(PROFILES_simple$HORIZON)
sel2.r <- grep("^R", PROFILES_simple$HORIZON, ignore.case=FALSE, fixed=FALSE)
#sel2.r2 <- which(SITES_simple$SOURCEID %in% PROFILES_simple$SOURCEID[grep("2R", PROFILES_simple$HORIZON, ignore.case=FALSE, fixed=FALSE)])
PROFILES_simple$BDRICM <- NA
PROFILES_simple[sel2.r,"BDRICM"] <- pmax(PROFILES_simple$UHDICM[sel2.r], PROFILES_simple$LHDICM[sel2.r])
bdr2.d <- aggregate(PROFILES_simple$BDRICM, list(PROFILES_simple$SOURCEID), max, na.rm=TRUE)
names(bdr2.d) <- c("SOURCEID", "BDRICM")
BDR.CUFCDB <- join(bdr2.d, SITES_simple[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
BDR.CUFCDB$BDRICM <- ifelse(is.infinite(BDR.CUFCDB$BDRICM), 250, BDR.CUFCDB$BDRICM)
summary(BDR.CUFCDB$BDRICM<250)
## 303 points
str(BDR.CUFCDB)
save(BDR.CUFCDB, file="BDR.CUFCDB.rda")
