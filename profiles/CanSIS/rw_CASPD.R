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
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}

## read exported tables:
pedon_physical <- read.csv("CAN_pedon_physical.csv")
pedon_physical$SAMPLEID <- make.unique(paste("CAN", pedon_physical$PEDON, pedon_physical$LAYERNO, sep="_"))
pedon_chemical <- read.csv("CAN_pedon_chemical.csv")
pedon_chemical$SAMPLEID <- make.unique(paste("CAN", pedon_chemical$PEDON, pedon_chemical$LAYERNO, sep="_"))
#pedon_morphology <- read.csv("CAN_pedon_morphology.csv")
horizons <- plyr::join(pedon_physical, pedon_chemical[,c("SAMPLEID","PH1","PH2","ORGCARB","CEC_BUFF")], by="SAMPLEID")
horizons$SOURCEID <- as.factor(paste("CAN", horizons$PEDON, sep="_"))
## replace "-99.9":
for(j in c("PS_TSAND", "PS_50_2U", "PS_2UCLAY", "PH1", "PH2", "UDEPTH", "LDEPTH", "BD", "PS_VCSAND", "PS_NO10", "PS_NO4", "ORGCARB", "CEC_BUFF")){ 
  horizons[,j] <- ifelse(horizons[,j] < -9, NA, horizons[,j])
}
horizons$CRFVOL <- horizons$PS_VCSAND ## rowSums(horizons[,c("PS_VCSAND", "PS_NO10", "PS_NO4")], na.rm=TRUE)
summary(horizons$CRFVOL)
horizons$ORCDRC <- signif(horizons$ORGCARB * 10/1.724, 3)
horizons <- rename(horizons, c("PS_TSAND"="SNDPPT", "PS_50_2U"="SLTPPT", "PS_2UCLAY"="CLYPPT", "PH1"="PHIHOX", "PH2"="PHIKCL", "UDEPTH"="UHDICM", "LDEPTH"="LHDICM", "BD"="BLD", "CEC_BUFF"="CECSUM"))
horizons$BLD <- horizons$BLD*1000
## many missing lower depths:
horizons$LHDICM <- ifelse(is.na(horizons$LHDICM)&!is.na(horizons$UHDICM), horizons$UHDICM+30, horizons$LHDICM)
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2

## Classification:
CAN_cl <- read.csv("CAN_classes.csv")
CAN_cor <- read.csv("CAN_correlation.csv")
pedon_cor <- join(CAN_cor, CAN_cl, by="CSSC_Great_Groups")
pedon_id <- read.csv("CAN_pedon_id.csv")
pedon_class <- read.csv("CAN_pedon_class.csv")

pedon_id$LONWGS84 <- cols2dms(pedon_id$LOCLONGDEG, pedon_id$LOCLONGMIN, pedon_id$LOCLONGSEC, rep("W", length(pedon_id$LOCLONGDEG)))
pedon_id$LONWGS84 <- ifelse(pedon_id$LONWGS84>-1, NA, pedon_id$LONWGS84) 
pedon_id$LATWGS84 <- cols2dms(pedon_id$LOCLATDEG, pedon_id$LOCLATMIN, pedon_id$LOCLATSEC, rep("N", length(pedon_id$LOCLATDEG)))
pedon_id$LATWGS84 <- ifelse(pedon_id$LATWGS84<25, NA, pedon_id$LATWGS84)
pedon_id$TIMESTRR <- as.Date(paste(pedon_id$YEAR_), format="%Y")
#plot(pedon_id$LONWGS84, pedon_id$LATWGS84)
str(pedon_id)
pedon_id$TAX_GTGRP <- join(pedon_id, pedon_class, by="PEDON", type="left")$TAX_GTGRP
pedon_id[,c("TAXNWRB","TAXOUSDA")] <- join(pedon_id, pedon_cor, by="TAX_GTGRP", type="left", match="first")[,c("WRB_subgroup","USDA_suborder")]
pedon_id$SOURCEDB = "CanSIS"
pedon_id$SOURCEID = paste("CAN", pedon_id$PEDON, sep="_")
pedon_id <- pedon_id[!is.na(pedon_id$LATWGS84)&!is.na(pedon_id$LONWGS84),]

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXOUSDA.CanSIS <- pedon_id[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXOUSDA")]
TAXOUSDA.CanSIS <- TAXOUSDA.CanSIS[!is.na(TAXOUSDA.CanSIS$TAXOUSDA)&!is.na(TAXOUSDA.CanSIS$LONWGS84)&nchar(paste(TAXOUSDA.CanSIS$TAXOUSDA))>0,]
str(TAXOUSDA.CanSIS)
## 3000 profiles
coordinates(TAXOUSDA.CanSIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.CanSIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.CanSIS["TAXOUSDA"])
save(TAXOUSDA.CanSIS, file="TAXOUSDA.CanSIS.rda")

TAXNWRB.CanSIS <- pedon_id[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.CanSIS <- TAXNWRB.CanSIS[!is.na(TAXNWRB.CanSIS$TAXNWRB)&!is.na(TAXNWRB.CanSIS$LONWGS84)&nchar(paste(TAXNWRB.CanSIS$TAXNWRB))>0,]
TAXNWRB.CanSIS$TAXNWRB <- iconv(paste(TAXNWRB.CanSIS$TAXNWRB), "latin1", "UTF-8", "")
str(TAXNWRB.CanSIS)
## 3360 profiles
coordinates(TAXNWRB.CanSIS) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.CanSIS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.CanSIS["TAXNWRB"])
save(TAXNWRB.CanSIS, file="TAXNWRB.CanSIS.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.CanSIS <- plyr::join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","PHIKCL","BLD","ORCDRC","CECSUM")], pedon_id[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], type="left")
SPROPS.CanSIS <- SPROPS.CanSIS[!is.na(SPROPS.CanSIS$LONWGS84) & !is.na(SPROPS.CanSIS$LATWGS84) & !is.na(SPROPS.CanSIS$DEPTH),]
View(SPROPS.CanSIS)
## 17,535 horizons
save(SPROPS.CanSIS, file="SPROPS.CanSIS.rda")

