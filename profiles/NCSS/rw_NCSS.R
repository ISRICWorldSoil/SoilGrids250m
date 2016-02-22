## SoilGrids.org project
## NCSS Characterization Database [http://ncsslabdatamart.sc.egov.usda.gov/];
## "NCSS_Soil_Characterization_Database_3_10_2015.mdb"
## by: Tom.Hengl@isric.org

#library(RODBC)
library(aqp)
library(plyr)
library(stringr)
library(sp)
library(GSIF)

# ------------------------------------------------------------
# Main tables - NCSS
# ------------------------------------------------------------

#cNCSS <- odbcConnect(dsn="NCSS")
#sqlTables(cNCSS)$TABLE_NAME
# get tables:
#site <- sqlFetch(cNCSS, "site", stringsAsFactors=FALSE)
site <- read.csv("NCSS_Site_Location.csv")
str(site)  # 62,836 profiles
layer <- read.csv("NCSS_Layer.csv")
bdm <- read.csv("NCSS_Bulk_Density_and_Moisture.csv")
summary(bdm$db_od)
## multiple measurements hence average
bdm.0 <- bdm[!duplicated(bdm$labsampnum, fromLast=TRUE),]
summary(bdm.0$db_od)
carb <- read.csv("NCSS_Carbon_and_Extractions.csv")
pH <- read.csv("NCSS_pH_and_Carbonates.csv")
PSDA <- read.csv("NCSS_PSDA_and_Rock_Fragments.csv")
Phosphorus <- read.csv("NCSS_Phosphorus.csv")
water <- read.csv("NCSS_Water_Content.csv")
MajorElements <- read.csv("NCSS_Major_Elements.csv")
Salt <- read.csv("NCSS_Salt.csv")
CEC <- read.csv("NCSS_CEC_and_Bases.csv")
horizons <- plyr::join_all(list(layer, bdm.0, carb, pH, PSDA, Phosphorus, water, MajorElements, Salt, CEC), by="labsampnum")
head(horizons)
nrow(horizons)
## 416,732 obs. of  165 variables

#tax <- sqlFetch(cNCSS, "Taxonomy_New", stringsAsFactors=FALSE)
tax <- read.csv("NCSS_Pedon_Taxonomy.csv")
str(tax)  # 64,050
## Joining by: site_key, latitude_decimal_degrees, longitude_decimal_degrees 
NCSS <- join(site, tax, type="left")
str(NCSS)
taxL <- ifelse(is.na(NCSS$corr_taxsubgrp), paste(NCSS$samp_taxsubgrp), paste(NCSS$corr_taxsubgrp))
s <- summary(as.factor(taxL), maxsum=210)
soiltype <- data.frame(corr_taxsubgrp=attr(s, "names"), count=s)
write.csv(soiltype, "soiltype_count.csv")
ggroup <- data.frame(Great_Group=ifelse(is.na(NCSS$corr_taxgrtgroup), paste(NCSS$samp_taxgrtgroup), paste(NCSS$corr_taxgrtgroup)), Suborder=ifelse(is.na(NCSS$corr_taxsuborder), paste(NCSS$samp_taxsuborder), paste(NCSS$corr_taxsuborder)), Order=ifelse(is.na(NCSS$corr_taxorder), paste(NCSS$samp_taxorder), paste(NCSS$corr_taxorder)))
taxGG <- levels(as.factor(tolower(paste(ggroup[,1], ggroup[,2], ggroup[,3], sep="_"))))
write.csv(taxGG, file="taxGG.csv")

# ------------------------------------------------------------
# NCSS Taxonomy
# ------------------------------------------------------------

NCSS$TIMESTRR <- as.Date(NCSS$site_obsdate, format="%m/%d/%Y")
NCSS$SOURCEID <- make.unique(paste(NCSS$upedonid))
TAXOUSDA.NCSS <- NCSS[,c("SOURCEID","longitude_decimal_degrees", "latitude_decimal_degrees","samp_classification_name","corr_classification_name","corr_taxsuborder","TIMESTRR")]
## rename:
TAXOUSDA.NCSS <- rename(TAXOUSDA.NCSS, replace=c("corr_taxsuborder"="TAXOUSDA", "corr_classification_name"="TAXNUSDA", "longitude_decimal_degrees"="LONWGS84", "latitude_decimal_degrees"="LATWGS84"))
## where 'corr' is missing take at least the 'samp' (described during sampling)
TAXOUSDA.NCSS$SOURCEDB <- "NCSS"
TAXOUSDA.NCSS$TAXNUSDA <- ifelse(is.na(TAXOUSDA.NCSS$TAXNUSDA), paste(TAXOUSDA.NCSS$samp_classification_name), paste(TAXOUSDA.NCSS$TAXNUSDA))  
TAXOUSDA.NCSS <- TAXOUSDA.NCSS[!is.na(TAXOUSDA.NCSS$LONWGS84)&!is.na(TAXOUSDA.NCSS$TAXOUSDA)&nchar(paste(TAXOUSDA.NCSS$TAXOUSDA))>0,c("SOURCEID","SOURCEDB","TIMESTRR","TAXNUSDA","TAXOUSDA","LONWGS84","LATWGS84")]
str(TAXOUSDA.NCSS)
## 31,232 profiles with classification
## convert to SPDF:
coordinates(TAXOUSDA.NCSS) <- ~ LONWGS84 + LATWGS84
## Not necessary any more??
#proj4string(TAXOUSDA.NCSS) <- CRS("+proj=longlat +datum=NAD83")
#TAXOUSDA.NCSS <- spTransform(TAXOUSDA.NCSS, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
proj4string(TAXOUSDA.NCSS) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(TAXOUSDA.NCSS["TAXOUSDA"])
save(TAXOUSDA.NCSS, file="TAXOUSDA.NCSS.rda")

# ------------------------------------------------------------
# Organic soils
# ------------------------------------------------------------

org <- read.csv("NCSS_Organic.csv")
str(org)  # 4362 samples
NCSS.org <- join(join(org, horizons[,c("labsampnum","site_key","pedon_key","hzn_top","hzn_bot","ph_h2o","caco3")], type="left"), NCSS, type="left")
NCSS.org[which(NCSS.org$SOURCEID=="U06MI083-026"),]
PEAT.NCSS <- NCSS.org[,c("SOURCEID", "longitude_decimal_degrees", "latitude_decimal_degrees","corr_taxsuborder", "hzn_top", "hzn_bot", "decomp_state", "oc", "c_tot", "db_od", "ph_h2o", "caco3")]
PEAT.NCSS$thickness <- PEAT.NCSS$hzn_bot-PEAT.NCSS$hzn_top
## correction for OC:
PEAT.NCSS$OC_h <- signif(ifelse(!is.na(PEAT.NCSS$c_tot), ifelse((PEAT.NCSS$ph_h2o > 7)&!is.na(PEAT.NCSS$caco3), PEAT.NCSS$c_tot - .12 * PEAT.NCSS$caco3, PEAT.NCSS$c_tot), PEAT.NCSS$oc), 3)
#summary(as.factor(PEAT.NCSS$corr_taxsuborder))
## Select only profiles with high OC or histosols:
sel <- !is.na(PEAT.NCSS$longitude_decimal_degrees) & (PEAT.NCSS$corr_taxsuborder %in% c("saprists", "hemists", "folists", "fibrists") | (!is.na(PEAT.NCSS$OC_h) & PEAT.NCSS$OC_h > 10 & !is.na(PEAT.NCSS$thickness) & PEAT.NCSS$thickness > 10))
PEAT.NCSS <- PEAT.NCSS[sel,]
str(PEAT.NCSS)
summary(PEAT.NCSS$longitude_decimal_degrees)
depths(PEAT.NCSS) <- SOURCEID ~ hzn_top + hzn_bot
site(PEAT.NCSS) <- ~ longitude_decimal_degrees + latitude_decimal_degrees + corr_taxsuborder
coordinates(PEAT.NCSS) <- ~ longitude_decimal_degrees + latitude_decimal_degrees
#proj4string(PEAT.NCSS) <- CRS("+proj=longlat +datum=NAD83")
#PEAT.NCSS@sp <- spTransform(PEAT.NCSS@sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
proj4string(PEAT.NCSS) <- CRS("+proj=longlat +datum=WGS84")
## convert to a simple table:
PEAT.NCSS.df <- GSIF::as.data.frame(PEAT.NCSS)
PEAT.NCSS.df$Mean_peat_depth_cm <- NA; PEAT.NCSS.df$Max_peat_depth_cm <- NA
write.csv(PEAT.NCSS.df[,c("SOURCEID","longitude_decimal_degrees","latitude_decimal_degrees","corr_taxsuborder","decomp_state_A","Mean_peat_depth_cm", "Max_peat_depth_cm", as.vector(sapply(LETTERS[1:6], function(x){paste0(c("hzn_top_","hzn_bot_","OC_h_","db_od_"), x)})))], file="PEAT.NCSS.csv", na = "")
## "SOURCEID",	"modelling_x", "modelling_y", "TAXOUSDA", "Land_use",	"Mean_peat_depth_cm", "Max_peat_depth_cm", "Upper", "Lower", "BD", "OC"

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

hs <- plyr::join(horizons[,c("labsampnum","site_key","layer_sequence","hzn_desgn","hzn_top","hzn_bot","oc","ph_h2o","ph_kcl","c_tot","caco3","db_od","cec_sum","clay_tot_psa","sand_tot_psa","silt_tot_psa","wpg2")], site[,c("site_key", "usiteid", "site_obsdate", "longitude_decimal_degrees", "latitude_decimal_degrees")], by="site_key", type="left")
hs$HONU <- paste(hs$usiteid, hs$layer_sequence, sep="_")
## Many horizons without depth!
sel0 <- which(hs$hzn_top==0&hs$hzn_bot==0)
horF <- sapply(as.character(hs$HONU[sel0]), function(x){strsplit(x, "_")[[1]][1]})
horS <- as.numeric(sapply(as.character(hs$HONU[sel0]), function(x){strsplit(x, "_")[[1]][2]}))-1
horn <- ifelse(horS<=0|is.na(horS), paste0(horF,'_',1), paste0(horF,'_',horS))
selF <- match(horn, hs$HONU)
hs$hzn_top[sel0] <- hs$hzn_bot[selF]
hs$hzn_bot[which(hs$hzn_bot < hs$hzn_top)] <- NA
hs$DEPTH <- hs$hzn_top + (hs$hzn_bot - hs$hzn_top)/2
hs <- hs[!is.na(hs$longitude_decimal_degrees),]
hs$TIMESTRR <- as.Date(hs$site_obsdate, format="%m/%d/%Y")
hs$SOURCEID <- paste(hs$usiteid)
hs$SAMPLEID <- make.unique(paste(hs$labsampnum))
hs$ORCDRC <- signif(10*ifelse(!is.na(hs$c_tot), ifelse((hs$ph_h2o > 7)&!is.na(hs$caco3), hs$c_tot - .12 * hs$caco3, hs$c_tot), hs$oc), 4)
hs$ORCDRC <- ifelse(hs$ORCDRC < 0, 0, hs$ORCDRC)
#hist(log1p(hs$ORCDRC))
## fix some columns:
hs$wpg2 <- ifelse(hs$wpg2 > 100|hs$wpg2<0 , NA, hs$wpg2)
## merge BulkDens and PSDA and derive GRAVEL content volumetric?
mBD <- hs$wpg2/100 * 2.6 + (1-hs$wpg2/100) * ifelse(is.na(hs$db_od), mean(hs$db_od, na.rm=TRUE), hs$db_od)
hs$CRFVOL <- hs$wpg2*mBD/2.6
hs$BLD <- hs$db_od*1000
## Fix coordinates:
#coordinates(hs) <- ~ longitude_decimal_degrees + latitude_decimal_degrees
#proj4string(hs) <- CRS("+proj=longlat +datum=NAD83")
#hs <- spTransform(hs, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
SPROPS.NCSS <- plyr::rename(hs, replace=c("hzn_top"="UHDICM", "hzn_bot"="LHDICM", "longitude_decimal_degrees"="LONWGS84", "latitude_decimal_degrees"="LATWGS84", "ph_h2o"="PHIHOX", "cec_sum"="CECSUM", "clay_tot_psa"="CLYPPT", "sand_tot_psa"="SNDPPT", "silt_tot_psa"="SLTPPT", "ph_kcl"="PHIKCL"))
SPROPS.NCSS <- SPROPS.NCSS[!is.na(SPROPS.NCSS$DEPTH),c("SOURCEID","SAMPLEID","LONWGS84","LATWGS84","TIMESTRR","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","PHIKCL","BLD","ORCDRC","CECSUM")]
SPROPS.NCSS$SOURCEDB <- "NCSS"
## Fix textures where necessary e.g. negative CLY content:
SPROPS.NCSS[SPROPS.NCSS$SOURCEID=="VPI-Acred-2",]
SPROPS.NCSS[SPROPS.NCSS$SOURCEID=="VPI-Acred-3",]
SPROPS.NCSS$CLYPPT <- ifelse(SPROPS.NCSS$CLYPPT<0, NA, SPROPS.NCSS$CLYPPT)
SPROPS.NCSS$SLTPPT <- ifelse(SPROPS.NCSS$SLTPPT<0, NA, SPROPS.NCSS$SLTPPT)
SPROPS.NCSS$SNDPPT <- ifelse(SPROPS.NCSS$SNDPPT<0, NA, SPROPS.NCSS$SNDPPT)
TexSum <- rowSums(SPROPS.NCSS[,c("CLYPPT","SNDPPT","SLTPPT")], na.rm=TRUE)
hist(TexSum)
selTex <- TexSum > 120 | TexSum < 80
summary(selTex)
SPROPS.NCSS$SLTPPT[selTex] <- NA
## filter pH typos:
SPROPS.NCSS$PHIKCL <- ifelse(SPROPS.NCSS$PHIKCL>11|SPROPS.NCSS$PHIKCL<2.5, NA, SPROPS.NCSS$PHIKCL)
SPROPS.NCSS$PHIHOX <- ifelse(SPROPS.NCSS$PHIHOX>11.5|SPROPS.NCSS$PHIHOX<2.5, NA, SPROPS.NCSS$PHIHOX)
summary(SPROPS.NCSS$PHIHOX)
hist(SPROPS.NCSS$CECSUM)
SPROPS.NCSS$CECSUM <- ifelse(SPROPS.NCSS$CECSUM>500, NA, SPROPS.NCSS$CECSUM)

str(SPROPS.NCSS)
## 348400
SPROPS.NCSS$BLD <- ifelse(SPROPS.NCSS$BLD<100, NA, SPROPS.NCSS$BLD)
summary(SPROPS.NCSS$BLD) ## 1502
summary(SPROPS.NCSS$ORCDRC)  ## 1.6%
## 322,446 records!
save(SPROPS.NCSS, file="SPROPS.NCSS.rda")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- hs[!is.na(hs$hzn_desgn)&nchar(paste(hs$hzn_desgn))>0,]
summary(horizons.s$hzn_desgn)
sel.r <- grep(pattern="^R", as.character(horizons.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", as.character(horizons.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", as.character(horizons.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$hzn_top[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$hzn_bot[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$hzn_top[sel.r3]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$site_key), max, na.rm=TRUE)
names(bdr.d) <- c("site_key", "BDRICM")
BDR.NCSS <- join(bdr.d, site[,c("site_key", "longitude_decimal_degrees", "latitude_decimal_degrees")], type="left")
BDR.NCSS$BDRICM <- ifelse(is.infinite(BDR.NCSS$BDRICM), 250, BDR.NCSS$BDRICM)
BDR.NCSS <- BDR.NCSS[!is.na(BDR.NCSS$BDRICM),]
str(BDR.NCSS)
summary(BDR.NCSS$BDRICM<250)
BDR.NCSS <- plyr::rename(BDR.NCSS, c("site_key"="SOURCEID", "longitude_decimal_degrees"="LONWGS84", "latitude_decimal_degrees"="LATWGS84"))
BDR.NCSS$SOURCEID <- paste0("site_key_", BDR.NCSS$SOURCEID)
BDR.NCSS$SOURCEDB <- "NCSS"
## 1325 points
save(BDR.NCSS, file="BDR.NCSS.rda") 

