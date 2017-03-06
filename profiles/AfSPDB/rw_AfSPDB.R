# Reading and writing of the African SPDB (ca 18,000 profiles);
# Download data from [http://www.isric.org/data/africa-soil-profiles-database]
# Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(sp)
library(zoo)
library(plotKML)
library(tools)

Layers <- read.csv("AfSP01301Qry_Layers.csv", stringsAsFactors=FALSE, fileEncoding="latin1")
Profiles <- read.csv("AfSP01301Qry_Profiles.csv", stringsAsFactors=FALSE, fileEncoding="latin1")
## select columns of interest:
s.lst <- c("ProfileID", "Reliab", "X_LonDD", "Y_LatDD", "XYAccur", "T_Year", "WRB06", "WRB06rg", "FAO88", "USDA", "LocalCls", "Location", "Drain", "RockDpth")
h.lst <- c("ProfileID", "LayerID", "LayerNr", "UpDpth", "LowDpth", "HorDes", "ColorM", "ColorD", "FldTxtr", "CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "EC", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC", "ExCa", "ExNa", "ExMg", "ExK", "ExBases", "ExAl", "ExAcid")

## clean up horizons:
Layers$LayerNr <- as.integer(Layers$LayerNr)
Layers$FldTxtr[which(Layers$FldTxtr=="NA")] <- NA
Layers$HorDes[which(Layers$HorDes=="NA")] <- NA
Layers$ColorM[which(Layers$ColorM=="NA")] <- NA
Layers$ColorD[which(Layers$ColorD=="NA")] <- NA
Layers$ColorM[which(Layers$ColorM=="<Null>")] <- NA
Layers$ColorD[which(Layers$ColorD=="<Null>")] <- NA
Layers$FldTxtr[which(Layers$FldTxtr=="NA")] <- NA
# fix the Munsell color codes?
# replace "-9999":
for(j in c("CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "EC", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC", "ExCa", "ExNa", "ExMg", "ExK", "ExBases", "ExAl", "ExAcid")){
  if(is.numeric(Layers[,j])){
    Layers[,j] <- signif(ifelse(Layers[,j]==-9999, NA, Layers[,j]), 2)
  } else {
    Layers[,j] <- ifelse(Layers[,j]==-9999, NA, Layers[,j])
  }
}
horizons <- rename(Layers[,c("ProfileID","LayerID","UpDpth","LowDpth","CfPc","Sand","Silt","Clay","BlkDens","PHH2O","PHKCl","OrgC","CecSoil")], c("ProfileID"="SOURCEID", "LayerID"="SAMPLEID", "UpDpth"="UHDICM", "LowDpth"="LHDICM", "CfPc"="CRFVOL", "Sand"="SNDPPT", "Silt"="SLTPPT", "Clay"="CLYPPT", "BlkDens"="BLD", "PHH2O"="PHIHOX", "PHKCl"="PHIKCL", "OrgC"="ORCDRC", "CecSoil"="CECSUM"))
horizons$SAMPLEID <- make.unique(as.character(horizons$SAMPLEID))
horizons$LHDICM <- as.numeric(horizons$LHDICM)
horizons$CRFVOL <- as.numeric(horizons$CRFVOL)
## 79,531 layer

## hits on visual check:
## 39 negative UHDICM values (likely 'mulch' layer records as LHDICM = 0)
hist(log1p(horizons$UHDICM))
str(horizons[which(horizons$UHDICM<0),])
## 7 UHDICM values < 800cm 
str(horizons[which(horizons$UHDICM>800),])
## 16 LHDICM values < 800cm 
str(horizons[which(horizons$LHDICM>800),])
## 2 obs with deeper 'upper' horizons
str(horizons[which(horizons$LHDICM-horizons$UHDICM<=0),])
##13 obs with messy texture fractions  
TexSum <- rowSums(horizons[,c("CLYPPT","SNDPPT","SLTPPT")], na.rm=TRUE)
hist(TexSum)
selTex <- TexSum > 110 | TexSum < 90
str(horizons[which(TexSum>0&TexSum<90),])
## 61 obs above pH 10.5 - rare potential for EXTREME sodic soils (but above 10.5 is questionable calibration?? Or an arid mineral spring ;P) check spatial location for clustering.
str(horizons[which(horizons$PHIHOX>=10.5),])
## as above but 46 obs. Check spatial location and with soil chemist.
str(horizons[which(horizons$PHIKCL>=9.5),])
#94 obs with 0 OC but values for CEC 
str(horizons[which(horizons$ORCDRC<=0),])
View(horizons[which(horizons$ORCDRC<=0),])
##3 obs of CEC = 0
str(horizons[which(horizons$CECSUM<=0),])
## issues with BLD in gm per cm3 rather than kg per m3? - check

## clean up sites:
sites <- Profiles[!duplicated(Profiles$ProfileID),]
sites <- sites[!is.na(sites$X_LonDD)&!is.na(sites$Y_LatDD), s.lst]
sites$XYAccur <- ifelse(sites$XYAccur==-9999, NA, sites$XYAccur)
sites$T_Year <- ifelse(sites$T_Year==-9999, NA, sites$T_Year)
sites$T_Year <- as.Date(paste(sites$T_Year), format="%Y")
sites$WRB06[which(sites$WRB06=="NA")] <- NA
sites$WRB06rg[which(sites$WRB06rg=="NA")] <- NA
sites$FAO88[which(sites$FAO88=="NA")] <- NA
sites$USDA[which(sites$USDA=="NA")] <- NA
sites$LocalCls[which(sites$LocalCls=="NA")] <- NA
sites$Location[which(sites$Location=="NA")] <- NA
sites$Drain[which(sites$Drain=="NA")] <- NA
sites$SOURCEDB = "AfSPDB"
sites <- rename(sites, c("ProfileID"="SOURCEID", "X_LonDD"="LONWGS84", "Y_LatDD"="LATWGS84", "T_Year"="TIMESTRR", "WRB06"="TAXNWRB", "USDA"="TAXNUSDA"))

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

## Filtered by Johan Leernaars:
RockDpth.sel <- grep(pattern=">", sites$RockDpth, fixed=FALSE)
rd = as.numeric(sapply(as.character(sites$RockDpth[RockDpth.sel]), function(x){strsplit(x, ">")[[1]][2]}))
sel.n <- unique(c(grep(pattern="LEPT", ignore.case=TRUE, paste(sites$TAXNWRB)), grep(pattern="LEPT", ignore.case=TRUE, paste(sites$TAXNWRB))))
sites$BDRICM <- as.numeric(as.character(sites$RockDpth))*100
sites$BDRICM <- ifelse(is.na(sites$BDRICM), 250, sites$BDRICM)
sites[RockDpth.sel, "BDRICM"] <- 100 * rd
sites[-sel.n, "BDRICM"] <- ifelse(sites[-sel.n, "BDRICM"]<100, NA, sites[-sel.n, "BDRICM"])
summary(sites$BDRICM)

BDR.AfSPDB <- sites[!is.na(sites$BDRICM),c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","BDRICM")]
summary(BDR.AfSPDB$BDRICM<250)
## 2382 points
str(BDR.AfSPDB)
save(BDR.AfSPDB, file="BDR.AfSPDB.rda")

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXOUSDA.AfSPDB <- sites[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNUSDA")]
TAXOUSDA.AfSPDB <- TAXOUSDA.AfSPDB[!is.na(TAXOUSDA.AfSPDB$TAXNUSDA)&!is.na(TAXOUSDA.AfSPDB$LONWGS84)&nchar(paste(TAXOUSDA.AfSPDB$TAXNUSDA))>0,]
str(TAXOUSDA.AfSPDB)
## 2577 profiles
coordinates(TAXOUSDA.AfSPDB) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.AfSPDB) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.AfSPDB["TAXNUSDA"])
save(TAXOUSDA.AfSPDB, file="TAXOUSDA.AfSPDB.rda")

TAXNWRB.AfSPDB <- sites[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.AfSPDB <- TAXNWRB.AfSPDB[!is.na(TAXNWRB.AfSPDB$TAXNWRB)&!is.na(TAXNWRB.AfSPDB$LONWGS84)&nchar(paste(TAXNWRB.AfSPDB$TAXNWRB))>0,]
TAXNWRB.AfSPDB$TAXNWRB <- iconv(paste(TAXNWRB.AfSPDB$TAXNWRB), "latin1", "UTF-8", "")
str(TAXNWRB.AfSPDB)
## 3360 profiles
coordinates(TAXNWRB.AfSPDB) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.AfSPDB) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.AfSPDB["TAXNWRB"])
save(TAXNWRB.AfSPDB, file="TAXNWRB.AfSPDB.rda")

## visual check 
View(TAXNWRB.AfSPDB$TAXNWRB)
## 22 obs with potential non-ascii characters
print(showNonASCII(TAXNWRB.AfSPDB$TAXNWRB))
## no other detected issues 
View(TAXOUSDA.AfSPDB$TAXNUSDA)
print(showNonASCII(as.character(TAXOUSDA.AfSPDB$TAXNUSDA)))

# ------------------------------------------------------------
# Organic soils
# ------------------------------------------------------------

sel.hist <- which(horizons$SOURCEID %in% as.character(TAXNWRB.AfSPDB$SOURCEID[grep(TAXNWRB.AfSPDB$TAXNWRB, pattern="ist")]))
sel.OC <- which(horizons$ORCDRC > 120)
PEAT.AfSPDB <- join(horizons[c(sel.hist,sel.OC),], as.data.frame(TAXNWRB.AfSPDB), type="left")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.AfSPDB <- plyr::join(horizons, sites[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], type="left")
SPROPS.AfSPDB$DEPTH <- SPROPS.AfSPDB$UHDICM + (SPROPS.AfSPDB$LHDICM - SPROPS.AfSPDB$UHDICM)/2
SPROPS.AfSPDB <- SPROPS.AfSPDB[!is.na(SPROPS.AfSPDB$LONWGS84) & !is.na(SPROPS.AfSPDB$LATWGS84) & !SPROPS.AfSPDB$LONWGS84==0 & !is.na(SPROPS.AfSPDB$DEPTH),]
SPROPS.AfSPDB$BLD = SPROPS.AfSPDB$BLD * 1000
str(SPROPS.AfSPDB)
## 73,332
save(SPROPS.AfSPDB, file="SPROPS.AfSPDB.rda")

## v1.3
profsv13 = SPROPS.AfSPDB[which(SPROPS.AfSPDB$SOURCEID %in% Profiles$ProfileID[Profiles$DbVersion==1.3]),]
saveRDS(profsv13, "AfPDB_profsv13.rds")
plot(profsv13[,c("LONWGS84","LATWGS84")])
pnts.xy = profsv13
#pnts.xy = pnts.xy[!duplicated(pnts.xy$SOURCEID),]
pnts.xy = pnts.xy[!is.na(pnts.xy$ORCDRC),]
coordinates(pnts.xy) = ~ LONWGS84 + LATWGS84
proj4string(pnts.xy) = "+proj=longlat +datum=WGS84"
unlink("AfPDB_profsv13.shp")
writeOGR(pnts.xy, "AfPDB_profsv13.shp", "AfPDB_profsv13", "ESRI Shapefile")
plotKML(pnts.xy)
