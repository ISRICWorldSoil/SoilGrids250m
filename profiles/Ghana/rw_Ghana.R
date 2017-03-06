## Ghana legacy data; contact: Kwabena A. Nketia <kanj241@gmail.com>
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(sp)
library(plotKML)
library(utils)
load(".RData")

locs <- read.csv("Geopoints.csv")
locs$SOURCEID = make.unique(paste(locs$Pr_Key))
locs = rename(locs, c("X_LonDD"="LONWGS84", "Y_LatDD"="LATWGS84"))

SITE <- read.csv("Profiles.csv")
str(SITE)
## 1403 points
SITE$SOURCEID = make.unique(paste(SITE$Pr_Key))
SITE$TIMESTRR = as.Date(SITE$Obs_Year, format="%d-%b-%y")
SITE$SOURCEDB = "GhanaDB"
SITE$TAXNWRB <- ifelse(is.na(SITE$WRB06), paste(SITE$Ori_IntCls), paste(SITE$WRB06))
summary(as.factor(SITE$TAXNWRB))

horizons = read.csv("Layers.csv")
horizons$SAMPLEID = make.unique(paste(horizons$LayerID))
horizons$LowDpth = as.numeric(paste(horizons$LowDpth))
horizons$UpDpth = as.numeric(paste(horizons$UpDpth))
for(j in c("LowDpth","UpDpth","PH","TotOrgC","CECsoil","Silt","Sand","Clay","CECsoil")){ horizons[,j] = ifelse(horizons[,j]< -20, NA, horizons[,j]) }
horizons.s = rename(horizons, c("Pr_Key"="SOURCEID", "UpDpth"="UHDICM", "LowDpth"="LHDICM", "Sand"="SNDPPT","Silt"="SLTPPT","Clay"="CLYPPT","PH"="PHIHOX","CECsoil"="CECSUM"))
summary(horizons.s$TotOrgC)
horizons.s$ORCDRC <- horizons.s$TotOrgC*10
summary(horizons.s$ORCDRC)
summary(horizons.s$pH_Method)
horizons.s$PHIKCL = ifelse(horizons.s$pH_Method=="pHKCl_1:1", horizons.s$PHIHOX, NA)
horizons.s$PHIHOX = ifelse(horizons.s$pH_Method=="pHKCl_1:1", NA, horizons.s$PHIHOX)
horizons.s$PHIHOX = ifelse(horizons.s$PHIHOX<2|horizons.s$PHIHOX>11, NA, horizons.s$PHIHOX)
summary(horizons.s$PHIHOX)
summary(horizons.s$PHIKCL)

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.Ghana <- join(SITE[,c("SOURCEID","SOURCEDB","TAXNWRB")], locs[,c("SOURCEID","LONWGS84","LATWGS84")], match="first")
TAXNWRB.Ghana <- TAXNWRB.Ghana[!is.na(TAXNWRB.Ghana$TAXNWRB)&!TAXNWRB.Ghana$TAXNWRB=="NA"&!is.na(TAXNWRB.Ghana$LONWGS84)&nchar(paste(TAXNWRB.Ghana$TAXNWRB))>0,]
str(TAXNWRB.Ghana)
## 505 profiles
coordinates(TAXNWRB.Ghana) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.Ghana) <- "+proj=longlat +datum=WGS84"
#plotKML(TAXNWRB.Ghana["TAXNWRB"])
save(TAXNWRB.Ghana, file="TAXNWRB.Ghana.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons.s$SLTPPT <- ifelse(is.na(horizons.s$SLTPPT), 100-(horizons.s$CLYPPT+horizons.s$SNDPPT), horizons.s$SLTPPT)
horizons.s$SNDPPT <- ifelse(is.na(horizons.s$SNDPPT), 100-(horizons.s$CLYPPT+horizons.s$SLTPPT), horizons.s$SNDPPT)
horizons.s$DEPTH <- horizons.s$UHDICM + (horizons.s$LHDICM - horizons.s$UHDICM)/2
horizons.s$SAMPLEID = horizons.s$LayerID

SPROPS.Ghana <- join(horizons.s[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","PHIHOX","PHIKCL","ORCDRC","CECSUM")], join(SITE[,c("SOURCEID","SOURCEDB","TIMESTRR")], locs[,c("SOURCEID","LONWGS84","LATWGS84")]), type="left")
SPROPS.Ghana <- SPROPS.Ghana[!is.na(SPROPS.Ghana$LONWGS84) & !is.na(SPROPS.Ghana$LATWGS84) & !is.na(SPROPS.Ghana$DEPTH),]
str(SPROPS.Ghana)
summary(SPROPS.Ghana$ORCDRC) ## mean = 1.7%
hist(SPROPS.Ghana$CLYPPT, col="gray")
## 9245
summary(SPROPS.Ghana$TIMESTRR)
save(SPROPS.Ghana, file="SPROPS.Ghana.rda")
plot(SPROPS.Ghana$LONWGS84, SPROPS.Ghana$LATWGS84, pch="+")
