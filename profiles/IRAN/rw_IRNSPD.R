## Reading and writing of IRNSPD data;
## Mohammad, H.B., (2000) Soil resources and use potentiality map of Iran. Soil and Water Research Institute, Teheran, Iran. (1398 profiles)
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(plotKML)
library(sp)

# horizon chemistry:
horizons <- read.table("iran_sdbana.txt", header=FALSE, sep=",", na.strings = c("?","","?.",-2147483647))[,1:12]
names(horizons) <- c("SOURCEID", "SANO", "UHDICM", "LHDICM", "PHIHOX", "EC", "OC", "CACO", "PBS", "SNDPPT", "SLTPPT", "CLYPPT")
horizons$SOURCEID <- as.factor(paste("IR", horizons$SOURCEID, sep="_"))
str(horizons)
horizons <- horizons[!is.na(horizons$SOURCEID),]  # 5078 horizons!
horizons$UHDICM <- ifelse(horizons$UHDICM==-1e+308&horizons$SANO=="A", 0, horizons$UHDICM)
for(j in 5:length(names(horizons))){ horizons[,j] <- ifelse(horizons[,j]<0, NA, horizons[,j]) }

# horizon descriptions: 
iran_sdbhor <- read.table("iran_sdbhor.txt", header=FALSE, sep=",", na.strings = c("?","","?.","??",-2147483647))[,1:8]
names(iran_sdbhor) <- c("SOURCEID", "HRNO", "DESI", "UHDICM", "LHDICM", "COL1", "TEX1", "SANO") 
str(iran_sdbhor)
levels(iran_sdbhor$DESI)
iran_sdbhor$SOURCEID <- as.factor(paste("IR", iran_sdbhor$SOURCEID, sep="_"))

# coordinates:
iran_xy <- read.table("iran_sgdb.txt", header=FALSE, sep=",", na.strings = c("?","","?.","??",-1.00e+308))
names(iran_xy)[1:6] <- c("pnt", "LATWGS84", "LONWGS84", "FAOR", "TAXGWRB", "SOURCEID")
iran_xy$SOURCEID <- as.factor(paste("IR", iran_xy$SOURCEID, sep="_"))
iran_xy <- iran_xy[!is.na(iran_xy$SOURCEID)&!is.na(iran_xy$LATWGS84)&!is.na(iran_xy$LONWGS84),]
# generate horizon / site tables:
site <- iran_xy[!duplicated(iran_xy$SOURCEID),]
site$KEY <- toupper(site$FAOR)
## classification key:
legFAO_90 <- read.csv("FAO_90_names.csv")
legFAO_90$KEY <- toupper(legFAO_90$KEY)
site$TAXNWRB <- join(site, legFAO_90, type="left")$FAO_90_names
str(site)
site$SOURCEDB <- "IRSPDB"

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.Iran <- site[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.Iran <- TAXNWRB.Iran[!is.na(TAXNWRB.Iran$TAXNWRB)&!is.na(TAXNWRB.Iran$LONWGS84)&nchar(paste(TAXNWRB.Iran$TAXNWRB))>0,]
str(TAXNWRB.Iran)
## 3547 profiles
coordinates(TAXNWRB.Iran) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.Iran) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.Iran["TAXNWRB"])
save(TAXNWRB.Iran, file="TAXNWRB.Iran.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

# prepare horizons table:
str(horizons)
horizons$ORCDRC <- horizons$OC*10    ## OC in permilles
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
horizons$SAMPLEID <- make.unique(paste(horizons$SOURCEID, horizons$SANO, sep="_"))

SPROPS.Iran <- join(horizons[!is.na(horizons$DEPTH),c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIHOX","ORCDRC")], site[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.Iran <- SPROPS.Iran[!is.na(SPROPS.Iran$LONWGS84) & !is.na(SPROPS.Iran$LATWGS84) & !is.na(SPROPS.Iran$DEPTH),]
str(SPROPS.Iran)
## 4938
save(SPROPS.Iran, file="SPROPS.Iran.rda")
plot(SPROPS.Iran$LONWGS84, SPROPS.Iran$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- iran_sdbhor[!is.na(iran_sdbhor$DESI)&nchar(paste(iran_sdbhor$DESI))>0,]
sel.r <- grep(pattern="^R", horizons.s$DESI, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$DESI, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", horizons.s$DESI, ignore.case=FALSE, fixed=FALSE)
sel.r4 <- unique(c(which(horizons.s$SOURCEID %in% site$SOURCEID[grep(pattern="LEPT", ignore.case=TRUE, paste(site$TAXNWRB))]), which(horizons.s$SOURCEID %in% site$SOURCEID[grep(pattern="LIT", ignore.case=TRUE, paste(site$TAXNWRB))])))
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
horizons.s$BDRICM[sel.r4] <- horizons.s$LHDICM[sel.r4]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.Iran <- join(site[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.Iran$BDRICM <- ifelse(is.infinite(BDR.Iran$BDRICM), 250, BDR.Iran$BDRICM)
BDR.Iran <- BDR.Iran[!is.na(BDR.Iran$BDRICM),]
str(BDR.Iran)
summary(BDR.Iran$BDRICM<250)
## 427 points
save(BDR.Iran, file="BDR.Iran.rda")

## end of script;