## Reading and writing Namibian soil profile data (2157 profiles);
## data received from NAMSOTER 2002 (https://library.wur.nl/isric/fulltext/isricu_i27858_002.pdf)
## by Tom.Hengl@isric.org


library(rgdal)
library(plotKML)
library(plyr)

profs <- read.csv("Namibia_all_profiles.csv")
str(profs)
SITE <- data.frame(profs[!duplicated(profs$PRID),c("PRID","LONG","LATI","DATE","CLAF")])
s <- summary(SITE$CLAF)
soiltype <- data.frame(CLAF=attr(s, "names"), count=s)
write.csv(soiltype, "soiltype_count.csv")
SITE$SOURCEDB = "NAMSOTER"
SITE$TIMESTRR <- as.Date(paste(SITE$DATE), format="%d-%m-%Y")
SITE <- rename(SITE, c("PRID"="SOURCEID", "LONG"="LONWGS84", "LATI"="LATWGS84", "CLAF"="TAXNWRB"))
View(SITE)

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.NAMSOTER <- SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB","TIMESTRR")]
TAXNWRB.NAMSOTER <- TAXNWRB.NAMSOTER[!is.na(TAXNWRB.NAMSOTER$TAXNWRB)&!is.na(TAXNWRB.NAMSOTER$LONWGS84)&nchar(paste(TAXNWRB.NAMSOTER$TAXNWRB))>0,]
str(TAXNWRB.NAMSOTER)
## 2138 profiles
coordinates(TAXNWRB.NAMSOTER) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.NAMSOTER) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.NAMSOTER["TAXNWRB"])
save(TAXNWRB.NAMSOTER, file="TAXNWRB.NAMSOTER.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons <- read.csv("Namibia_all_horizons.csv")
horizons$SAMPLEID <- make.unique(paste(horizons$PRID, horizons$HONU, sep="_"))
horizons <- rename(horizons, c("PRID"="SOURCEID","HBDE"="LHDICM","SDTO"="SNDPPT","STPC"="SLTPPT","CLPC"="CLYPPT","BULK"="BLD","PHAQ"="PHIHOX","TOTC"="ORCDRC","CECS"="CECSUM"))
horizons$ORCDRC <- horizons$ORCDRC*10
summary(horizons$ORCDRC)
## upper horizon boundary missing:
horizons$UHDICM <- NA
horizons$UHDICM <- ifelse(horizons$HONU==1, 0, horizons$UHDICM)
summary(horizons$HONU)
h.lst <- lapply(1:7, function(x){which(horizons$HONU==x)})
for(i in 2:7){
  sel <- match(horizons$SOURCEID[h.lst[[i]]], horizons$SOURCEID[h.lst[[i-1]]])
  horizons$UHDICM[h.lst[[i]]] <- horizons$LHDICM[h.lst[[i-1]]][sel]
}
horizons$LHDICM <- ifelse(is.na(horizons$LHDICM), horizons$UHDICM+50, horizons$LHDICM)
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
View(horizons[,c("SOURCEID","UHDICM","LHDICM","DEPTH","ORCDRC","BLD")])
summary(horizons$DEPTH)
horizons$BLD <- ifelse(horizons$BLD == 9999, NA, horizons$BLD*1000)
summary(horizons$BLD)
## 6208
## negative values for silt?
horizons[horizons$SOURCEID=="NAM-531",]
horizons$SLTPPT <- ifelse(horizons$SLTPPT <0, NA, horizons$SLTPPT)

SPROPS.NAMSOTER <- join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","PHIHOX","ORCDRC","CECSUM","BLD")], SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.NAMSOTER <- SPROPS.NAMSOTER[!is.na(SPROPS.NAMSOTER$LONWGS84) & !is.na(SPROPS.NAMSOTER$LATWGS84) & !is.na(SPROPS.NAMSOTER$DEPTH),]
str(SPROPS.NAMSOTER)
## 5911
save(SPROPS.NAMSOTER, file="SPROPS.NAMSOTER.rda")
plot(SPROPS.NAMSOTER$LONWGS84, SPROPS.NAMSOTER$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- horizons[!is.na(horizons$HODE)&nchar(paste(horizons$HODE))>0,]
sel.r <- grep(pattern="^R", horizons.s$HODE, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$HODE, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", horizons.s$HODE, ignore.case=FALSE, fixed=FALSE)
sel.r4 <- unique(c(which(horizons.s$SOURCEID %in% SITE$SOURCEID[grep(pattern="LEPT", ignore.case=TRUE, paste(SITE$TAXNWRB))]), which(horizons.s$SOURCEID %in% SITE$SOURCEID[grep(pattern="LIT", ignore.case=TRUE, paste(SITE$TAXNWRB))])))
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
horizons.s$BDRICM[sel.r4] <- horizons.s$LHDICM[sel.r4]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.NAMSOTER <- join(SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.NAMSOTER$BDRICM <- ifelse(is.infinite(BDR.NAMSOTER$BDRICM), 250, BDR.NAMSOTER$BDRICM)
BDR.NAMSOTER <- BDR.NAMSOTER[!is.na(BDR.NAMSOTER$BDRICM),]
str(BDR.NAMSOTER)
summary(BDR.NAMSOTER$BDRICM<250)
## 431 points
save(BDR.NAMSOTER, file="BDR.NAMSOTER.rda")