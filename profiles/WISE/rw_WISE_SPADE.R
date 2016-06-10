## Reading and writing of ISRIC WISE (8189 profiles) and SPADE (472) data;
## ISRIC-WISE international soil profile DB [http://www.isric.org/data/isric-wise-international-soil-profile-dataset]
## Tom.Hengl@isric.org

## Download ISRIC WISE database:
download.file("http://www.isric.org/sites/default/files/private/datasets/ISRIC-WISE3r1.zip", destfile=paste(getwd(), "ISRIC-WISE3r1.zip", sep="/"))
unzip("ISRIC-WISE3r1.zip")
## The Soil Profile Analytical Database of Europe [http://eusoils.jrc.ec.europa.eu/projects/spade/spadeM.html]:
download.file("http://eusoils.jrc.ec.europa.eu/projects/spade/Data/SPADE_M_v2.zip", destfile=paste(getwd(), "SPADE_M_v2.zip", sep="/"))
unzip("SPADE_M_v2.zip")

library(RODBC)
library(aqp)
library(plyr)
library(sp)
# define a new function to merge the degree, min, sec columns:
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}
m.lst <- tolower(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

# ------------------------------------------------------------
# Fetch tables
# ------------------------------------------------------------

## import to R:
#cGSPD <- odbcConnect(dsn="ISRIC-WISE")
#odbcGetInfo(cGSPD)
#sqlTables(cGSPD)$TABLE_NAME
# get horizon table:
#HORIZON <- sqlFetch(cGSPD, "WISE3_HORIZON", stringsAsFactors=FALSE, as.is=TRUE)  
HORIZON <- read.csv("WISE3_HORIZON.csv")
str(HORIZON)
# number of samples:
length(HORIZON$ORGC[!is.na(HORIZON$ORGC)])
mean(HORIZON$ORGC, na.rm=TRUE) ## Orgacic Carbon in promilles
summary(HORIZON$VMC1)
#SITE <- sqlFetch(cGSPD, "WISE3_SITE", stringsAsFactors=FALSE) 
SITE <- read.csv("WISE3_SITE.csv", stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

# HORIZON$DEPTH <- (HORIZON$BOTDEP - HORIZON$TOPDEP)/2
SITE$LATSEC <- ifelse(is.na(SITE$LATSEC), 0, SITE$LATSEC)
SITE$LONSEC <- ifelse(is.na(SITE$LONSEC), 0, SITE$LONSEC)
SITE$LATMIN <- ifelse(is.na(SITE$LATMIN), 0, SITE$LATMIN)
SITE$LONMIN <- ifelse(is.na(SITE$LONMIN), 0, SITE$LONMIN)
SITE.s <- subset(SITE, !is.na(SITE$LATDEG)|!is.na(SITE$LONDEG))  # no point in using profiles that have no geo-reference!
# 8189 profiles
# get WGS84 coordinates:
SITE.s$LATWGS84 <- cols2dms(SITE.s$LATDEG, SITE.s$LATMIN, SITE.s$LATSEC, SITE.s$LATIT)
SITE.s$LONWGS84 <- cols2dms(SITE.s$LONDEG, SITE.s$LONMIN, SITE.s$LONSEC, SITE.s$LONGI)
summary(SITE.s$LATWGS84); summary(SITE.s$LONWGS84)
# add columns:
SITE.s$SOURCEID <- as.factor(SITE.s$WISE3_id)
SITE.s$TAXGWRB <- as.factor(SITE.s$WRB2006)
SITE.s$TAXOUSDA <- as.factor(SITE.s$USCL)
SITE.s$TIMESTRR <- as.Date(paste0(SITE.s$DATEYR, "-", ifelse(is.na(SITE.s$DATEMON), m.lst[1], m.lst[SITE.s$DATEMON]), "-1"), format="%Y-%b-%d")
## classification key:
legFAO_90 <- read.csv("cleanup_SU_SYM90.csv")
SITE.s$TAXNWRB <- join(SITE.s, legFAO_90, type="left")$Name
SITE.s$TAXNWRB[SITE.s$TAXNWRB=="#N/A"] <- NA
SITE.s$SOURCEDB = "WISE"
str(SITE.s)

HORIZON$SOURCEID <- as.factor(HORIZON$WISE3_ID)
## rename columns:
HORIZON <- rename(HORIZON, c("TOPDEP"="UHDICM", "BOTDEP"="LHDICM", "ORGC"="ORCDRC", "PHH2O"="PHIHOX", "PHKCL"="PHIKCL", "SAND"="SNDPPT", "SILT"="SLTPPT", "CLAY"="CLYPPT", "GRAVEL"="CRFVOL", "BULKDENS"="BLD", "CECSOIL"="CECSUM"))
## subset to complete data?
sel.h <- !is.na(HORIZON$LHDICM)&!is.na(HORIZON$UHDICM)
HORIZON <- HORIZON[sel.h,]

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

HORIZON$SAMPLEID <- make.unique(paste(HORIZON$SOURCEID, HORIZON$HONU, sep="_"))
HORIZON$SOURCEDB = "WISE"
HORIZON$DEPTH <- HORIZON$UHDICM + (HORIZON$LHDICM - HORIZON$UHDICM)/2
HORIZON$BLD <- HORIZON$BLD * 1000
summary(HORIZON$BLD)
SPROPS.WISE <- join(HORIZON[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CRFVOL","CECSUM","SNDPPT","CLYPPT","BLD","SLTPPT","PHIHOX","PHIKCL","ORCDRC")], SITE.s[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB")])
str(SPROPS.WISE)
summary(SPROPS.WISE$ORCDRC) ## mean=1.4%
summary(SPROPS.WISE$BLD) ## mean=1381
hist(SPROPS.WISE$BLD, breaks=30)
hist(SPROPS.WISE$CECSUM, breaks=30)
hist(SPROPS.WISE$PHIHOX, breaks=30)
## 47,833

save(SPROPS.WISE, file="SPROPS.WISE.rda")
plot(SPROPS.WISE$LONWGS84, SPROPS.WISE$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

HORIZON.s <- HORIZON[!is.na(HORIZON$DESIG)&nchar(paste(HORIZON$DESIG))>0,]
sel.r <- grep(pattern="^R", HORIZON.s$DESIG, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", HORIZON.s$DESIG, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", HORIZON.s$DESIG, ignore.case=FALSE, fixed=FALSE)
HORIZON.s$BDRICM <- NA
HORIZON.s$BDRICM[sel.r] <- HORIZON.s$UHDICM[sel.r]
HORIZON.s$BDRICM[sel.r2] <- HORIZON.s$LHDICM[sel.r2]
HORIZON.s$BDRICM[sel.r3] <- HORIZON.s$LHDICM[sel.r3]
bdr.d <- aggregate(HORIZON.s$BDRICM, list(HORIZON.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.WISE <- join(SITE.s[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.WISE$BDRICM <- ifelse(is.infinite(BDR.WISE$BDRICM), 250, BDR.WISE$BDRICM)
BDR.WISE <- BDR.WISE[!is.na(BDR.WISE$BDRICM),]
str(BDR.WISE)
## 7679 points
save(BDR.WISE, file="BDR.WISE.rda") 


# ------------------------------------------------------------
# SPADE:
# ------------------------------------------------------------

c5GSPD <- odbcConnect(dsn="SPADE_M_v2")
sqlTables(c5GSPD)$TABLE_NAME
# get horizon table:
DAT_HOR <- sqlFetch(c5GSPD, "DAT_HOR", stringsAsFactors=FALSE, as.is=TRUE)
str(DAT_HOR)
mean(DAT_HOR$OC_V, na.rm=TRUE)  ## Organic carbon in %!
DAT_HOR$ORCDRC <- DAT_HOR$OC_V*10
# hist(DAT_HOR$HOR_ID, breaks=25)
summary(as.factor(DAT_HOR$PH_M))
DIC_METH <- sqlFetch(c5GSPD, "DIC_METH")

#DAT_PLOT <- sqlFetch(c5GSPD, "DAT_PLOT", stringsAsFactors=FALSE, as.is=TRUE)
DAT_PLOT <- read.csv("DAT_PLOT.csv")
str(DAT_PLOT)
## rename column:
names(DAT_PLOT)[which(names(DAT_PLOT)=="SOIL_NAME")] = "Local_name"
## classification key:
DIC_SOIL <- read.csv("DIC_SOIL.csv")
DAT_PLOT$TAXNWRB <- join(DAT_PLOT, DIC_SOIL[DIC_SOIL$LEGEND==1990,], type="left", by="SOIL_C")$SOIL_NAME
# clean-up the SITE table:
DAT_PLOT <- DAT_PLOT[!duplicated(DAT_PLOT$PLOT_ID),]

# add columns:
DAT_PLOT$SOURCEID <- as.factor(paste("SPADE", DAT_PLOT$PLOT_ID, sep="_"))
DAT_PLOT$TAXGWRB <- as.factor(DAT_PLOT$SOIL_C)
DAT_PLOT$SOURCEDB = "SPADE"

DAT_HOR$SOURCEID <- as.factor(paste("SPADE", DAT_HOR$PLOT_ID, sep="_"))
DAT_HOR$SILT1_V <- ifelse(is.na(DAT_HOR$SILT1_V), 0, DAT_HOR$SILT1_V)
DAT_HOR$SILT2_V <- ifelse(is.na(DAT_HOR$SILT2_V), 0, DAT_HOR$SILT2_V)
DAT_HOR$SLTPPT <- DAT_HOR$SILT1_V + DAT_HOR$SILT2_V
DAT_HOR$SAND1_V <- ifelse(is.na(DAT_HOR$SAND1_V), 0, DAT_HOR$SAND1_V)
DAT_HOR$SAND2_V <- ifelse(is.na(DAT_HOR$SAND2_V), 0, DAT_HOR$SAND2_V)
DAT_HOR$SAND3_V <- ifelse(is.na(DAT_HOR$SAND3_V), 0, DAT_HOR$SAND3_V)
DAT_HOR$SNDPPT <- DAT_HOR$SAND1_V + DAT_HOR$SAND2_V + DAT_HOR$SAND3_V
DAT_HOR$PHIKCL <- NA
DAT_HOR$PHIKCL[which(DAT_HOR$PH_M %in% "A14")] <- DAT_HOR$PH_V[which(DAT_HOR$PH_M %in% "A14")]
DAT_HOR$PHIHO5 <- NA
DAT_HOR$PHIHO5[which(DAT_HOR$PH_M %in% "A12")] <- DAT_HOR$PH_V[which(DAT_HOR$PH_M %in% "A12")]
# rename some columns:
names(DAT_PLOT)[which(names(DAT_PLOT)=="LON_COOR_V")] <- "LONWGS84"
names(DAT_PLOT)[which(names(DAT_PLOT)=="LAT_COOR_V")] <- "LATWGS84"
names(DAT_HOR)[which(names(DAT_HOR)=="HOR_BEG_V")] <- "UHDICM"
names(DAT_HOR)[which(names(DAT_HOR)=="HOR_END_V")] <- "LHDICM"
names(DAT_HOR)[which(names(DAT_HOR)=="ORGC")] <- "ORCDRC"
names(DAT_HOR)[which(names(DAT_HOR)=="CLAY_V")] <- "CLYPPT"
names(DAT_HOR)[which(names(DAT_HOR)=="GRAV_C")] <- "CRFVOL"
names(DAT_HOR)[which(names(DAT_HOR)=="BD_V")] <- "BLD"
# check the textures:
DAT_HOR[1:3,c("SNDPPT", "SLTPPT", "CLYPPT")]
# subset to complete data:
sel.h <- !is.na(DAT_HOR$LHDICM)&!is.na(DAT_HOR$UHDICM)
sel.s <- !is.na(DAT_PLOT$LONWGS84)&!is.na(DAT_PLOT$LATWGS84)
## depth to R horizon:
sel.r4 <- grep(pattern="^R", DAT_HOR$HOR_NAME, ignore.case=FALSE, fixed=FALSE)
sel.r5 <- unique(c(sel.r4, grep(pattern="999", DAT_HOR$LHDICM, ignore.case=FALSE, fixed=FALSE)))
DAT_HOR$BDRICM <- NA
DAT_HOR$BDRICM[sel.r5] <- DAT_HOR$UHDICM[sel.r5]
bdr.d2 <- aggregate(DAT_HOR$BDRICM, list(DAT_HOR$SOURCEID), max, na.rm=TRUE)
names(bdr.d2) <- c("SOURCEID", "BDRICM")
DAT_PLOTm <- merge(DAT_PLOT, bdr.d2, all.y=FALSE)
DAT_PLOTm$BDRICM <- ifelse(DAT_PLOTm$BDRICM<0, NA, DAT_PLOTm$BDRICM)
summary(DAT_PLOTm$BDRICM) 

profs.f2 <- join(DAT_PLOTm[sel.s,c("SOURCEID", "LONWGS84", "LATWGS84", "TAXGWRB","BDRICM")], DAT_HOR[sel.h,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT", "SLTPPT", "CRFVOL","PHIHO5","PHIKCL","BLD","ORCDRC")], type='inner')
depths(profs.f2) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f2) <- ~ LONWGS84 + LATWGS84 + TAXGWRB + BDRICM
## check if there are still some duplicates
profs.f2@site$SOURCEID[duplicated(profs.f2@site$SOURCEID)]
## spatial duplicates:
sp2 <- profs.f2@site
coordinates(sp2) <- ~ LONWGS84 + LATWGS84
sel2 <- sp::zerodist(sp2)
str(sel2)
## 12 duplicate points

# prepare a SITE and hor table table:
wsp9 <- list(sites=profs.f2@site, horizons=profs.f2@horizons)
str(wsp9)
wsp9$sites$TAXGWRB <- as.character(wsp9$sites$TAXGWRB)
lapply(as.list(wsp9$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp9$horizons), function(x){sum(!is.na(x))})
save(wsp9, file="../wsp9.rda", compress="xz")

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

## Universal soil classification data:
USCL <- read.csv("USCL_points.csv", stringsAsFactors=FALSE)
USCL <- rename(USCL[,c("WRB2006_subgroup","WISE_id","LONdd","LATdd","SOURCE_ID")], c("WRB2006_subgroup"="TAXNWRB", "WISE_id"="SOURCEID", "LONdd"="LONWGS84", "LATdd"="LATWGS84", "SOURCE_ID"="SOURCEDB"))

TAXNWRB.WISE <- rbind.fill(list(SITE.s[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")], DAT_PLOT[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB")], USCL))
TAXNWRB.WISE <- TAXNWRB.WISE[!is.na(TAXNWRB.WISE$TAXNWRB)&!is.na(TAXNWRB.WISE$LONWGS84)&nchar(paste(TAXNWRB.WISE$TAXNWRB))>0,]
str(TAXNWRB.WISE)
## 16,475 profiles
coordinates(TAXNWRB.WISE) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.WISE) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.WISE["TAXNWRB"])
save(TAXNWRB.WISE, file="TAXNWRB.WISE.rda")

TAXOUSDA.WISE <- SITE.s[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXOUSDA")]
TAXOUSDA.WISE <- TAXOUSDA.WISE[!is.na(TAXOUSDA.WISE$TAXOUSDA)&nchar(paste(TAXOUSDA.WISE$TAXOUSDA))>0,]
str(TAXOUSDA.WISE)
## 3491 profiles
coordinates(TAXOUSDA.WISE) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.WISE) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.WISE["TAXOUSDA"])
save(TAXOUSDA.WISE, file="TAXOUSDA.WISE.rda")

## export taxonomy data for correlation:
WISE_tax <- SITE.s[,c("SOURCEID", "LATWGS84", "LONWGS84", "TAXNWRB", "TAXOUSDA", "LFORM", "LANDUS")]
summary(as.factor(WISE_tax$TAXOUSDA))
save(WISE_tax, file="WISE_tax.rda")

# end of script;