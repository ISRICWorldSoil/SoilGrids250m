# title         : rw_NCSS.R
# purpose       : Reading and writing of NSCD data;
# reference     : NCSS Characterization Database [http://ncsslabdatamart.sc.egov.usda.gov/] (25,802 profiles);
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : "Repository2010.mdb" MS Access dbs
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : LARGE dataset!

## Download the database:
# download.file("http://globalsoilmap.net/data/Repository2010.7z", destfile=paste(getwd(),"Repository2010.7z",sep="/"))

library(RODBC)
library(aqp)
library(plyr)
library(sp)
# define a new function to merge the degree, min, sec columns:
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}
load("Repository2010.RData")

# ------------------------------------------------------------
# Fetch tables - NCSS
# ------------------------------------------------------------

cNCSS <- odbcConnect(dsn="NCSS")
sqlTables(cNCSS)$TABLE_NAME
# get tables:
site <- sqlFetch(cNCSS, "site", stringsAsFactors=FALSE)
str(site)  # 37,333 profiles
pedon <- sqlFetch(cNCSS, "pedon", stringsAsFactors=FALSE)
str(pedon)
pedon$observation_date <- as.character(as.Date(pedon$observation_date))
# has to remove Date format otherwise it gets problems with NA's
PSDA <- sqlFetch(cNCSS, "PSDA and Rock Fragments", stringsAsFactors=FALSE)   ## gravel content ("wpg2") is in mass percentage and needs to be converted to volume %
str(PSDA)
Organic <- sqlFetch(cNCSS, "Organic", stringsAsFactors=FALSE)
str(Organic)
CEC <- sqlFetch(cNCSS, "CEC and Bases", stringsAsFactors=FALSE)
str(CEC)
Carbon <- sqlFetch(cNCSS, "Carbon and Extractions", stringsAsFactors=FALSE)
str(Carbon)
BulkDens <- sqlFetch(cNCSS, "Bulk Density and Moisture", stringsAsFactors=FALSE)
str(BulkDens)  ## we need "Bulk Density, <2mm Fraction, Ovendry"
pH <- sqlFetch(cNCSS, "pH and Carbonates", stringsAsFactors=FALSE)
str(pH)
layer <- sqlFetch(cNCSS, "layer", stringsAsFactors=FALSE)
str(layer)
tax <- sqlFetch(cNCSS, "Taxonomy_New", stringsAsFactors=FALSE)
tax$last_correlated_date <- as.character(as.Date(tax$last_correlated_date))
Phosphorus <- sqlFetch(cNCSS, "Phosphorus", stringsAsFactors=FALSE)
MajorElements <- sqlFetch(cNCSS, "Major Elements", stringsAsFactors=FALSE)
str(MajorElements)
Salt <- sqlFetch(cNCSS, "Salt", stringsAsFactors=FALSE)
str(Salt)
TaxonomyE <- sqlFetch(cNCSS, "Taxonomy_Error", stringsAsFactors=FALSE)
str(TaxonomyE)
TaxonomyE$last_correlated_date <- as.character(as.Date(TaxonomyE$last_correlated_date))
procs <- sqlFetch(cNCSS, "analysis_procedure", stringsAsFactors=FALSE)
# save(list=c("Organic", "site", "pedon", "PSDA", "CEC", "Carbon", "BulkDens", "pH", "layer", "tax", "Phosphorus", "MajorElements", "Salt", "TaxonomyE"), file="cNCSS.RData")


# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

# add missing columns:
layer$DEPTH <- layer$hzn_top + (layer$hzn_bot - layer$hzn_top)/2
# mask out profiles with missing coordinates:
site <- site[!is.na(site$longitude_degrees)&!is.na(site$latitude_degrees),]
site$longitude_it <- ifelse(site$longitude_direction=="west", "W", "E")
site$latitude_it <- ifelse(site$latitude_direction=="north"|site$latitude_direction=="North", "N", "S")
site$LAT <- cols2dms(site$latitude_degrees, site$latitude_minutes, site$latitude_seconds, site$latitude_it)
site$LON <- cols2dms(site$longitude_degrees, site$longitude_minutes, site$longitude_seconds, site$longitude_it)
summary(site$LAT) # 253 NA's
summary(site$LON)
summary(site$latitude_seconds)
str(site)

# ------------------------------------------------------------
# create the horizon and site tables (takes time!)
# ------------------------------------------------------------
 
h1 <- merge(PSDA, CEC[,!(names(CEC) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h2 <- merge(Organic[,!(names(Organic) %in% c("result_source_key", "prep_code", "c_tot", "oc", "n_tot", "c_n_ra"))], Carbon, by=c("natural_key"), all=TRUE)
h3 <- merge(h1, h2[,!(names(h2) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h4 <- merge(h3, layer, by=c("natural_key"), all=TRUE)
h5 <- merge(h4[,!(names(h4) %in% c("db_od"))], BulkDens[,!(names(BulkDens) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
horizon <- merge(h5[,!(names(h5) %in% c("ph_h2o"))], pH[,!(names(pH) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
# names(horizon)
## fix some columns:
horizon$wpg2 <- ifelse(horizon$wpg2 > 100|horizon$wpg2 <0 , NA, horizon$wpg2)
## merge BulkDens and PSDA and derive GRAVEL content:
mBD <- (horizon$wpg2/100 * 2.6 + (1-horizon$wpg2/100) * horizon$db_od)
horizon$GRAVEL <- horizon$wpg2*mBD/2.6
## check visually:
# plot(y=horizon$GRAVEL, x=horizon$wpg2, xlim=c(0,100), ylab="GRAVEL (vol %)", xlab="GRAVEL (mass %)", main="NCSS (250K measurements)", pch="+", cex=.7)

pedon$site_key <- as.factor(pedon$site_key)
# there are also pedons with multiple site IDs!?
s1 <- merge(pedon, tax[,!(names(tax) %in% c("natural_key"))], by=c("pedon_key"), all=TRUE)
# summary(s1$site_key)  
s1 <- s1[!duplicated(s1$site_key),]
str(s1)
# merge the site table:
site$site_key <- as.factor(site$site_key)
site.s <- merge(site, s1, by=c("site_key"), all.x=TRUE)
site.s <- site.s[!is.na(site.s$LAT)&!is.na(site.s$LON), c("site_key","user_site_id","LAT","LON","horizontal_datum_name","observation_date","sampled_taxon_name","correlated_taxon_name","correlated_taxon_kind","correlated_class_name","SSL_class_name")]
str(site.s)

# subset tables:
hs <- subset(horizon[,c("site_key", "natural_key", "layer_sequence", "hzn_desgn", "hzn_bot", "hzn_top", "hzn_vert_subdvn", "clay_tot_psa", "sand_tot_psa", "silt_tot_psa", "pyr_col", "wpg2", "caco3", "ph_hist", "oc", "c_tot", "base_sum", "cec_sum", "cec_nh4", "ph_h2o", "ph_kcl", "db_od", "hzn_desgn")], !is.na(horizon$DEPTH)&!is.na(horizon$site_key)&horizon$layer_sequence<15)
# strange - there are layer_sequence numbers up to 259?!
str(hs)
# summary(as.factor(horizon$layer_sequence))
hs$site_key <- as.factor(hs$site_key)
# remove duplicate site keys (there are many!!):
hs$layer_ID <- as.factor(paste(hs$site_key, hs$layer_sequence, sep="_"))
hs <- hs[!duplicated(hs$layer_ID),]
# remove IDs that do not exist in the site table:
sel <- hs$site_key %in% site.s$site_key
length(levels(as.factor(hs$layer_ID)))
mean(hs$oc, na.rm=TRUE) ## Organic Carbon in percent!
## estimate OC by correcting for CaCO3:
hs$ORCDRC <- signif(10*ifelse(!is.na(hs$c_tot), ifelse((hs$ph_h2o > 7)&!is.na(hs$caco3), hs$c_tot - .12 * hs$caco3, hs$c_tot), hs$oc), 4)
hs$ORCDRC <- ifelse(hs$ORCDRC < 0, 0, hs$ORCDRC) 
## stats:
# summary(hs$ORCDRC)


# ------------------------------------------------------------
# SoilProfileCollection
# ------------------------------------------------------------

NCSS_all <- list(sites=site.s, horizons=hs[sel,])
str(NCSS_all)
# specify classes:
NCSS_all$sites$user_site_id <- as.factor(NCSS_all$sites$user_site_id)
NCSS_all$sites$horizontal_datum_name <- as.factor(NCSS_all$sites$horizontal_datum_name)
NCSS_all$sites$correlated_taxon_name <- as.factor(NCSS_all$sites$correlated_taxon_name)
NCSS_all$sites$sampled_taxon_name <- as.factor(NCSS_all$sites$sampled_taxon_name)
NCSS_all$sites$correlated_taxon_kind <- as.factor(NCSS_all$sites$correlated_taxon_kind)
NCSS_all$sites$correlated_class_name <- as.factor(NCSS_all$sites$correlated_class_name)
NCSS_all$sites$SSL_class_name <- as.factor(NCSS_all$sites$SSL_class_name)
NCSS_all$horizons$hzn_desgn <- as.factor(NCSS_all$horizons$hzn_desgn)
NCSS_all$horizons$natural_key <- as.factor(NCSS_all$horizons$natural_key)
NCSS_all$horizons$pyr_col <- NULL
# round up the numbers:
NCSS_all$horizons$ORCDRC <- round(NCSS_all$horizons$ORCDRC, 1)
NCSS_all$horizons$clay_tot_psa <- round(NCSS_all$horizons$clay_tot_psa, 0)
NCSS_all$horizons$sand_tot_psa <- round(NCSS_all$horizons$sand_tot_psa, 0)
NCSS_all$horizons$silt_tot_psa <- round(NCSS_all$horizons$silt_tot_psa, 0)
NCSS_all$horizons$oc <- round(NCSS_all$horizons$oc, 2)
NCSS_all$horizons$c_tot <- round(NCSS_all$horizons$c_tot, 2)
NCSS_all$horizons$base_sum <- round(NCSS_all$horizons$base_sum, 1)
NCSS_all$horizons$caco3 <- round(NCSS_all$horizons$caco3, 2)
NCSS_all$horizons$cec_sum <- round(NCSS_all$horizons$cec_sum, 1)
NCSS_all$horizons$cec_nh4 <- round(NCSS_all$horizons$cec_nh4, 1)
NCSS_all$horizons$ph_h2o <- round(NCSS_all$horizons$ph_h2o, 1)
NCSS_all$horizons$ph_kcl <- round(NCSS_all$horizons$ph_kcl, 1)
NCSS_all$horizons$db_od <- round(NCSS_all$horizons$db_od, 2)
save(NCSS_all, file="NCSS_all.rda", compress="xz")

## Prepare data for WSP:
h.s <- NCSS_all$horizons[,c("site_key","hzn_top","hzn_bot","sand_tot_psa","silt_tot_psa","clay_tot_psa","wpg2","ORCDRC","ph_h2o","ph_kcl","db_od","hzn_desgn")]
s.s <- NCSS_all$sites[,c("site_key","user_site_id","LON","LAT","correlated_class_name","observation_date")]
h.s <- merge(h.s, s.s[,c("site_key", "user_site_id")], by="site_key")
# reproject to WGS84:
coordinates(s.s) <- ~LON+LAT
proj4string(s.s) <- CRS("+proj=longlat +datum=NAD83")
s.s <- data.frame(spTransform(s.s, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
# rename columns:
str(s.s)
names(s.s)[which(names(s.s)=="LON")] <- "LONWGS84"
names(s.s)[which(names(s.s)=="LAT")] <- "LATWGS84"
names(s.s)[which(names(s.s)=="user_site_id")] <- "SOURCEID"
names(s.s)[which(names(s.s)=="correlated_class_name")] <- "TAXNUSDA"
s.s$TIMESTRR <- as.Date(s.s$observation_date, format="%Y-%m-%d")

# horizons table / rename columns:
str(h.s)
names(h.s)[which(names(h.s)=="user_site_id")] <- "SOURCEID"
names(h.s)[which(names(h.s)=="hzn_top")] <- "UHDICM"
names(h.s)[which(names(h.s)=="hzn_bot")] <- "LHDICM"
names(h.s)[which(names(h.s)=="clay_tot_psa")] <- "CLYPPT"
names(h.s)[which(names(h.s)=="sand_tot_psa")] <- "SNDPPT"
names(h.s)[which(names(h.s)=="silt_tot_psa")] <- "SLTPPT"
names(h.s)[which(names(h.s)=="wpg2")] <- "CRFVOL"
names(h.s)[which(names(h.s)=="ph_h2o")] <- "PHIHO5"
names(h.s)[which(names(h.s)=="ph_kcl")] <- "PHIKCL"
names(h.s)[which(names(h.s)=="db_od")] <- "BLD"
summary(h.s$PHIHO5)
## Depth to bedrock [http://www.swac.umn.edu/classes/soil2125/doc/s3chap2.htm]:
sel.r <- grep("r", h.s$hzn_desgn, ignore.case=TRUE, fixed=FALSE)
h.s$BDRICM <- NA
h.s$BDRICM[sel.r] <- pmax(h.s$UHDICM[sel.r], h.s$LHDICM[sel.r], na.rm=TRUE)
bdr.d <- aggregate(h.s$BDRICM, list(h.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
s.sm <- merge(s.s, bdr.d, all.y=FALSE)
s.sm$BDRICM <- ifelse(s.sm$BDRICM<0|s.sm$BDRICM>400, NA, s.sm$BDRICM)
summary(s.sm$BDRICM)
## 15 profiles need to be fixed manually :((
s.sm$BDRICM[which(s.sm$site_key==25915)] <- 70
s.sm$BDRICM[which(s.sm$site_key==1068)] <- 114
s.sm$BDRICM[which(s.sm$site_key==1056)] <- 157
s.sm$BDRICM[which(s.sm$site_key==1039)] <- 229
s.sm$BDRICM[which(s.sm$site_key==5614)] <- 158
s.sm$BDRICM[which(s.sm$site_key==5615)] <- 65
s.sm$BDRICM[which(s.sm$site_key==5720)] <- 142
s.sm$BDRICM[which(s.sm$site_key==5726)] <- 152
s.sm$BDRICM[which(s.sm$site_key==7020)] <- 152
s.sm$BDRICM[which(s.sm$site_key==9928)] <- 74
s.sm$BDRICM[which(s.sm$site_key==9272)] <- 107
s.sm$BDRICM[which(s.sm$site_key==16289)] <- 160
s.sm$BDRICM[which(s.sm$site_key==16290)] <- 150
s.sm$BDRICM[which(s.sm$site_key==16292)] <- 145
s.sm$BDRICM[which(s.sm$site_key==17420)] <- 16


## final merge:
profs.f <- join(h.s[,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT", "SLTPPT", "CRFVOL","PHIHO5","PHIKCL","BLD","ORCDRC")], s.sm[,c("SOURCEID", "LATWGS84", "LONWGS84", "TIMESTRR", "TAXNUSDA", "BDRICM")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
## check if there are still some duplicates
selD = profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
summary(selD)
#profs.f@site <- profs.f@site[!(profs.f@site$SOURCEID %in% selD),]
#profs.f@horizons <- profs.f@horizons[!(profs.f@horizons$SOURCEID %in% selD),]
sptID <- as.factor(paste(profs.f@horizons$SOURCEID, profs.f@horizons$TAXNUSDA, profs.f@horizons$TIMESTRR, profs.f@horizons$BDRICM, sep="__"))  ## round(profs.f@horizons$LONWGS84, 4), round(profs.f@horizons$LATWGS84, 4),
selS <- sptID %in% unique(sptID)
## extract site data
## THE FOLLOWING OPERATION TAKES >10 mins and about 6GB RAM!!!
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TAXNUSDA + TIMESTRR + BDRICM
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
str(sp::zerodist(sp))
## 3072 duplicate points!

profs.f@site$SOURCEDB = "NCSS"
## prepare a SITE and hor table table:
#wsp7 <- list(sites=profs.f@site, horizons=profs.f@horizons)
wsp7 <- list(sites=s.sm[,c("SOURCEID", "LATWGS84", "LONWGS84", "TIMESTRR", "TAXNUSDA", "BDRICM")], horizons=h.s[,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT", "SLTPPT", "CRFVOL","PHIHO5","PHIKCL","BLD","ORCDRC")])
str(wsp7)
lapply(as.list(wsp7$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp7$horizons), function(x){sum(!is.na(x))})
wsp7$sites$TAXNUSDA <- as.character(wsp7$sites$TAXNUSDA)
wsp7$sites$SOURCEDB = "NCSS"
save(wsp7, file="../wsp7.rda", compress="xz")




save.image("Repository2010.RData")

# end of script;