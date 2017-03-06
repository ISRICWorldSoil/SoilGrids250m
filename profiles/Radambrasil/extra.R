# ------------------------------------------------------------
# another copy of the same thing...
# ------------------------------------------------------------

##The National Soil Profile Database for Brazil Available to International Scientists [https://www.soils.org/publications/sssaj/articles/69/3/0649], contact: Miguel Cooper;
## Download the database:
#download.file("http://www.esalq.usp.br/gerd/BrazilSoilDB_08VI05.xls", destfile="BrazilSoilDB_08VI05.xls")

library(aqp)
library(plyr)
library(GSIF)
library(maptools)
library(gdata)
#perl <- gdata:::findPerl("perl")
#perl = "C:/Perl64/bin/perl.exe"
#tmp <- read.xls("BrazilSoilDB_08VI05.xls", perl=perl, sheet=2)
## takes 2-3 mins to read!
#profs <- tmp
# str(profs)

profs <- read.csv("BrazilSoilDB_08VI05.csv")
## ID column:
profs$SOURCEID <- as.factor(paste(iconv(profs$Source, to="ASCII", sub="byte"), profs$OrgProfID, sep="_"))
# check if there are duplicates:
length(levels(profs$SOURCEID))
## 5858 profile IDs
plyr:::nunique(paste(profs$Long, profs$Lat, sep="_"))
## 5772 points

## remove "_" sign from the end:
profs$SoilClass <- gsub(pattern="*_$", replacement="", profs$SoilClass)
profs$taxon_BR <- as.factor(gsub("_", " ", profs$SoilClass))
## WRB classes added by Alessandro Rosa / Eliana de Souza:
taxa <- read.xls("BR_taxa_cleanup.xls", perl=perl, sheet=8)
profs <- merge(profs, taxa[,c("taxon_BR","taxon_WRB")], sort=FALSE, all.x=TRUE)
profs$taxon_WRB[which(profs$taxon_WRB=="#N/A!")] = NA
profs$taxon_WRB <- as.factor(paste(profs$taxon_WRB))
summary(profs$taxon_WRB)
summary(profs$SoilClass=="")
## Estimate "DBRICM" -> observed soil depth where word "Litólico" is in "SoilClass" and/or "HzSimb" has letter "R/C";
profs$SoilDepth_R <- NA
sel.d <- grep("Lit", profs$taxon_BR, ignore.case=TRUE, fixed=FALSE)
sel.r <- grep("R", profs$HzSimb, ignore.case=FALSE, fixed=FALSE)
profs$SoilDepth_R[sel.d] <- profs$SoilDepth[sel.d]
profs$SoilDepth_R[sel.r] <- profs$SoilDepth[sel.r]

# rename columns:
profs.f <- rename(profs, c("taxon_WRB"="TAXGWRB", "Long"="LONWGS84", "Lat"="LATWGS84", "HzDeIn"="UHDICM", "HzDeFn"="LHDICM", "pH_H2O"="PHIHO5", "pH_KCL"="PHIKCL", "C"="ORCDRC", "Sand"="SNDPPT", "Silt"="SLTPPT", "Clay"="CLYPPT", "CEC_pH7"="CEC", "SoilDepth_R"="BDRICM"))
summary(profs.f$LHDICM)
summary(profs.f$PHIHO5)
profs.f$PHIHO5 <- ifelse(profs.f$PHIHO5 < 2 | profs.f$PHIHO5 > 10, NA, profs.f$PHIHO5)
profs.f$PHIKCL <- ifelse(profs.f$PHIKCL < 2 | profs.f$PHIKCL > 10, NA, profs.f$PHIKCL)
summary(profs.f$ORCDRC)
# Organic carbon in permilles:
profs.f$ORCDRC <- profs.f$ORCDRC*10
profs.f$CRFVOL <- (profs.f$CG + as.integer(profs.f$FG))/2
# subset to complete data:
sel.c <- !is.na(profs.f$LHDICM)&!is.na(profs.f$UHDICM)&!is.na(profs.f$LONWGS84)&!is.na(profs.f$LATWGS84)
summary(sel.c)
profs.spc <- profs.f[sel.c,]

## convert to SPC class:
depths(profs.spc) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.spc) <- ~ LONWGS84 + LATWGS84
## [1] "pedons (5820) rows of site data (5859)"
## TH: there are about 30-40 profiles with the same ID but different coordinates
plyr:::nunique(profs.spc@horizons$SOURCEID)
plyr:::nunique(paste(profs.spc@horizons$LONWGS84, profs.spc@horizons$LATWGS84, sep="_"))
## requires a check...
sel <- sapply(profs.spc@site$SOURCEID, function(x){sd(profs.spc@horizons[profs.spc@horizons$SOURCEID==x,"LONWGS84"])>0|sd(profs.spc@horizons[profs.spc@horizons$SOURCEID==x,"LATWGS84"])>0|length(unique(profs.spc@horizons[profs.spc@horizons$SOURCEID==x,"TAXGWRB"]))>1})
## takes few minutes...
summary(sel)
sel <- sel[!is.na(sel)]
selp <- names(sel)[sel]
profs.spc <- profs.f[sel.c & !(profs.f$SOURCEID %in% selp),]
depths(profs.spc) <- SOURCEID ~ UHDICM + LHDICM
site(profs.spc) <- ~ LONWGS84 + LATWGS84 + TAXGWRB
## spatial duplicates:
sp <- profs.spc@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
seld <- sp::zerodist(sp)
str(seld) ## 143 spatial duplicates!

## add missing columns:
profs.spc@site$SOURCEDB <- "BrazilSoilDB"

# prepare the dataset:
wsp2 <- list(sites=profs.spc@site, horizons=profs.spc@horizons[,c("SOURCEID","UHDICM","LHDICM","CRFVOL","PHIHO5","PHIKCL","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC","BDRICM")])
wsp2$sites$TAXGWRB <- as.character(wsp2$sites$TAXGWRB)
str(wsp2)
lapply(as.list(wsp2$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp2$horizons), function(x){sum(!is.na(x))})
save(wsp2, file="../wsp2.rda", compress="xz")

## export sites only:
profs.csv <- profs.f[sel.c & !(profs.f$SOURCEID %in% sel),]
str(profs.csv)
## clean up non-unique entries for taxonomy:
for(j in levels(profs.csv$SOURCEID)){
  x <- profs.csv[profs.csv$SOURCEID==j,c("TAXGWRB")]
  if(!summary(x)[[1]]==length(x)){
    profs.csv[profs.csv$SOURCEID==j,c("TAXGWRB")] <- names(summary(x)[1])
  }
} ## takes 3-4 minutes!
depths(profs.csv) <- SOURCEID ~ UHDICM + LHDICM
site(profs.csv) <- ~ LONWGS84 + LATWGS84 + TAXGWRB
coordinates(profs.csv) <- ~ LONWGS84 + LATWGS84
proj4string(profs.csv@sp) <- "+proj=longlat +datum=WGS84"
profs.csv@horizons <- profs.csv@horizons[,c("SOURCEID", "SoilClass", "UHDICM", "LHDICM", "HzSimb", "ColorMunsell", "SLTPPT", "CLYPPT", "SNDPPT", "PHIHO5", "PHIKCL", "ORCDRC", "CEC_pH7", "BSum", "BSat", "ALSat", "BDRICM")]
library(GSIF)
x <- as.data.frame(profs.csv)
str(x)
write.csv(x, "BrazilSoilDB_08VI05_filtered.csv")

## end of script;