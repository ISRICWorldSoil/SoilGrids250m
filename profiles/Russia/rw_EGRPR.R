## Russian profile data obtained from: http://egrpr.esoil.ru/
## by: Tom.Hengl@isric.org

#require(gdata)
library(XLConnect)
library(plyr)
Sys.setlocale('LC_ALL', 'russian')
library(sp)
library(plotKML)

#ru <- read.csv("Russia_EGRPR_soil_data.csv")
wb <- loadWorkbook("soil_data.xlsx")
ru <- readWorksheet(wb, sheet="soil_data", header=TRUE)
str(ru)
sel <- which(sapply(ru, is.character))

f.lst <- list(NULL)
for(j in 1:length(sel)){
  f.lst[[j]] <- levels(as.factor(ru[,sel[j]]))
}
names(f.lst) = attr(sel, "names")
str(f.lst)
max_length <- max(sapply(f.lst,length))
f <- data.frame(sapply(f.lst, function(x){ c(x, rep(NA, max_length - length(x))) }))
createSheet(wb, name = "factors")
writeWorksheet(wb, f, sheet = "factors")
saveWorkbook(wb)

## Translated levels:
tr <- read.csv("factors_translated.csv")
#ru.en <- join(ru, tr)

ru$SOURCEID <- as.factor(paste(ru$CardID, ru$SOIL_ID, sep="_"))
ru.xy <- ru[!duplicated(ru$SOURCEID),c("SOURCEID", "WRB06", "LONG", "LAT")]
nrow(ru.xy)
#coordinates(ru.xy) <- ~LONG+LAT
#proj4string(ru.xy) = "+proj=longlat +datum=WGS84"
#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#kml(ru.xy, file="profiles_Russia_EGRPR.kml", shape=shape, colour=WRB06, size=1, labels=ru.xy$WRB06, kmz=TRUE, balloon=TRUE)
#write.csv(as.data.frame(ru.xy), file="ru_profiles_ll.csv")

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.EGRPR <- rename(ru.xy, replace=c("LONG"="LONWGS84","LAT"="LATWGS84","WRB06"="TAXNWRB"))
TAXNWRB.EGRPR$SOURCEDB = "Russia_EGRPR"
TAXNWRB.EGRPR <- TAXNWRB.EGRPR[!is.na(TAXNWRB.EGRPR$TAXNWRB)&!is.na(TAXNWRB.EGRPR$LONWGS84)&nchar(paste(TAXNWRB.EGRPR$TAXNWRB))>0,]
str(TAXNWRB.EGRPR)
## 862 profiles
coordinates(TAXNWRB.EGRPR) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.EGRPR) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.EGRPR["TAXNWRB"])
save(TAXNWRB.EGRPR, file="TAXNWRB.EGRPR.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons <- ru[,c("SOURCEID","HIAUTH","HORTOP","HORBOT","ORGMAT","CECST","PHH2O","PHSLT","TEXTSAF","TEXTSIC","TEXSCM","TEXTSIM","TEXTSIF","TEXTCL","TEXTPHC","DVOL","GRVDEG")]
horizons$SAMPLEID <- make.unique(paste(horizons$SOURCEID, horizons$HIAUTH, sep="_"))
horizons <- rename(horizons, c("HORTOP"="UHDICM","HORBOT"="LHDICM","PHH2O"="PHIHOX","PHSLT"="PHIKCL","ORGMAT"="ORCDRC","CECST"="CECSUM","DVOL"="BLD"))
horizons$ORCDRC <- horizons$ORCDRC*10/1.724
summary(horizons$ORCDRC)
horizons$BLD <- horizons$BLD * 1000
levels(as.factor(horizons$GRVDEG))
horizons$CRFVOL <- 0
horizons$CRFVOL[horizons$GRVDEG=="???"] <- 0
horizons$CRFVOL[horizons$GRVDEG=="????????????"] <- 0.5
horizons$CRFVOL[horizons$GRVDEG=="???????????????"] <- (0.5+5)/2
horizons$CRFVOL[horizons$GRVDEG=="????????????????"] <- (5+10)/2
horizons$CRFVOL[horizons$GRVDEG=="????????????????"] <- 15
horizons$UHDICM <- as.numeric(horizons$UHDICM)
summary(horizons$UHDICM)
horizons$LHDICM <- as.numeric(horizons$LHDICM)
## filter out all zeros!
summary(horizons$PHIHOX)
summary(horizons$PHIKCL)
## convert texture fractions to USDA system:
horizons$SNDPPT <- horizons$TEXTSAF + horizons$TEXSCM
horizons$SLTPPT <- horizons$TEXTSIC + horizons$TEXTSIM + 0.8 * horizons$TEXTSIF
horizons$CLYPPT <- horizons$TEXTCL + 0.2 * horizons$TEXTSIF
## Correct texture fractions:
sumTex <- rowSums(horizons[,c("SLTPPT","CLYPPT","SNDPPT")])
horizons$SNDPPT <- horizons$SNDPPT / ((sumTex - horizons$CLYPPT) /(100 - horizons$CLYPPT))
horizons$SLTPPT <- horizons$SLTPPT / ((sumTex - horizons$CLYPPT) /(100 - horizons$CLYPPT))
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
nrow(horizons)
## 4961

SPROPS.EGRPR <- join(horizons[!is.na(horizons$DEPTH),c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","CRFVOL","BLD","SNDPPT","SLTPPT","PHIHOX","ORCDRC","CECSUM")], as.data.frame(TAXNWRB.EGRPR)[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
str(SPROPS.EGRPR)
## 11,232
save(SPROPS.EGRPR, file="SPROPS.EGRPR.rda")
plot(SPROPS.EGRPR$LONWGS84, SPROPS.EGRPR$LATWGS84, pch="+")
