## Mexican soil profile DB from CONABIO (5781 profiles);
## Instituto Nacional de Estad?stica y Geograf?a (INEGI), 2000, Conjunto de Datos de Perfiles de Suelos, Escala 1: 250 000 Serie II (Continuo Nacional). [http://www.inegi.org.mx/geo/contenidos/recnat/edafologia/PerfilesSuelo.aspx];
## Prepared by Mario A. Guevara (mguevara@conabio.gob.mx) and T. Hengl
## Mexican coordinate system (http://spatialreference.org/ref/sr-org/39/proj4/); point data are from 1980-2007.

load(".RData")
library(aqp)
library(GSIF)
library(plotKML)
library(rgdal)
library(sp)
library(maptools)
library(plyr)
mx.csy <- "+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs"

## Conjunto de Datos de Perfiles de Suelos, Escala 1:250 000 Serie II (Continuo Nacional)
edaf <- readOGR("perf_edaf_sii/edaf_puntos_sii.shp", "edaf_puntos_sii")
str(edaf, max.level=2)
proj4string(edaf) = mx.csy
edaf <- spTransform(edaf, CRS("+proj=longlat +datum=WGS84"))
## 16,820 points!
edaf$SOURCEID <- substr(edaf$ID_PERFIL, 1, 6)
write.csv(as.data.frame(edaf), "edaf_puntos_sii.csv")
SITE <- data.frame(edaf[!duplicated(edaf$SOURCEID),c("SOURCEID","ID_PERFIL","FECHA","GPO_SUELO","CALIF_PRIM")])
s <- summary(SITE$CALIF_PRIM)
soiltype <- data.frame(CALIF_PRIM=attr(s, "names"), count=s)
write.csv(soiltype, "soiltype_count.csv")
SITE$SOURCEDB = "MX_CDPS"
legFAO_90 <- read.csv("cleanup_MX_FAO.csv", fileEncoding="UTF-8")
SITE$TAXNWRB <- paste(join(SITE, legFAO_90, type="left")$WRB_2nd, SITE$GPO_SUELO)
SITE$TIMESTRR <- as.Date(paste(SITE$FECHA), format="%d/%m/%Y")
SITE <- rename(SITE, c("coords.x1"="LONWGS84", "coords.x2"="LATWGS84"))
#View(SITE)

perfilv12 <- readOGR("perf_edaf_si/perfilv12.shp", "perfilv12")
proj4string(perfilv12) = mx.csy
perfilv12 <- spTransform(perfilv12, CRS("+proj=longlat +datum=WGS84"))
perfilv12$SOURCEID <- paste(perfilv12$IDENTIFI, perfilv12$CLAVE_250, sep="_")
#SITE <- data.frame(perfilv12[!duplicated(perfilv12$SOURCEID),c("SOURCEID",)])
#write.csv(as.data.frame(perfilv12), "perfilv12.csv")

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.MX_CDPS <- SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB","TIMESTRR")]
TAXNWRB.MX_CDPS <- TAXNWRB.MX_CDPS[!is.na(TAXNWRB.MX_CDPS$TAXNWRB)&!is.na(TAXNWRB.MX_CDPS$LONWGS84)&nchar(paste(TAXNWRB.MX_CDPS$TAXNWRB))>0,]
str(TAXNWRB.MX_CDPS)
## 4,428 profiles
coordinates(TAXNWRB.MX_CDPS) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.MX_CDPS) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.MX_CDPS["TAXNWRB"])
save(TAXNWRB.MX_CDPS, file="TAXNWRB.MX_CDPS.rda")


# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons <- as.data.frame(edaf)[,c("ID_PERFIL","LIM_SUP","LIM_INF","R","L","A","PH","CO","CIC","coords.x1","coords.x2","SOURCEID","FECHA")]
horizons <- rename(horizons, c("ID_PERFIL"="SAMPLEID","LIM_SUP"="UHDICM","LIM_INF"="LHDICM","R"="SNDPPT","L"="SLTPPT","A"="CLYPPT","PH"="PHIHOX","CO"="ORCDRC","CIC"="CECSUM","coords.x1"="LONWGS84","coords.x2"="LATWGS84"))
horizons$ORCDRC <- horizons$ORCDRC*10
summary(horizons$ORCDRC)
horizons$TIMESTRR <- as.Date(paste(horizons$FECHA), format="%d/%m/%Y")
summary(horizons$TIMESTRR)
## filter out all zeros!
horizons <- horizons[!horizons$PHIHOX==0&!horizons$SNDPPT==0,]
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
## 13,939
str(horizons)
horizons$SOURCEDB = "MX_edaf_puntos"

horizons2 <- as.data.frame(perfilv12)[,c("NHORIZON","HLIMSUPE","HLIMINFE","ARCILLA","LIMO","ARENA","PH","MO","CIC","coords.x1","coords.x2","SOURCEID")]
horizons2 <- rename(horizons2, c("HLIMSUPE"="UHDICM","HLIMINFE"="LHDICM","ARCILLA"="CLYPPT","LIMO"="SLTPPT","ARENA"="SNDPPT","PH"="PHIHOX","MO"="ORCDRC","CIC"="CECSUM","coords.x1"="LONWGS84","coords.x2"="LATWGS84"))
horizons2$ORCDRC <- horizons2$ORCDRC*10/1.724
summary(horizons2$ORCDRC)
## filter out all zeros!
horizons2$DEPTH <- horizons2$UHDICM + (horizons2$LHDICM - horizons2$UHDICM)/2
horizons2 <- horizons2[!is.na(horizons2$DEPTH),]
## 10,890
horizons2$SAMPLEID <- make.unique(paste(horizons2$SOURCEID, horizons2$NHORIZON, sep="_"))
horizons2$SOURCEDB = "MX_perfilv12"
horizons2$TIMESTRR = NA
hist(horizons2$CLYPPT, col="gray")
## points from 'edaf' have a higher clay content
hist(horizons$CLYPPT, col="gray")
hist(horizons$SNDPPT, col="gray")

sel.n <- c("SOURCEID","SOURCEDB","TIMESTRR","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","PHIHOX","ORCDRC","CECSUM","LONWGS84","LATWGS84")
SPROPS.MX_PdS <- rbind(horizons[,sel.n], horizons2[,sel.n])
str(SPROPS.MX_PdS)
## 24,829
save(SPROPS.MX_PdS, file="SPROPS.MX_PdS.rda")
plot(SPROPS.MX_PdS$LONWGS84, SPROPS.MX_PdS$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

## Depth to bedrock i.e. 'R' horizon:
sel.r1 <- grep(pattern="x", perfilv12$LIM_ROCA, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="x", perfilv12$LIM_REGO, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="^R", perfilv12$HSIMBOLO, ignore.case=FALSE, fixed=FALSE)
perfilv12$BDRICM <- 250
perfilv12$BDRICM[unique(c(sel.r1, sel.r2, sel.r3))] <- perfilv12$PROFUNDI[unique(c(sel.r1, sel.r2, sel.r3))]
summary(perfilv12$BDRICM<250)
BDR.perfilv12 <- as.data.frame(perfilv12[c("SOURCEID","BDRICM")])
BDR.perfilv12 <- rename(BDR.perfilv12, c("coords.x1"="LONWGS84", "coords.x2"="LATWGS84"))
BDR.perfilv12$SOURCEDB <- "MX_INEGI"

summary(edaf$CALIF_PRIM) ## "Lítico", "squelético"
summary(edaf$GPO_SUELO)  ## LEPTOSOL, REGOSOL
summary(edaf$NOMEN_HTE)  ## CR?, R
sel2.r1 <- grep(pattern="Lítico", edaf$CALIF_PRIM, ignore.case=FALSE, fixed=FALSE)
sel2.r2 <- grep(pattern="LEPT", edaf$GPO_SUELO, ignore.case=FALSE, fixed=FALSE)
sel2.r3 <- grep(pattern="R", edaf$NOMEN_HTE, ignore.case=FALSE, fixed=FALSE)
sel.t <- unique(c(sel2.r1, sel2.r2, sel2.r3))
edaf$BDRICM <- NA
edaf[sel.t,"BDRICM"] <- pmax(edaf$LIM_SUP[sel.t], edaf$LIM_INF[sel.t])
bdr.d <- aggregate(edaf$BDRICM, list(edaf$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.MX_PdS <- join(bdr.d, SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
BDR.MX_PdS$BDRICM <- ifelse(is.infinite(BDR.MX_PdS$BDRICM), 250, BDR.MX_PdS$BDRICM)
BDR.MX_PdS <- rbind.fill(BDR.MX_PdS, BDR.perfilv12)
str(BDR.MX_PdS)
summary(BDR.MX_PdS$BDRICM<250)
## 7906 points
save(BDR.MX_PdS, file="BDR.MX_PdS.rda")