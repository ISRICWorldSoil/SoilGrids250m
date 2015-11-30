## Prepare callibration points for TAXNWRB (SoilGrids250m)
## By Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(GSIF)
library(rgdal)
library(lattice)
library(scales)
library(plotKML)
library(maps)
library(stringr)

## list of input data sets:
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("TAXNWRB.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst ## 14
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})

## Desert soils from African soil atlas:
afsoil <- readOGR("G:\\soilstorage\\PolygonMaps\\Africa\\afticasoilmap.shp", "afticasoilmap")
levels(as.factor(afsoil$SU_WRB1_PH))
afsoil <- afsoil[!is.na(afsoil$SU_WRB1_PH),]
Protic_Arenosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="ARpr",], type="random", n=300)
Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LP",], type="random", n=200)
Lithic_Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LPli",], type="random", n=100)
Rendzic_Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LPrz",], type="random", n=100)

n.lst <- c("Protic Arenosols", "Leptosols", "Lithic Leptosols", "Rendzic Leptosols")
TAXNWRB.sim <- list(Protic_Arenosol.sim, Leptosol.sim, Lithic_Leptosol.sim, Rendzic_Leptosol.sim)
TAXNWRB.sim.df <- list(NULL)
for(j in 1:length(TAXNWRB.sim)){
   TAXNWRB.sim.df[[j]] <- cbind(as.data.frame(TAXNWRB.sim[[j]]), TAXNWRB=rep(n.lst[j], length(TAXNWRB.sim[[j]])))
}
TAXNWRB.sim <- do.call(rbind, TAXNWRB.sim.df)
str(TAXNWRB.sim)
TAXNWRB.sim <- plyr::rename(TAXNWRB.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXNWRB.sim$SOURCEID <- paste("SIM", 1:nrow(TAXNWRB.sim), sep="_")
TAXNWRB.sim$SOURCEDB = "Simulated"
in.lst[[length(in.lst)+1]] <- TAXNWRB.sim

## Bind everything together:
all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
## 46,488
all.pnts <- as.data.frame(all.pnts)
coordinates(all.pnts) <- ~ LONWGS84+LATWGS84
proj4string(all.pnts) <- "+proj=longlat +datum=WGS84"
## Remove spatial duplicates (except for the NCSS data):
all.pnts$LOC_ID <- as.factor(paste(all.pnts@coords[,1], all.pnts@coords[,2], sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## 2603 duplicate points
## data.frame(all.pnts[all.pnts$LOC_ID=="25.944444_-24.561111",])
#TAXNWRB.pnts <- all.pnts[!duplicated(all.pnts$LOC_ID),]
## TH: We keep the duplicates because there are indeed many points with approximate geocoordinates
length(all.pnts)

## clean-up names:
FAO_levs <- read.csv("WRB_versions.csv")
str(FAO_levs)
levs <- levels(FAO_levs$WRB_2006_NAME)
str(levs)
# 511 subgroups, TOO MANY for SoilGrids250m!
levels(as.factor(all.pnts$SOURCEDB))
## Correlation tables:
levs.c.1974 <- join(data.frame(WRB_2006_NAME=levs), FAO_levs[,c("WRB_2006_NAME","FAO_1974_NAME")], type="left", match="first")
w2 <- sapply(gregexpr("[[:alpha:]]+", levs.c.1974$FAO_1974_NAME), function(x) sum(x > 0))
levs.c.1974 <- levs.c.1974[!(as.character(levs.c.1974$WRB_2006_NAME)==as.character(levs.c.1974$FAO_1974_NAME))&w2==2,]
write.csv(levs.c.1974, file="FAO1974_to_WRB2006.csv")
levs.c.1994 <- join(data.frame(WRB_2006_NAME=levs), FAO_levs[,c("WRB_2006_NAME","WRB_1994_NAME")], type="left", match="first")
w2t <- sapply(gregexpr("[[:alpha:]]+", levs.c.1994$WRB_1994_NAME), function(x) sum(x > 0))
levs.c.1994 <- levs.c.1994[!(as.character(levs.c.1994$WRB_2006_NAME)==as.character(levs.c.1994$WRB_1994_NAME))&w2t==2,]
write.csv(levs.c.1994, file="WRB1994_to_WRB2006.csv")

x <- summary(as.factor(all.pnts$TAXNWRB), maxsum = 1000)
soiltype <- data.frame(TAXNWRB=attr(x, "names"), count=x)
soiltype$WRB_2006_NAME <- NA
for(j in 1:length(levs)){
  pat1 <- ifelse(substr(levs[j], nchar(levs[j]), nchar(levs[j]))=="s", substr(levs[j], 1, nchar(levs[j])-1), levs[j])
  sel1 <- grep(pat1, paste(soiltype$TAXNWRB), ignore.case=TRUE)
  lev2 <- as.character(levs.c.1974$FAO_1974_NAME[match(levs[j], FAO_levs$WRB_2006_NAME)])
  if(!is.na(lev2)){
    pat2 <- ifelse(substr(lev2, nchar(lev2), nchar(lev2))=="s", substr(lev2, 1, nchar(lev2)-1), lev2)
    sel2 <- grep(pat2, paste(soiltype$TAXNWRB), ignore.case=TRUE)
  } else {
    sel2 <- NULL
  }
  soiltype[c(sel1, sel2),"WRB_2006_NAME"] = levs[j]
}
write.csv(soiltype, "soiltype_count.csv")

## Import the clean-up tables (clean-up by Maria and Tom)
legFAO <- read.csv("TAXNWRB_cleanup.csv", fileEncoding="UTF-8")
soiltype$WRB_2006_NAMEf <- as.character(join(soiltype[,c("TAXNWRB","count")], legFAO, type="left", by="TAXNWRB")$WRB_2006_NAME)
soiltype$WRB_2006_NAMEf <- ifelse(is.na(soiltype$WRB_2006_NAMEf), soiltype$WRB_2006_NAME, soiltype$WRB_2006_NAMEf)
all.pnts$WRB_2006_NAMEf <- join(all.pnts@data["TAXNWRB"], soiltype, type="left", by="TAXNWRB")$WRB_2006_NAMEf
soiltype.leg <- aggregate(soiltype[,c("count")], by=list(soiltype$WRB_2006_NAMEf), FUN=sum)
str(soiltype.leg)
names(soiltype.leg) = c("WRB_2006_NAMEf","count")
soiltype.leg <- soiltype.leg[rank(1/soiltype.leg$count)<100,]
soiltype.leg$Group <- sapply(as.character(soiltype.leg$WRB_2006_NAMEf), function(x){strsplit(x, " ")[[1]][2]})
write.csv(soiltype.leg, file="soiltype_legend.csv")

## Selection of soil types for SoilGrids250m
## THE FINAL LEGEND
levsf0 <- levels(as.factor(soiltype.leg$WRB_2006_NAMEf))
levsf <- c(levsf0, "Sapric Histosols", "Hemic Histosols", "Alic Nitisols", "Haplic Albeluvisols", "Cutanic Alisols", "Regic Anthrosols", "Petric Durisols", "Cryic Histosols", "Leptic Umbrisols", "Haplic Umbrisols", "Luvic Stagnosols", "Lixic Plinthosols", "Vitric Cryosols", "Histic Albeluvisols")
levsf <- unique(str_trim(as.character(unlist(sapply(levsf, function(x){strsplit(x, "/")[[1]]})))))
levsf <- levsf[-c(grep(pattern="Anthrosol", levsf), grep(pattern="Technic", levsf))]
## remove Anthrosols because there are still too difficult to map in automated manner
str(levsf)
soiltype.leg2 <- data.frame(WRB_2006_NAMEf=levsf)
soiltype.leg2$Group <- sapply(as.character(soiltype.leg2$WRB_2006_NAMEf), function(x){strsplit(x, " ")[[1]][2]})
write.csv(soiltype.leg2, file="TAXNWRB_legend2.csv")
## 108 classes on the end

## One by one subgroup name --> try to located them in the raw names
tax.lst <- list(NULL)
for(j in 1:length(levsf)){
  ## remove "s" if at the end of the class name:
  pat <- gsub("chaks", "chak", gsub("zems", "zem", gsub("sols", "sol", levsf[j])))  
  sel1 <- which(paste(all.pnts$WRB_2006_NAMEf) %in% pat)
  sel2 <- which(paste(all.pnts$WRB_2006_NAMEf) %in% levsf[j])
  sel3 <- grep(pat, paste(all.pnts$TAXNWRB), ignore.case=TRUE)
  sel <- unique(c(sel1, sel2, sel3))
  ## there can be multiple soil taxa at the same location
  if(length(sel)>0){
    tax.lst[[j]] <- data.frame(all.pnts[sel,])
    tax.lst[[j]]$TAXNWRB.f <- levsf[j]
  }
}
TAXNWRB.pnts <- do.call(rbind, tax.lst)
TAXNWRB.pnts$TAXNWRB.f <- as.factor(TAXNWRB.pnts$TAXNWRB.f)
summary(TAXNWRB.pnts$TAXNWRB.f)
length(TAXNWRB.pnts$TAXNWRB.f)
## FINAL NUMBER OF POINTS: 42,518 points
summary(as.factor(TAXNWRB.pnts$SOURCEDB))
#      AfSPDB       eSOTER  RadamBrasil     CN-SOTER       HRSPDB       IRSPDB   ISCN2012/N         ISIS 
#        1863         1854         5359         1272         1812         4412         7830          530 
#     MX_CDPS     NAMSOTER        SPADE         WISE Russia_EGRPR         OFRA        CIFOR   HILATS2014 
#        2567         1807           82         7096          439           28          196          334 
#      CanSIS    Simulated 
#        5140          500 

coordinates(TAXNWRB.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.pnts) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.pnts["TAXNWRB.f"], file.name="TAXNWRB_Nov_4_2015.kml", kmz=TRUE)
str(TAXNWRB.pnts@data)
TAXNWRB.pnts <- TAXNWRB.pnts[c("SOURCEDB","SOURCEID","TAXNWRB.f")]
save(TAXNWRB.pnts, file="TAXNWRB.pnts.rda")

## world plot - overlay and plot points and maps:
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
no.plt <- TAXNWRB.pnts@coords[,2]>-65&TAXNWRB.pnts@coords[,2]<85
png(file = "Fig_global_distribution_TAXNWRB.png", res = 150, width = 2000, height = 900)
windows(width = 20, height = 9)
dev.off()
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(TAXNWRB.pnts[!TAXNWRB.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
points(TAXNWRB.pnts[TAXNWRB.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()