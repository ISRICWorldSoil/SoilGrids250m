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
library(stringr)

## list of input data sets:
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("TAXNWRB.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst ## 14
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})

###################################
## Pseudo-points
###################################

## Desert soils from African soil atlas:
afsoil <- readOGR("G:\\soilstorage\\PolygonMaps\\Africa\\afticasoilmap.shp", "afticasoilmap")
levels(as.factor(afsoil$SU_WRB1_PH))
afsoil <- afsoil[!is.na(afsoil$SU_WRB1_PH),]
## GH: Make a regular sample with size proportional to area
## TH: because the areas are patchy, this is difficult to implement as the spsample kicks out many points - random sample seems to be more suited;
Protic_Arenosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="ARpr",], type="random", n=300)
Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LP",], type="random", n=200)
Lithic_Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LPli",], type="random", n=100)
Rendzic_Leptosol.sim <- spsample(afsoil[afsoil$SU_WRB1_PH=="LPrz",], type="random", n=100)

## Indonesia:
load("../TAXOUSDA/Indonesia.ov.rda")
cleanup_Indonesia <- read.csv("cleanup_Indonesia.csv")
Indonesia.ov <- join(Indonesia.ov, cleanup_Indonesia, type="left")

## Glacier areas from (http://www.glims.org/download/)
## Circum-polar_regions points (http://bolin.su.se/data/ncscd/)
x = lapply(paste0("../TAXOUSDA/", c("turbels.pnt.rda", "orthels.pnt.rda", "histels.pnt.rda", "gelands.pnt.rda", "gelods.pnt.rda", "turbels2.pnt.rda", "orthels2.pnt.rda")), load, .GlobalEnv)

n.lst <- c("Protic Arenosols", "Leptosols", "Lithic Leptosols", "Rendzic Leptosols", "Haplic Cryosols", "Haplic Cryosols", "Cryic Histosols", "Haplic Cryosols", "Haplic Cryosols")
TAXNWRB.sim <- lapply(list(Protic_Arenosol.sim, Leptosol.sim, Lithic_Leptosol.sim, Rendzic_Leptosol.sim, turbels.pnt, orthels.pnt, histels.pnt, turbels2.pnt, orthels2.pnt), spTransform, CRS("+proj=longlat +datum=WGS84"))
TAXNWRB.sim.df <- list(NULL)
for(j in 1:length(TAXNWRB.sim)){
   TAXNWRB.sim.df[[j]] <- cbind(as.data.frame(TAXNWRB.sim[[j]]), TAXNWRB=rep(n.lst[j], length(TAXNWRB.sim[[j]])))
}
TAXNWRB.sim.df[[length(n.lst)+1]] <- as.data.frame(Indonesia.ov[,c("TAXNWRB","x","y")])
TAXNWRB.sim <- do.call(rbind, TAXNWRB.sim.df)
str(TAXNWRB.sim) ## 2458
TAXNWRB.sim <- plyr::rename(TAXNWRB.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXNWRB.sim$SOURCEID <- paste("SIM", 1:nrow(TAXNWRB.sim), sep="_")
TAXNWRB.sim$SOURCEDB = "Simulated"
in.lst[[length(tax.lst)+1]] <- TAXNWRB.sim

###################################
## Bind together and clean up names
###################################

all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
all.pnts <- as.data.frame(all.pnts)
coordinates(all.pnts) <- ~ LONWGS84+LATWGS84
proj4string(all.pnts) <- "+proj=longlat +datum=WGS84"
## Remove spatial duplicates:
all.pnts$LOC_ID <- as.factor(paste(all.pnts@coords[,1], all.pnts@coords[,2], sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## >2000 duplicate points
#data.frame(all.pnts[all.pnts$LOC_ID=="25.944444_-24.561111",])
#all.pnts <- all.pnts[!duplicated(all.pnts$LOC_ID),]
## TH: We keep the duplicates because there are indeed many points with approximate geocoordinates
length(all.pnts)
## 52742

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

## Search by pattern:
x <- summary(as.factor(all.pnts$TAXNWRB), maxsum = 1000)
soiltype <- data.frame(TAXNWRB=attr(x, "names"), count=x)
soiltype$WRB_2006_NAME <- NA
for(j in 1:length(levs)){
  pat1 <- ifelse(substr(levs[j], nchar(levs[j]), nchar(levs[j]))=="s", substr(levs[j], 1, nchar(levs[j])-1), levs[j])
  sel1 <- grep(pat1, paste(soiltype$TAXNWRB), ignore.case=TRUE)
  lev2 <- as.character(levs.c.1974$FAO_1974_NAME[which(levs.c.1974$WRB_2006_NAME %in% levs[j])])
  if(any(!is.na(lev2))){
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
soiltype$WRB_2006_NAMEf <- as.character(plyr::join(soiltype[,c("TAXNWRB","count")], legFAO, type="left", by="TAXNWRB", match="first")$WRB_2006_NAME)
#soiltype$WRB_2006_NAMEf <- ifelse(is.na(soiltype$WRB_2006_NAMEf), soiltype$WRB_2006_NAME, soiltype$WRB_2006_NAMEf)
all.pnts$WRB_2006_NAMEf <- plyr::join(all.pnts@data["TAXNWRB"], soiltype, type="left", by="TAXNWRB")$WRB_2006_NAMEf
soiltype.leg <- aggregate(soiltype[,c("count")], by=list(soiltype$WRB_2006_NAMEf), FUN=sum)
str(soiltype.leg)
names(soiltype.leg) = c("WRB_2006_NAMEf","count")
#soiltype.leg <- soiltype.leg[rank(1/soiltype.leg$count)<100,]
soiltype.leg$Group <- sapply(as.character(soiltype.leg$WRB_2006_NAMEf), function(x){strsplit(x, " ")[[1]][2]})
write.csv(soiltype.leg, file="soiltype_legend.csv")

###################################
## THE FINAL LEGEND
###################################

## Selection of soil types for SoilGrids250m
levsf0 <- levels(as.factor(soiltype.leg$WRB_2006_NAMEf))
levsf <- c(levsf0, "Sapric Histosols", "Hemic Histosols", "Alic Nitisols", "Haplic Albeluvisols", "Cutanic Alisols", "Petric Durisols", "Cryic Histosols", "Leptic Umbrisols", "Haplic Umbrisols", "Luvic Stagnosols", "Lixic Plinthosols", "Vitric Cryosols", "Histic Albeluvisols")
levsf <- unique(str_trim(as.character(unlist(sapply(levsf, function(x){strsplit(x, "/")[[1]]})))))
#levsf <- levsf[-c(grep(pattern="Anthrosol", levsf), grep(pattern="Technic", levsf))]
## remove Anthrosols because there are still too difficult to map in automated manner
str(levsf)
soiltype.leg2 <- data.frame(WRB_2006_NAMEf=levsf)
soiltype.leg2$Group <- sapply(as.character(soiltype.leg2$WRB_2006_NAMEf), function(x){strsplit(x, " ")[[1]][2]})
write.csv(soiltype.leg2, file="TAXNWRB_legend2.csv")
## 118 classes on the end

## One by one subgroup name --> try to located them in the raw names
tax.lst <- list(NULL)
for(j in 1:length(levsf)){
  ## remove "s" if at the end of the class name:
  pat <- gsub("chaks", "chak", gsub("zems", "zem", gsub("sols", "sol", levsf[j])))  
  sel1 <- grep(pat, all.pnts$WRB_2006_NAMEf, ignore.case=TRUE)
  sel2 <- which(all.pnts$WRB_2006_NAMEf==levsf[j])
  sel3 <- grep(pat, all.pnts$TAXNWRB, ignore.case=TRUE)
  sel4 <- which(all.pnts$TAXNWRB==levsf[j])
  sel <- unique(c(sel1, sel2, sel3, sel4))
  ## there can be multiple soil taxa at the same location
  if(length(sel)>0){
    tax.lst[[j]] <- data.frame(all.pnts[sel,])
    tax.lst[[j]]$TAXNWRB.f <- levsf[j]
  }
}
TAXNWRB.pnts <- do.call(rbind, tax.lst)
TAXNWRB.pnts$TAXNWRB.f <- as.factor(TAXNWRB.pnts$TAXNWRB.f)
TAXNWRB.pnts <- TAXNWRB.pnts[,]
summary(TAXNWRB.pnts$TAXNWRB.f)
xS <- summary(as.factor(TAXNWRB.pnts$TAXNWRB.f), maxsum = 118)
soiltype_count <- data.frame(Group=attr(xS, "names"), count=xS)
write.csv(soiltype_count, file="soiltype_count2.csv")
length(TAXNWRB.pnts$TAXNWRB.f)
## FINAL NUMBER OF POINTS: 64,735 points

###################################
## Export points
###################################

summary(as.factor(TAXNWRB.pnts$SOURCEDB))
#          AfSPDB            eSOTER       RadamBrasil              ISIS 
#             3184              2132              6448               809 
#             WISE          Can_FECD            CanSIS          CN-SOTER 
#            10831              1878              8765              1626 
#           HRSPDB           MX_CDPS      Russia_EGRPR Alterra-BODEMDATA 
#             2854              5018               869               402 
#       ISCN2012/N        HILATS2014            IRSPDB          NAMSOTER 
#             9025               334              5223              2431 
#             OFRA             SPADE             CIFOR         Simulated 
#               98               134               207              2467

coordinates(TAXNWRB.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.pnts) <- "+proj=longlat +datum=WGS84"
TAXNWRB.pnts <- TAXNWRB.pnts[!(TAXNWRB.pnts@coords[,2]==0&TAXNWRB.pnts@coords[,1]==0),] ## c("SOURCEDB","SOURCEID","TAXNWRB.f")
unlink("TAXNWRB.pnts.shp")
writeOGR(TAXNWRB.pnts, "TAXNWRB.pnts.shp", "TAXNWRB.pnts", "ESRI Shapefile")
#plotKML(TAXNWRB.pnts["TAXNWRB.f"], file.name="TAXNWRB_Feb_15_2016.kml", kmz=TRUE)
str(TAXNWRB.pnts@data)
save(TAXNWRB.pnts, file="TAXNWRB.pnts.rda")
