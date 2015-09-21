## Prepare covariates for SoilGrids250m
## By: Tom.Hengl@isric.org

library(utils)
library(R.utils)
library(snowfall)
library(raster)
library(rgdal)
library(RSAGA)
library(plotKML)

gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
system('/usr/local/bin/saga_cmd --version')
saga_cmd <- "/usr/local/bin/saga_cmd"
load("equi7t3.rda")
source("tiler.R")
source("prepareCovs.R")
source("plotCovsSoilGrids250m.R")
source("make.covsRDA.R")

## Months of interest:
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
ms.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
msk.lst <- list.files(pattern=glob2rx("LMK_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
## 2594 tiles

## Data locations:
DEM.out.lst <- c("DEMMRG5", "SLPMRG5", "CRVMRG5", "VBFMRG5", "DVMMRG5", "VDPMRG5", "POSMRG5", "NEGMRG5", "TWIMRG5")
DEM.lst <- paste0("/data/MDEM/MDEM_", names(equi7t3), "_250m.tif")
SLP.lst <- paste0("/data/MDEM/MSLP_", names(equi7t3), "_250m.sdat")
CRV.lst <- paste0("/data/MDEM/MCRV_", names(equi7t3), "_250m.sdat")
DVM.lst <- paste0("/data/MDEM/MDVM_", names(equi7t3), "_250m.sdat")
NEG.lst <- paste0("/data/MDEM/MNEG_", names(equi7t3), "_250m.sdat")
POS.lst <- paste0("/data/MDEM/MPOS_", names(equi7t3), "_250m.sdat")
VBF.lst <- paste0("/data/MDEM/MVBF_", names(equi7t3), "_250m.sdat")
VDP.lst <- paste0("/data/MDEM/MVDP_", names(equi7t3), "_250m.sdat")
TWI.lst <- paste0("/data/MDEM/MTWI_", names(equi7t3), "_250m.sdat")
DEM.in.lst <- list(DEM.lst, SLP.lst, CRV.lst, VBF.lst, DVM.lst, VDP.lst, POS.lst, NEG.lst, TWI.lst)
EVIM.out.lst <- paste0("EX",1:6, "MOD5")
EVIS.out.lst <- paste0("ES",1:6, "MOD5")
EVIM.lst <- paste0("/data/MOD13Q1/SD_liste", ms.lst, ".vrt")
EVIS.lst <- paste0("/data/MOD13Q1/M_liste", ms.lst, ".vrt")
NBR4.out.lst <- c(paste0("I0",1:9, "MOD4"), paste0("I",10:12, "MOD4"))
NBR7.out.lst <- c(paste0("M0",1:9, "MOD4"), paste0("M",10:12, "MOD4"))
NBR4.lst <- paste0("/data/MCD43A4/NRB4_Ms_", m.lst, "_500m.tif")
NBR7.lst <- paste0("/data/MCD43A4/NRB7_Ms_", m.lst, "_500m.tif")
LSTD.out.lst <- c(paste0("T0",1:9, "MOD3"), paste0("T",10:12, "MOD3"), "TMDMOD3", "TSDMOD3")
LSTN.out.lst <- c(paste0("N0",1:9, "MOD3"), paste0("N",10:12, "MOD3"), "TMNMOD3", "TSNMOD3")
LSTD.lst <- c(paste0("/data/MOD11A2/M_liste", m.lst, "_D.vrt"), "/data/MOD11A2/TMDMOD3.tif", "/data/MOD11A2/TSDMOD3.tif")
LSTN.lst <- c(paste0("/data/MOD11A2/M_liste", m.lst, "_N.vrt"), "/data/MOD11A2/TMNMOD3.tif", "/data/MOD11A2/TSNMOD3.tif")
SNW.out.lst <- paste0("SN",1:6, "MOD4")
SNW.lst <- paste0("/data/MOD10A2/SNW_Ms_", ms.lst, "_500m.sdat")
PREm.out.lst <- c(paste0("P0",1:9, "MRG3"), paste0("P",10:12, "MRG3"))
PREm.lst <- paste0("/data/PREm/PREm_", 1:12, "_1km_sum.sdat")
lits <- c(paste0("0", 1:9), paste0(c(10:11,13:16)))
LIT.out.lst <- c(paste0("L0",1:9, "USG5"), paste0("L",c(10:11,13:16), "USG5"))
LFO.out.lst <- paste0("F0",1:7, "USG5")
LIT.lst <- lapply(lits, function(x){paste0("/data/EcoTapestry/equi7/L", x, "USG_", names(equi7t3), "_250m.tif")})
LFO.lst <- lapply(1:7, function(x){paste0("/data/EcoTapestry/equi7/F0", x, "USG_", names(equi7t3), "_250m.tif")})

## test preparing covs for 8 sample areas (single tile):
names(equi7t3)
#prepareCovsSoilGrids250m(s.zone=5, s.lst=145)
prepareCovsSoilGrids250m(s.zone=1, s.lst=rownames(equi7t3[[1]]@data)[which(equi7t3[[1]]$TILE=="072_048")])
prepareCovsSoilGrids250m(s.zone=1, s.lst=rownames(equi7t3[[1]]@data)[which(equi7t3[[1]]$TILE=="012_063")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=6, s.lst=rownames(equi7t3[[6]]@data)[which(equi7t3[[6]]$TILE=="087_063")])
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="072_087")])
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="048_003")])
prepareCovsSoilGrids250m(s.zone=4, s.lst=rownames(equi7t3[[4]]@data)[which(equi7t3[[4]]$TILE=="039_039")], close.gap=FALSE)
## Corrupt file?? Fill gaps operation hangs in loops
rownames(equi7t3[[4]]@data)[which(equi7t3[[4]]$TILE=="039_039")]
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="060_030")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="054_069")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="099_057")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="084_057")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="087_072")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="102_018")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="105_069")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="102_069")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="102_084")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="084_057")], close.gap=FALSE)
## Corrupt files??
prepareCovsSoilGrids250m(s.zone=7, s.lst=rownames(equi7t3[[7]]@data)[which(equi7t3[[7]]$TILE=="063_030")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=7, s.lst=rownames(equi7t3[[7]]@data)[which(equi7t3[[7]]$TILE=="063_012")], close.gap=FALSE)
## Corrupt files??
prepareCovsSoilGrids250m(s.zone=7, s.lst=rownames(equi7t3[[7]]@data)[which(equi7t3[[7]]$TILE=="090_048")])

## all together (TAKES >3 DAYS!):
for(i in 1:7){
  sfInit(parallel=TRUE, cpus=40)
  sfLibrary(RSAGA)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfLibrary(R.utils)
  sfExportAll()
  lst <- 1:length(equi7t3[[i]])
  x <- sfLapply(lst, function(x){ prepareCovsSoilGrids250m(s.zone=i, s.lst=x, close.gap=TRUE) })
  sfStop()
}

#prepareCovsSoilGrids250m(s.zone=1, s.lst=1:length(equi7t3[[1]]))
#prepareCovsSoilGrids250m(s.zone=3, s.lst=1:length(equi7t3[[3]]))
#prepareCovsSoilGrids250m(s.zone=4, s.lst=1:length(equi7t3[[4]]))
#prepareCovsSoilGrids250m(s.zone=5, s.lst=1:length(equi7t3[[5]]))
#prepareCovsSoilGrids250m(s.zone=6, s.lst=1:length(equi7t3[[6]]))
#prepareCovsSoilGrids250m(s.zone=7, s.lst=1:length(equi7t3[[7]]))

## plot in GE:
tif.lst <- list.files(path="./NA_060_036/", pattern=glob2rx("*_*_*_*.tif$"), full.names=TRUE)
des <- read.csv("SoilGrids250m_COVS250m.csv")
varn <- sapply(basename(tif.lst), function(x){strsplit(x, "_")[[1]][1]})
no <- match(varn, des$WORLDGRIDS_CODE)
plotKML.env(convert="convert", show.env=FALSE)
## this takes some time...
s <- raster::stack(tif.lst)
plotCovsSoilGrids250m(s, path="./KML/", ATTRIBUTE_TITLE=paste(des$ATTRIBUTE_TITLE[no]), DESCRIPTION=paste(des$DESCRIPTION[no]), EVI_range=c(120,2300), ES_range=c(0,4200), TD_range=c(265,295), TN_range=c(265,295), NBR4_range=c(300,7500), NBR7_range=c(60,1500), SNW_range=c(100,630))

## Prepare RDA files (one file per tile):
make.covsRDA(in.path="/data/covs", i="NA_060_036")

pr.dirs <- basename(list.dirs("/data/predicted")[-1])
sfInit(parallel=TRUE, cpus=35)
sfExport("make.covsRDA")
sfLibrary(rgdal)
sfLibrary(sp)
x <- sfLapply(pr.dirs, fun=function(i){ try(make.covsRDA(i, in.path="/data/covs") ) })
sfStop()

## clean-up:
#del.lst <- list.files(pattern=glob2rx("VDPMRG5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
