## Prepare covariates to generate SoilGrids250m
## By: Tom.Hengl@isric.org

setwd("/data/models")
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
names(equi7t3)
source("tiler.R")
source("prepareCovs.R")
source("plotCovsSoilGrids250m.R")
source("make.covsRDA.R")
## Months of interest:
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
ms.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
msk.lst <- list.files(path="/data/covs", pattern=glob2rx("LMK_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
## 2594 tiles

## remove all tiles with all missing pixels in the LMK
check.LMK <- function(i){
  r <- raster(i)
  na.count <- sum(getValues(is.na(r)|r==0))
  if(na.count==ncell(r)){
    return(0)
  } else {
    if(na.count==0){
      return(100)
    } else {
      return(signif((1-na.count/ncell(r))*100, 3))
    } 
  }
}
sfInit(parallel=TRUE, cpus=35)
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("check.LMK", "msk.lst")
selL <- sfLapply(msk.lst, check.LMK)
sfStop()
## 243 empty tiles on 1th November 2015
selD <- data.frame(name=msk.lst[which(selL==0)])
write.csv(selD, "empty_LMK_tiles.csv")
## remove all directories with empty landmask (CAREFULL!)
x = sapply(selD$name, function(x){unlink(dirname(as.character(x)), recursive = TRUE, force = TRUE)})
## 2356 dirs left
## Move land cover images:
GLC.out.lst <- paste0("C0", 1:9, "GLC5")
for(i in 1:length(msk.lst)){
  dirn <- basename(dirname(msk.lst[i]))
  for(j in 1:length(GLC.out.lst)){
    filen <- paste0(GLC.out.lst[j], "_", dirn, ".tif")
    from <- paste0("/data/GlobCover30/Mtiled/GLC", j*10, "_", dirn, ".tif")
    to <- paste0("/data/covs/", dirn, "/", filen)
    if(file.exists(from)){ file.rename(from=from, to=to) } 
  }
}

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
LSTD.out.lst <- c(paste0("T0",1:9, "MOD3"), paste0("T",10:12, "MOD3"), "TMDMOD3")
LSTN.out.lst <- c(paste0("N0",1:9, "MOD3"), paste0("N",10:12, "MOD3"), "TMNMOD3")
LSTD.lst <- c(paste0("/data/MOD11A2/LSTD_M_", m.lst, "_1km.tif"), "/data/MOD11A2/TMDMOD3.tif")
LSTN.lst <- c(paste0("/data/MOD11A2/LSTN_M_", m.lst, "_1km.tif"), "/data/MOD11A2/TMNMOD3.tif")
LSSD.out.lst <- c(paste0("T0",1:9, "MSD3"), paste0("T",10:12, "MSD3"))
LSSN.out.lst <- c(paste0("N0",1:9, "MSD3"), paste0("N",10:12, "MSD3"))
LSSD.lst <- paste0("/data/MOD11A2/LSTD_SD_", m.lst, "_1km.tif")
LSSN.lst <- paste0("/data/MOD11A2/LSTN_SD_", m.lst, "_1kmf.tif") ## filtered images!
SNW.out.lst <- paste0("SN",1:6,"MOD4")
SNW.lst <- paste0("/data/MOD10A2/Ms_liste", ms.lst, ".vrt")
PREm.out.lst <- c(paste0("P0",1:9, "MRG3"), paste0("P",10:12, "MRG3"))
PREm.lst <- paste0("/data/PREm/PREm_", 1:12, "_1km_sum.sdat")
lits <- c(paste0("0", 1:9), paste0(c(10:11,13:16)))
LIT.out.lst <- c(paste0("L0",1:9, "USG5"), paste0("L",c(10:11,13:16), "USG5"))
LFO.out.lst <- paste0("F0",1:7, "USG5")
LIT.lst <- lapply(lits, function(x){paste0("/data/EcoTapestry/equi7/L", x, "USG_", names(equi7t3), "_250m.tif")})
LFO.lst <- lapply(1:7, function(x){paste0("/data/EcoTapestry/equi7/F0", x, "USG_", names(equi7t3), "_250m.tif")})
NEO.out.lst <- paste0("VW",1:6,"MOD1")
NEO.lst <- paste0("/data/NEO/SKY_WV_M_", ms.lst, "_10km.tif")
FW.out.lst <- paste0("FW",c(4:5),"MOD5") ## Not all months are available
FW.lst <- list(paste0("/data/floods/FW", 4, "MOD5_", names(equi7t3), ".tif"), paste0("/data/floods/FW", 5, "MOD5_", names(equi7t3), ".tif"))
  

## test preparing covs for 8 sample areas (single tile)
#prepareCovsSoilGrids250m(s.zone=5, s.lst=145)
prepareCovsSoilGrids250m(s.zone=1, s.lst=rownames(equi7t3[[1]]@data)[which(equi7t3[[1]]$TILE=="072_048")])
prepareCovsSoilGrids250m(s.zone=1, s.lst=rownames(equi7t3[[1]]@data)[which(equi7t3[[1]]$TILE=="012_063")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=6, s.lst=rownames(equi7t3[[6]]@data)[which(equi7t3[[6]]$TILE=="087_063")])
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="072_087")])
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="048_003")])
prepareCovsSoilGrids250m(s.zone=4, s.lst=rownames(equi7t3[[4]]@data)[which(equi7t3[[4]]$TILE=="039_039")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=4, s.lst=rownames(equi7t3[[4]]@data)[which(equi7t3[[4]]$TILE=="030_015")], close.gap=FALSE)
## Corrupt file?? Fill gaps operation hangs in loops
#rownames(equi7t3[[4]]@data)[which(equi7t3[[4]]$TILE=="039_039")]
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
## Corrupt files??
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="072_090")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=3, s.lst=rownames(equi7t3[[3]]@data)[which(equi7t3[[3]]$TILE=="075_087")], close.gap=FALSE)

## all together (TAKES 3-4hrs):
for(i in 1:7){
  sfInit(parallel=TRUE, cpus=40)
  sfLibrary(RSAGA)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfLibrary(R.utils)
  sfExportAll()
  lst <- 1:length(equi7t3[[i]])
  x <- sfLapply(lst, function(x){ try( prepareCovsSoilGrids250m(s.zone=i, s.lst=x, close.gap=TRUE) ) })
  sfStop()
}

#prepareCovsSoilGrids250m(s.zone=1, s.lst=1:length(equi7t3[[1]]))
#prepareCovsSoilGrids250m(s.zone=3, s.lst=1:length(equi7t3[[3]]))
#prepareCovsSoilGrids250m(s.zone=4, s.lst=1:length(equi7t3[[4]]))
#prepareCovsSoilGrids250m(s.zone=5, s.lst=1:length(equi7t3[[5]]))
#prepareCovsSoilGrids250m(s.zone=6, s.lst=1:length(equi7t3[[6]]))
#prepareCovsSoilGrids250m(s.zone=7, s.lst=1:length(equi7t3[[7]]))

## plot in GE:
# tif.lst <- list.files(path="/data/covs/NA_060_036/", pattern=glob2rx("*_*_*_*.tif$"), full.names=TRUE)
# des <- read.csv("SoilGrids250m_COVS250m.csv")
# varn <- sapply(basename(tif.lst), function(x){strsplit(x, "_")[[1]][1]})
# no <- match(varn, des$WORLDGRIDS_CODE)
# plotKML.env(convert="convert", show.env=FALSE)
# ## this takes some time...
# s <- raster::stack(tif.lst)
# plotCovsSoilGrids250m(s, path="/data/KML/", ATTRIBUTE_TITLE=paste(des$ATTRIBUTE_TITLE[no]), DESCRIPTION=paste(des$DESCRIPTION[no]), EVI_range=c(120,2300), ES_range=c(0,4200), TD_range=c(265,295), TN_range=c(265,295), NBR4_range=c(300,7500), NBR7_range=c(60,1500), SNW_range=c(100,630))

## Prepare RDA files (one file per tile):
make.covsRDA(in.path="/data/covs", i="AS_075_087")
#test <- readRDS("/data/covs/AS_075_087/AS_075_087.rds")
#plot(stack(test[15:25]),col=SAGA_pal[[1]])

## TAKES ca 2 hrs to generate all RDS files
pr.dirs <- basename(list.dirs("/data/covs")[-1])
## 2356 dirs
ok.lst <- basename(dirname(list.files(path="/data/predicted", pattern=glob2rx("TAXOUSDA_??_???_???.tif$"), full.names=TRUE, recursive=TRUE)))
del.lst <- pr.dirs[which(!pr.dirs %in% ok.lst)]
del.lst <- paste0('/data/covs/', del.lst, '/', del.lst,'.rds')
unlink(del.lst)
#x = paste0("/data/covs/", basename(dirname(del.lst)), "/N10MOD3_", basename(dirname(del.lst)),".tif")
#unlink(x)

sfInit(parallel=TRUE, cpus=40)
sfExport("make.covsRDA")
sfLibrary(rgdal)
sfLibrary(sp)
x <- sfLapply(pr.dirs, fun=function(i){ try(make.covsRDA(i, in.path="/data/covs") ) })
sfStop()

## create prediction dirs:
#x <- lapply(gsub("covs", "predicted", pr.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)

## clean-up:
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("C??GLC5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("ES?MOD5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("T??MSD3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("SN?MOD4_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("GTDHYS3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("*_*_*.rds"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)