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
length(msk.lst)
## 2599 tiles

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
sfInit(parallel=TRUE, cpus=48)
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("check.LMK", "msk.lst")
selL <- sfLapply(msk.lst, check.LMK)
sfStop()
## 243 empty tiles on 20th February 2016
selD <- data.frame(name=msk.lst[which(selL==0)])
write.csv(selD, "empty_LMK_tiles.csv")
## remove all directories with empty landmask (CAREFULL!)
x = sapply(selD$name, function(x){unlink(dirname(as.character(x)), recursive = TRUE, force = TRUE)})
## 2356 dirs left
## 3 tiles contain 1-2 pixels only
#rm.t = c("AS_084_045", "AS_081_075", "OC_069_102")

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
EVIM.lst <- paste0("/data/MOD13Q1/M_liste", ms.lst, ".vrt")
EVIS.lst <- paste0("/data/MOD13Q1/SD_liste", ms.lst, ".vrt")
NBR4.out.lst <- c(paste0("I0",1:9, "MOD4"), paste0("I",10:12, "MOD4"))
NBR7.out.lst <- c(paste0("M0",1:9, "MOD4"), paste0("M",10:12, "MOD4"))
NBR4.lst <- paste0("/data/MCD43A4/M_liste", m.lst, "_4.vrt")
NBR7.lst <- paste0("/data/MCD43A4/M_liste", m.lst, "_7.vrt")
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
## Precipitation:
PREm.out.lst <- c(paste0("P0",1:9, "MRG3"), paste0("P",10:12, "MRG3"), "PRSMRG3")
PREm.lst <- c(paste0("/data/PREm/P0",1:9, "MRG3a.tif"), paste0("/data/PREm/P",10:12, "MRG3a.tif"), "/data/PREm/PRSMRG3a.tif")
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
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="099_084")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="093_012")], close.gap=FALSE)
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="078_060")])
prepareCovsSoilGrids250m(s.zone=5, s.lst=rownames(equi7t3[[5]]@data)[which(equi7t3[[5]]$TILE=="075_072")])

## Tiles that usually need carefull checking:
c.ts <- c("NA_108_042", "NA_081_081", "NA_084_081", "NA_087_084", "NA_090_087", "OC_147_048", "EU_057_051")
m1km <- readRDS(paste0("/data/covs1km/", c.ts[1], "/", c.ts[1], ".rds"))

## all together (TAKES 3-4hrs):
for(i in 1:7){
  sfInit(parallel=TRUE, cpus=48)
  sfLibrary(RSAGA)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfLibrary(R.utils)
  sfExportAll()
  lst <- 1:length(equi7t3[[i]])
  x <- sfClusterApplyLB(lst, function(x){ try( prepareCovsSoilGrids250m(s.zone=i, s.lst=x) ) })
  sfStop()
}

#prepareCovsSoilGrids250m(s.zone=1, s.lst=1:length(equi7t3[[1]]))

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
make.covsRDA(in.path="/data/covs", i="NA_075_072")
#test <- readRDS("/data/covs/NA_075_072/NA_075_072.rds")
#plot(stack(test[29:35]),col=SAGA_pal[[1]], zlim=c(0,100))
#plot(stack(test[51:65]),col=SAGA_pal[[1]], zlim=c(0,100))
#plot(stack(test[145:160]),col=SAGA_pal[[1]])

## TAKES ca 2 hrs to generate all RDS files
pr.dirs <- basename(list.dirs("/data/covs")[-1])
## 2353 dirs
## Some tiles can not be used for prediction:
#ok.lst <- basename(dirname(list.files(path="/data/predicted", pattern=glob2rx("TAXOUSDA_??_???_???.tif$"), full.names=TRUE, recursive=TRUE)))
#del.lst <- pr.dirs[which(!pr.dirs %in% ok.lst)]
#del.lst <- paste0('/data/covs/', del.lst, '/', del.lst,'.rds')
#unlink(del.lst)
#x = paste0("/data/covs/", basename(dirname(del.lst)), "/N10MOD3_", basename(dirname(del.lst)),".tif")
#unlink(x)
sfInit(parallel=TRUE, cpus=48)
sfExport("make.covsRDA")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(R.utils)
x <- sfLapply(pr.dirs, fun=function(i){ try(make.covsRDA(i, in.path="/data/covs") ) })
sfStop()

## Alternative:
#gz files for h2o (ca 2 hrs to make):
#make.csv.gz(i="NA_060_036", in.path="/data/covs")

## create prediction dirs:
#x <- lapply(gsub("covs", "predicted", pr.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)

## ------------ 1 km ------------------

## Resample covs to 1km resolution (for testing purposes only!):
new.dirs <- list.dirs("/data/covs")[-1]
## create dirs for covs and predictions:
x <- lapply(gsub("covs", "covs1km", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
x <- lapply(gsub("covs", "predicted1km", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## function for fast resampling:
cov.lst <- as.vector(sapply(basename(list.files(path="/data/covs/NA_060_036", pattern=glob2rx("*.tif$"))), function(x){strsplit(x, "_")[[1]][1]}))
## 160 covs

resample_tif <- function(i, in.path="/data/covs", res=1000, out.dir="covs1km", cov.lst){
  in.file <- paste0(in.path, '/', i, '/', cov.lst, '_', i, '.tif')
  x <- lapply(in.file, function(x){ if(!file.exists(gsub("covs", out.dir, x))){ try( system(paste0(gdalwarp, ' ', x, ' ', gsub("covs", out.dir, x), ' -r \"near\" -tr ', res, ' ', res, ' -co \"COMPRESS=DEFLATE\"')) ) } }) 
}
#resample_tif(i="NA_060_036", cov.lst=cov.lst)
#resample_tif(i="NA_102_042", cov.lst=cov.lst)

## this takes only ca 20 mins with 'near' resampling
sfInit(parallel=TRUE, cpus=48)
sfExport("resample_tif", "cov.lst", "gdalwarp")
x <- sfClusterApplyLB(pr.dirs, fun=function(i){ resample_tif(i, cov.lst=cov.lst) })
sfStop()

## Some layers missing:
t.tif <- list.files(path="/data/covs1km", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
n.tifs <- lapply(cov.lst, grep, x=t.tif)
## Covariates with problems:
cov.incomplete = cov.lst[sapply(n.tifs, function(i){length(i)<2353})]
lapply(n.tifs[which(cov.lst %in% cov.incomplete)], length)

#make.covsRDA(i="NA_102_042", in.path="/data/covs1km")
## RDS files at 1 km
sfInit(parallel=TRUE, cpus=48)
sfExport("make.covsRDA")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(R.utils)
x <- sfClusterApplyLB(pr.dirs, fun=function(i){ try(make.covsRDA(i, in.path="/data/covs1km") ) })
sfStop()

## Check if everything is OK:
rds.lst <- list.files(path="/data/covs1km", pattern=glob2rx("*.rds$"), full.names=TRUE, recursive=TRUE)
pr.dirs[which(!pr.dirs %in% basename(dirname(rds.lst)))]
# Error in checkForRemoteErrors(val) : 
#   8 nodes produced errors; first error: Error in .local(.Object, ...) : 
#   `/data/covs1km/NA_099_000/M03MOD4_NA_099_000.tif' not recognised as a supported file format.

## ------------ Clean-up ------------------

## Total clean-up
#empty.lst <- list.files(path="/data/covs1km", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
#del.lst <- empty.lst[which(file.size(empty.lst)==0)]

#del.lst <- list.files(path="/data/covs", pattern=glob2rx("P??MRG3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("P??MRG3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("L??USG5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("L??USG5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("E??MOD5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("I??MOD4_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("M??MOD4_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("C??GLC5_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("TSDMOD3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("TSNMOD3_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
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
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("*_*_*.csv.gz"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("*_*_*.rds"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("*_*_*.csv.gz"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs", pattern=glob2rx("LMK_*_*_*.sdat$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/covs1km", pattern=glob2rx("LMK_*_*_*.*"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)