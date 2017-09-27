## Distribution of acid and saline/sodic soils based on SoilGrids250m
## Requested by the International Fertilizer Association (IFA)
## Code by Tom.Hengl@isric.org, decision rules by Niels Batjes (niels.batjes@isric.org) and Tom Hengl
## Compare with previous global assessments by: Wicke et al. "The global technical and economic potential of bioenergy from salt-affected soils" DOI: 10.1039/C1EE01029H (Analysis) Energy Environ. Sci., 2011, 4, 2669-2681

setwd("/data/models")
library(rgdal)
library(raster)
library(GSIF)
library(snowfall)

fao.lst <- c("Gleyic.Solonetz", "Mollic.Solonetz", "Solodic.Planosols", "Calcic.Solonetz", "Haplic.Solonchaks..Sodic.", "Calcic.Gypsisols", "Haplic.Cambisols..Sodic.", "Haplic.Calcisols..Sodic.", "Gypsic.Solonchaks", "Haplic.Regosols..Sodic.")
fao.grades = c(4, 3, 2, 4, 4, 4, 1, 3, 4, 2)
## https://www.blogs.nrcs.usda.gov/wps/PA_NRCSConsumption/download?cid=nrcseprd589210&ext=pdf
## Sodic Soils have maybe low levels of neutral soluble salts (ECe > 4.0 dS/m), they have relatively high levels of sodium on the exchange complex (ESP and SAR values are above 15 and 13, respectively). 
## The pH values of sodic soils exceed 8.5, rising to 10 or higher in some cases.

sodic_grade <- function(i, in.path, fao.lst, fao.grades, ph_t1=85, ph_t2=81){
  out.p <- paste0(in.path, "/", i, "/SLGWRB_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- c(paste0(in.path, "/", i, "/TAXNWRB_", fao.lst, "_", i, ".tif"), paste0(in.path, "/", i, "/PHIHOX_M_sl",1:5,"_", i, ".tif"))
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    names(s) <- c(fao.lst, paste0("PHIHOX_M_sl",1:5))
    gc()
    s$SLGWRB <- round(rowSums( data.frame(mapply('*', s@data[,fao.lst], fao.grades, SIMPLIFY=FALSE)) )/100)
    ph = rowMeans(s@data[,paste0("PHIHOX_M_sl",1:5)])
    s$SLGWRB <- ifelse(ph>ph_t1, 4, ifelse(ph>ph_t2, 3, s$SLGWRB))
    writeGDAL(s["SLGWRB"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

#sodic_grade(i="T33349", in.path="/data/tt/SoilGrids250m/predicted250m", fao.lst, fao.grades)
#sodic_grade(i="T10410", in.path="/data/tt/SoilGrids250m/predicted250m", fao.lst, fao.grades)

## Acidic subsoils:
PHIHOX_U = c(20, 45, 50, 55, 66, 73, 110)
Ca_deficiency = c("Extremely high", "Very high", "High", "Moderate", "Low", "Low")
## Classes (http://www.fao.org/soils-portal/soil-management/management-of-some-problem-soils/acid-soils/en/):
fao.acid = c("Haplic.Cambisols..Dystric.", "Haplic.Fluvisols..Dystric.", "Haplic.Gleysols..Dystric.", "Haplic.Planosols..Dystric.", "Haplic.Regosols..Dystric.", "Haplic.Alisols", "Cutanic.Alisols", "Haplic.Acrisols..Alumic.", "Haplic.Acrisols", "Haplic.Acrisols..Ferric.", "Haplic.Acrisols..Humic.", "Plinthic.Acrisols", "Vetic.Acrisols", "Gleyic.Podzols", "Haplic.Podzols")
fao.acid.grades = c(0.5,1.5,2,1,1,3,3,4,3,3,3,4,4,2.5,2)

## Derive acidity grade using soil pH and soil types:
wrapper.ACID <- function(i, in.path, fao.acid, fao.acid.grades, PHIHOX_U){
  out.p <- paste0(in.path, "/", i, "/ACDWRB_M_ss_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- c(paste0(in.path, "/", i, "/TAXNWRB_", fao.acid, "_", i, ".tif"), paste0(in.path, "/", i, "/PHIHOX_M_sl",3:5,"_", i, ".tif"))
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    names(s) <- c(fao.acid, paste0("PHIHOX_M_sl",3:5))
    wrb = rowSums( data.frame(mapply('*', s@data[,fao.acid], fao.acid.grades, SIMPLIFY=FALSE)) )/100
    ph = rowMeans(s@data[,paste0("PHIHOX_M_sl",3:5)])
    ph = ifelse(ph<PHIHOX_U[2], 4, ifelse(ph<PHIHOX_U[3], 3, ifelse(ph<PHIHOX_U[4], 2, ifelse(ph<PHIHOX_U[5], 1, 0))))
    ## round up to a higher number (to be on safe side)
    s$ACDWRB <- round(apply(data.frame(wrb, ph), 1, max))
    writeGDAL(s["ACDWRB"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
    gc(); gc()
  }
}

#wrapper.ACID(i="T38716", in.path="/data/tt/SoilGrids250m/predicted250m", fao.acid, fao.acid.grades, PHIHOX_U)
#wrapper.ACID(i="T10410", in.path="/data/tt/SoilGrids250m/predicted250m", fao.acid, fao.acid.grades, PHIHOX_U)

# clean-up:
#for(i in c("ACDWRB_M_ss", "SLGWRB")){  
#   del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
#}

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1])

sfInit(parallel=TRUE, cpus=48)
sfExport("sodic_grade", "fao.lst", "fao.grades")
sfLibrary(raster)
sfLibrary(rgdal)
out <- sfClusterApplyLB(pr.dirs, function(i){try( sodic_grade(i, in.path="/data/tt/SoilGrids250m/predicted250m", fao.lst, fao.grades) )})
sfStop()

sfInit(parallel=TRUE, cpus=48)
sfExport("wrapper.ACID", "fao.acid", "fao.acid.grades", "PHIHOX_U")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(GSIF)
out <- sfClusterApplyLB(pr.dirs, function(i){try( wrapper.ACID(i, in.path="/data/tt/SoilGrids250m/predicted250m", fao.acid, fao.acid.grades, PHIHOX_U) )})
sfStop()

## metadata:
metasd <- read.csv('/data/GEOG/META_GEOTIFF_1B.csv', stringsAsFactors = FALSE)
sel.metasd = names(metasd)[-sapply(c("FileName","VARIABLE_NAME"), function(x){grep(x, names(metasd))})]
source("mosaick_functions_ll.R")

## Make mosaics -----
t.vars = c("SLGWRB", "ACDWRB_M_ss")
r <- raster("/data/stacked250m/LCEE10.tif")
cellsize = res(r)[1]
te = paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")

library(snowfall)
sfInit(parallel=TRUE, cpus=ifelse(length(t.vars)>45, 45, length(t.vars)))
sfExport("t.vars", "make_mosaick_ll", "metasd", "sel.metasd", "cellsize", "te")
out <- sfClusterApplyLB(1:length(t.vars), function(x){ try( make_mosaick_ll(varn=t.vars[x], in.path="/data/tt/SoilGrids250m/predicted250m", tr=cellsize, te=te, ot="Byte", dstnodata=255, metadata=metasd[which(metasd$FileName == paste0(t.vars[x], "_250m_ll.tif")), sel.metasd]) )})
sfStop()
