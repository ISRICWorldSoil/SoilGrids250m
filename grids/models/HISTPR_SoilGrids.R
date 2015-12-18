## Distribution of organic soils (probability of histosols) based on SoilGrids250m
## Tom.Hengl@isric.org

library(rgdal)
library(raster)
library(GSIF)
fao.lst <- c("Sapric.Histosols", "Hemic.Histosols", "Fibric.Histosols", "Cryic.Histosols", "Histic.Albeluvisols")
usda.lst <- c("Saprists", "Hemists", "Folists", "Fibrists")

histosol.prob <- function(i, in.path, fao.lst, usda.lst){
  out.p <- paste0(in.path, "/", i, "/HISTPR_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- c(paste0(in.path, "/", i, "/TAXNWRB_", fao.lst, "_", i, ".tif"), paste0(in.path, "/", i, "/TAXOUSDA_", usda.lst, "_", i, ".tif"))
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    names(s) <- c(fao.lst, usda.lst)
    gc()
    s$HISTPR <- (s@data[,"Sapric.Histosols"] + s@data[,"Saprists"])/2 + (s@data[,"Hemic.Histosols"] + s@data[,"Hemists"])/2 + (s@data[,"Fibric.Histosols"] + s@data[,"Fibrists"])/2 + s@data[,"Cryic.Histosols"] + s@data[,"Histic.Albeluvisols"] + s@data[,"Folists"]
    writeGDAL(s["HISTPR"], out.p, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
  }
}

#histosol.prob(i="SA_051_069", in.path="/data/predicted", fao.lst, usda.lst)

## Organic carbon stock:
wrapper.OCSTHA <- function(i, in.path, n.lst=c("ORCDRC","BLD","CRFVOL"), ORCDRC.sd=5, BLD.sd=100, CRFVOL.sd=4){
  out.all <- paste0(in.path, "/", i, "/OCSTHA_M_sd", 1:6, "_", i,".tif")
  if(any(!file.exists(out.all))){
    for(d in 1:6){
      tif.lst <- paste0(in.path, "/", i, "/", n.lst, "_M_sd", d, "_", i, ".tif")
      s <- raster::stack(tif.lst)
      s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
      ## Predict organic carbon stock (in tones / ha):
      s$v <- round(as.vector(OCSKGM(ORCDRC=s@data[,1], BLD=s@data[,2], CRFVOL=s@data[,3], HSIZE=get("stsize", envir = GSIF.opts)[d]*100, ORCDRC.sd=ORCDRC.sd, BLD.sd=BLD.sd, CRFVOL.sd=CRFVOL.sd)*10))
      writeGDAL(s["v"], out.all[d], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
  }
}

#wrapper.OCSTHA(i="NA_060_036", in.path="/data/predicted")

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/predicted1km")[-1])

sfInit(parallel=TRUE, cpus=40)
sfExport("histosol.prob", "fao.lst", "usda.lst")
sfLibrary(raster)
sfLibrary(rgdal)
out <- sfClusterApplyLB(pr.dirs, function(i){try(histosol.prob(i, in.path="/data/predicted", fao.lst, usda.lst))})
sfStop()

sfInit(parallel=TRUE, cpus=40)
sfExport("wrapper.OCSTHA")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(GSIF)
out <- sfClusterApplyLB(pr.dirs, function(i){try(wrapper.OCSTHA(i, in.path="/data/predicted1km"))})
sfStop()
