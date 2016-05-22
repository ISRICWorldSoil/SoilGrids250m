## Distribution of organic soils (probability of histosols) based on SoilGrids250m
## Tom.Hengl@isric.org

library(rgdal)
library(raster)
library(GSIF)
fao.lst <- c("Calcic.Histosols", "Cryic.Histosols", "Fibric.Histosols", "Hemic.Histosols", "Sapric.Histosols", "Histic.Albeluvisols")
usda.lst <- c("Saprists", "Hemists", "Folists", "Fibrists")

histosol.prob <- function(i, in.path, fao.lst, usda.lst){
  out.p <- paste0(in.path, "/", i, "/HISTPR_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- c(paste0(in.path, "/", i, "/TAXNWRB_", fao.lst, "_", i, ".tif"), paste0(in.path, "/", i, "/TAXOUSDA_", usda.lst, "_", i, ".tif"))
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    names(s) <- c(fao.lst, usda.lst)
    gc()
    s$HISTPR <- (s@data[,"Sapric.Histosols"] + s@data[,"Saprists"])/2 + (s@data[,"Hemic.Histosols"] + s@data[,"Hemists"])/2 + (s@data[,"Fibric.Histosols"] + s@data[,"Fibrists"])/2 + s@data[,"Calcic.Histosols"] + s@data[,"Cryic.Histosols"] + s@data[,"Histic.Albeluvisols"] + s@data[,"Folists"]
    writeGDAL(s["HISTPR"], out.p, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
  }
}

#histosol.prob(i="SA_051_069", in.path="/data/predicted1km", fao.lst, usda.lst)

## Organic carbon stock (six standard layers):
wrapper.OCSTHA <- function(i, in.path, n.lst=c("ORCDRC","BLD","CRFVOL"), ORCDRC.sd=20, BLD.sd=100, CRFVOL.sd=5){
  ## six standard layers 0-5, 5-15, 15-30, 30-60, 60-100, 100-200:
  out.all <- paste0(in.path, "/", i, "/OCSTHA_M_sd", 1:6, "_", i,".tif")
  if(any(!file.exists(out.all))){
    for(d in 1:6){
      Utif.lst <- paste0(in.path, "/", i, "/", n.lst, "_M_sl", d, "_", i, ".tif")
      Ltif.lst <- paste0(in.path, "/", i, "/", n.lst, "_M_sl", d+1, "_", i, ".tif")
      s <- raster::stack(c(Utif.lst,Ltif.lst))
      s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
      s$ORCDRC <- rowMeans(s@data[,grep("ORCDRC", names(s))], na.rm = TRUE)
      s$BLD <- rowMeans(s@data[,grep("BLD", names(s))], na.rm = TRUE)
      s$CRFVOL <- rowMeans(s@data[,grep("CRFVOL", names(s))], na.rm = TRUE)
      ## Predict organic carbon stock (in tones / ha):
      s$v <- round(as.vector(OCSKGM(ORCDRC=s$ORCDRC, BLD=s$BLD, CRFVOL=s$CRFVOL, HSIZE=get("stsize", envir = GSIF.opts)[d]*100, ORCDRC.sd=ORCDRC.sd, BLD.sd=BLD.sd, CRFVOL.sd=CRFVOL.sd)*10))
      writeGDAL(s["v"], out.all[d], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      gc()
    }
  }
}

#wrapper.OCSTHA(i="NA_060_036", in.path="/data/predicted1km")

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/predicted")[-1])

sfInit(parallel=TRUE, cpus=48)
sfExport("histosol.prob", "fao.lst", "usda.lst")
sfLibrary(raster)
sfLibrary(rgdal)
out <- sfClusterApplyLB(pr.dirs, function(i){try( histosol.prob(i, in.path="/data/predicted", fao.lst, usda.lst) )})
sfStop()

sfInit(parallel=TRUE, cpus=48)
sfExport("wrapper.OCSTHA")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(GSIF)
out <- sfClusterApplyLB(pr.dirs, function(i){try( wrapper.OCSTHA(i, in.path="/data/predicted") )})
sfStop()

## clean-up:
# for(i in c("OCSTHA", "HISTPR")){  
#   del.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }
