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

#histosol.prob(i="SA_051_069", in.path="/data/tt/SoilGrids250m/predicted250m", fao.lst, usda.lst)

## Organic carbon stock (six standard layers) corrected for depth to bedrock:
wrapper.OCSTHA <- function(i, in.path, n.lst=c("ORCDRC","BLD","CRFVOL"), ORCDRC.sd=20, BLD.sd=100, CRFVOL.sd=5, BDR.lst=c("BDRICM","BDRLOG","BDTICM"), sdepth = c(0, 5, 15, 30, 60, 100, 200)){
  ## six standard layers 0-5, 5-15, 15-30, 30-60, 60-100, 100-200:
  out.all <- paste0(in.path, "/", i, "/OCSTHA_M_sd", 1:6, "_", i,".tif")
  if(any(!file.exists(out.all))){
    sD <- raster::stack(paste0(in.path, "/", i, "/", BDR.lst, "_M_", i, ".tif"))
    sD <- as(sD, "SpatialGridDataFrame")
    sD$BDRLOG <- ifelse(sD@data[,grep("BDRLOG_M", names(sD))]>50, 100, NA)
    ## parallel minimum:
    sD$BDRICM <- pmin(sD@data[,grep("BDRICM_M", names(sD))], sD@data[,grep("BDTICM_M", names(sD))], sD$BDRLOG, na.rm = TRUE)
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
      ## Correct for depth to bedrock:
      s$v <- ifelse(sD@data[s@grid.index,"BDRICM"] > sdepth[d+1], s$v, ifelse(sD@data[s@grid.index,"BDRICM"] > sdepth[d], s$v*(sD@data[s@grid.index,"BDRICM"]-sdepth[d])/(sdepth[d+1]-sdepth[d]), 0))
      writeGDAL(s["v"], out.all[d], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      writeGDAL(sD["BDRICM"], paste0(in.path, "/", i, "/BDRMIN_M_", i, ".tif"), type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      gc()
    }
  }
}

#wrapper.OCSTHA(i="NA_060_036", in.path="/data/tt/SoilGrids250m/predicted250m")

## derive OCS from organic carbon density (and depth to bedrock) -----
## this formula (average between two estimates) gives more or less best results although one should be very careful with using the predicted BLD
wrapper.ocd2OCSTHA <- function(i, in.path, BDR.lst=c("BDRICM","BDRLOG","BDTICM"), sdepth = c(0, 5, 15, 30, 60, 100, 200), st = c(30,100,200), n.lst=c("ORCDRC","BLD.f","CRFVOL"), ORCDRC.sd=20, BLD.sd=100, CRFVOL.sd=5){
  out.all <- paste0(in.path, "/", i, "/OCSTHA_M_sd", 1:6, "_", i,".tif")
  if(any(!file.exists(out.all))){
    ## depth to bedrock maps:
    sD <- raster::stack(paste0(in.path, "/", i, "/", BDR.lst, "_M_", i, ".tif"))
    sD <- as(sD, "SpatialGridDataFrame")
    sD$BDRLOG <- ifelse(sD@data[,grep("BDRLOG_M", names(sD))]>50, 100, 200)
    ## parallel minimum (use minimum because lower depths are typically over-estimated):
    sD$BDRICM <- pmin(sD@data[,grep("BDRICM_M", names(sD))], sD@data[,grep("BDTICM_M", names(sD))], sD$BDRLOG, na.rm = TRUE)
    ## Organic carbon density maps:
    s = stack(paste0(in.path, "/", i, "/OCDENS_M_sl", 1:7, "_", i,".tif"))
    s = as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    for(d in 1:(length(sdepth)-1)){
      ## Organic carbon content + BLD + CRF maps:
      Utif.lst <- paste0(in.path, "/", i, "/", n.lst, "_M_sl", d, "_", i, ".tif")
      Ltif.lst <- paste0(in.path, "/", i, "/", n.lst, "_M_sl", d+1, "_", i, ".tif")
      s0 <- raster::stack(c(Utif.lst,Ltif.lst))
      s0 <- as(s0, "SpatialGridDataFrame")
      s0$OC <- rowMeans(s0@data[,grep("ORCDRC", names(s0))], na.rm = TRUE)
      s0$BD <- rowMeans(s0@data[,grep("BLD.f", names(s0))], na.rm = TRUE)
      s0$CF <- rowMeans(s0@data[,grep("CRFVOL", names(s0))], na.rm = TRUE)
      ## tons per ha (average between OCS based on OCD and based on the OCS formula):
      v1 <- round(rowMeans(s@data[,paste0("OCDENS_M_sl", d:(d+1), "_", i)], na.rm=TRUE)*(sdepth[d+1]-sdepth[d])/100)
      v2 <- round(as.vector(OCSKGM(ORCDRC=s0$OC, BLD=s0$BD, CRFVOL=s0$CF, HSIZE=(sdepth[d+1]-sdepth[d]), ORCDRC.sd=ORCDRC.sd, BLD.sd=BLD.sd, CRFVOL.sd=CRFVOL.sd)*10))
      v = (v1+v2[s@grid.index])/2
      ## Correct (reduce OCS) for depth to bedrock:
      s@data[,paste0("OCSTHA_",d)] <- ifelse(sD@data[s@grid.index,"BDRICM"] > sdepth[d+1], v, ifelse(sD@data[s@grid.index,"BDRICM"] > sdepth[d], v*(sD@data[s@grid.index,"BDRICM"]-sdepth[d])/(sdepth[d+1]-sdepth[d]), 0))
      writeGDAL(s[paste0("OCSTHA_",d)], out.all[d], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
    ## cumulative OCS for standard depths 0-30, 0-100 and 0-200 cm:
    for(k in 1:length(st)){
      if(st[k]==200){
        s$SOCS = rowSums(s@data[,paste0("OCSTHA_",1:6)], na.rm=TRUE)
      }
      if(st[k]==100){
        s$SOCS = rowSums(s@data[,paste0("OCSTHA_",1:5)], na.rm=TRUE)
      }
      if(st[k]==30){
        s$SOCS = rowSums(s@data[,paste0("OCSTHA_",1:3)], na.rm=TRUE)
      }
      ## tones / ha
      writeGDAL(s["SOCS"], paste0(in.path, "/", i, "/OCSTHA_M_", st[k], "cm_", i,".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
  }
}

#wrapper.ocd2OCSTHA(i="T38715", in.path="/data/tt/SoilGrids250m/predicted250m")
#wrapper.ocd2OCSTHA(i="T10410", in.path="/data/tt/SoilGrids250m/predicted250m")
#wrapper.ocd2OCSTHA(i="T34387", in.path="/data/tt/SoilGrids250m/predicted250m")

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1])

sfInit(parallel=TRUE, cpus=56)
sfExport("histosol.prob", "fao.lst", "usda.lst")
sfLibrary(raster)
sfLibrary(rgdal)
out <- sfClusterApplyLB(pr.dirs, function(i){try( histosol.prob(i, in.path="/data/tt/SoilGrids250m/predicted250m", fao.lst, usda.lst) )})
sfStop()

sfInit(parallel=TRUE, cpus=56)
sfExport("wrapper.ocd2OCSTHA")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(GSIF)
out <- sfClusterApplyLB(pr.dirs, function(i){try( wrapper.ocd2OCSTHA(i, in.path="/data/tt/SoilGrids250m/predicted250m") )})
sfStop()

## clean-up:
# for(i in c("OCSTHA", "HISTPR")){  
#   del.lst <- list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx(paste0("^", i, "*.tif")), full.names=TRUE, recursive=TRUE)
#   unlink(del.lst)
# }

