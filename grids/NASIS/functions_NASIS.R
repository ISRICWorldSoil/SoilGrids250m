## Functions to prepare 100 m covariates and generate predictions
## Tom.Hengl@isric.org


make_RDS_tiles <- function(i, tile.tbl, in.path="/data/GEOG/NASIS", out.path="/data/NASIS/covs100m", mask="/data/GEOG/NASIS/COUNTY6.tif", covs.lst, LNDCOV6.leg, PMTGSS7.leg, DRNGSS7.leg, PVEGKT6.leg, COUNTY6.leg, GESUSG6.leg){
  out.rda <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/T", tile.tbl[i,"ID"], ".rds")
  if(!file.exists(out.rda)){
    m = readGDAL(fname=mask, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    names(m)[1] = "COUNTY6"
    m = as(m, "SpatialPixelsDataFrame")
    x = spTransform(m, CRS("+proj=longlat +datum=WGS84"))
    m$LONWGS84 <- x@coords[,1]
    m$LATWGS84 <- x@coords[,2]
    for(j in 1:length(covs.lst)){
      m@data[,covs.lst[j]] <- signif(readGDAL(paste0(in.path, "/", covs.lst[j], ".tif"), offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index], 4)
    }
    ## factors:
    m@data[,"LNDCOV6"] = plyr::join(data.frame(Value=m$LNDCOV6), LNDCOV6.leg, type="left")$LC_class
    m@data[,"PMTGSS7"] = plyr::join(data.frame(Value=m$PMTGSS7), PMTGSS7.leg, type="left")$pmaterial_class_f
    m@data[,"DRNGSS7"] = plyr::join(data.frame(Value=m$DRNGSS7), DRNGSS7.leg, type="left")$drainage_class
    m@data[,"PVEGKT6"] = plyr::join(data.frame(X=m$PVEGKT6), PVEGKT6.leg, type="left")$NAMES
    m@data[,"COUNTY6"] = plyr::join(data.frame(X=m$COUNTY6), COUNTY6.leg, type="left")$NAMES
    m@data[,"GESUSG6.f"] = plyr::join(data.frame(X=m$GESUSG6), GESUSG6.leg, type="left")$NAMES
    m = m[complete.cases(m@data[,c("DEMNED6","GESUSG6.f")]),]
    ## Fix missing pixels in parent material and drainage class:
    rf.pm = readRDS("/data/USA48/RF_parent_material.rds")
    rf.dr = readRDS("/data/USA48/RF_drainage_class.rds")
    out = predict.rfsrc(rf.pm, m@data[,rf.pm$xvar.names], na.action="na.impute", importance=FALSE, membership=FALSE)
    if(length(out$class)==nrow(m)){
      m$PMTGSS7.f = paste(out$class)
      g = m["PMTGSS7.f"]
      g$INT = plyr::join(data.frame(pmaterial_class_f=g$PMTGSS7.f), PMTGSS7.leg, type="left", match="first")$Value
      writeGDAL(g["INT"], paste0(out.path, "/T", tile.tbl[i,"ID"], "/PMTGSS7_predicted_", tile.tbl[i,"ID"], ".tif"), options="COMPRESS=DEFLATE")
      m$PMTGSS7 = ifelse(is.na(m$PMTGSS7), paste(m$PMTGSS7.f), paste(m$PMTGSS7))
    }
    out = predict.rfsrc(rf.dr, m@data[,rf.dr$xvar.names], na.action="na.impute", importance=FALSE, membership=FALSE)
    if(length(out$class)==nrow(m)){
      m$DRNGSS7.f = paste(out$class)
      g = m["DRNGSS7.f"]
      g$INT = plyr::join(data.frame(drainage_class=g$DRNGSS7.f), DRNGSS7.leg, type="left", match="first")$Value
      writeGDAL(g["INT"], paste0(out.path, "/T", tile.tbl[i,"ID"], "/DRNGSS7_predicted_", tile.tbl[i,"ID"], ".tif"), options="COMPRESS=DEFLATE")
      m$DRNGSS7 = ifelse(is.na(m$DRNGSS7), paste(m$DRNGSS7.f), paste(m$DRNGSS7))
    }
    ## Save output:
    saveRDS(m, file=out.rda, compress=TRUE)
    gc()
  }
}


extract.tiled <- function(x, tile.pol, path="/data/NASIS/covs100m", ID="ID", cpus=48){
  x$row.index <- 1:nrow(x)
  ov.c <- over(spTransform(x, CRS(proj4string(tile.pol))), tile.pol)
  ov.t <- which(!is.na(ov.c[,ID]))
  ## for each point get the tile name:
  ov.c <- data.frame(ID=ov.c[ov.t,ID], row.index=ov.t)
  tiles.lst <- basename(dirname(list.files(path=path, pattern=glob2rx("*.rds$"), recursive=TRUE)))
  ov.c <- ov.c[ov.c[,ID] %in% tiles.lst,]
  tiles <- levels(as.factor(paste(ov.c[,ID]))) 
  cov.c <- as.list(tiles)
  names(cov.c) <- tiles
  ## extract using snowfall
  sfInit(parallel=TRUE, cpus=cpus)
  sfExport("x", "path", "ov.c", "ID", "cov.c", ".extract.tile", "tile.pol")
  sfLibrary(raster)
  sfLibrary(rgdal)
  ov.lst <- sfLapply(1:length(cov.c), function(i){try(.extract.tile(i, x=x, ID=ID, path=path, ov.c=ov.c, cov.c=cov.c, tile.pol=tile.pol), silent = TRUE)}) 
  snowfall::sfStop()
  ## bind together:
  out <- dplyr::rbind_all(ov.lst)
  out <- plyr::join(x@data, as.data.frame(out), type="left", by="row.index")
  return(out)
}

.extract.tile <- function(i, x, ID, path, ov.c, cov.c, tile.pol){
  row.index <- ov.c$row.index[ov.c[,ID]==names(cov.c)[i]]
  pnts <- x[row.index,]
  pnts <- spTransform(pnts, CRS(proj4string(tile.pol)))
  m <- readRDS(paste0(path, "/", names(cov.c)[i], "/", names(cov.c)[i], ".rds"))
  out <- sp::over(y=m, x=pnts)
  out$band1 <- NULL
  out$row.index <- row.index
  xy <- data.frame(pnts@coords)
  names(xy) <- c("X","Y")
  out <- cbind(out, xy)
  return(out)
}


## Predict classes in loop:
predict_factor_tile <- function(i, mRF, varn, levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", PMTGSS7.leg, DRNGSS7.leg){
  options(rf.cores=15, mc.cores=15)
  gc(); gc()
  out.c <- paste0(out.path, "/T", i, "/", varn, "_", levs, "_T", i, ".tif")
  if(any(!file.exists(out.c))){
    m <- readRDS(paste0(in.path, "/T", i, "/T", i, ".rds"))
    if(is.character(mRF)){
      mRF = readRDS.gz(mRF)
    }
    m$PMTGSS7 = factor(m$PMTGSS7, levels=unique(PMTGSS7.leg$pmaterial_class_f))
    m$DRNGSS7 = factor(m$DRNGSS7, levels=unique(DRNGSS7.leg$drainage_class))
    m$ID = m@grid.index
    m.PVEGKT6 = data.frame(model.matrix(~PVEGKT6-1, m@data))
    m.PVEGKT6$ID = row.names(m.PVEGKT6)
    m.LNDCOV6 = data.frame(model.matrix(~LNDCOV6-1, m@data))
    m.LNDCOV6$ID = row.names(m.LNDCOV6)
    m.PMTGSS7 = data.frame(model.matrix(~PMTGSS7-1, m@data))
    m.PMTGSS7$ID = row.names(m.PMTGSS7)
    m.DRNGSS7 = data.frame(model.matrix(~DRNGSS7-1, m@data))
    m.DRNGSS7$ID = row.names(m.DRNGSS7)
    m@data = plyr::join_all(list(m@data, m.PVEGKT6, m.LNDCOV6, m.PMTGSS7, m.DRNGSS7), by="ID")[,mRF$xvar.names]
    v = predict.rfsrc(mRF, m@data, na.action="na.impute", importance=FALSE, membership=TRUE)
    m@data = data.frame(v$predicted)*100
    ## probabilities:
    for(j in 1:ncol(m)){
      out <- paste0(out.path, "/T", i, "/", varn, "_", names(m)[j], "_T", i, ".tif")
      writeGDAL(m[j], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
    }
    ## dominant class
    m@data[,varn] = as.integer(v$class)
    writeGDAL(m[varn], paste0(out.path, "/T", i, "/", varn, "_T", i, ".tif"), type="Int16", mvFlag=32767, options="COMPRESS=DEFLATE")
    gc(); gc()
  }  
}

## Create mosaics:
mosaic_tiles_100m <- function(j, in.path, varn, r="bilinear", ot="Byte", dstnodata=255, out.path, te){
  out.tif <- paste0(out.path, varn, "_", j, "_100m.tif")
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -r \"', r, '\" -te ', te))
  }
}
