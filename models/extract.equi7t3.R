
extract.equi7t3 <- function(x, y, path=".", equi7t3, ID="SOURCEID", cpus, silent=TRUE){
  x$row.index <- 1:nrow(x)
  ## get continents: 
  ov.c <- lapply(equi7t3, function(t){over(spTransform(x, CRS(proj4string(t))), t)})
  ov.t <- lapply(ov.c, function(t){which(!is.na(t$TILE))})
  ## get tile NAMES:
  ov.c <- lapply(1:length(equi7t3), function(i){cbind(ov.c[[i]][ov.t[[i]],c("SHORTNAME","TILE")], row.index=ov.t[[i]])})
  ov.c <- do.call(rbind, ov.c)
  ov.c$TILENAME <- paste0(substr(ov.c$SHORTNAME, 5, 6), "_", ov.c$TILE)
  ## list all available tifs:
  cov.lst <- as.vector(unlist(lapply(y, function(i){list.files(path=path, pattern=glob2rx(paste0(i, "_*_*_*.tif$")), recursive=TRUE, full.names=TRUE)})))
  ## remove duplicated points (due to overlap in Equi7t3):
  ov.c <- ov.c[!duplicated(ov.c$row.index),]
  tiles <- levels(as.factor(ov.c$TILENAME))    
  cov.c <- lapply(tiles, function(i){grep(i, cov.lst)})
  names(cov.c) <- tiles 
  cov.c <- cov.c[lapply(cov.c,length)>0]
  ## extract using snowfall
  snowfall::sfInit(parallel=TRUE, cpus=cpus)
  snowfall::sfExport("x", "ov.c", "cov.c", "cov.lst", ".extract.tile", "equi7t3")
  snowfall::sfLibrary(package="raster", character.only=TRUE)
  snowfall::sfLibrary(package="rgdal", character.only=TRUE)
  ov.lst <- snowfall::sfClusterApplyLB(1:length(cov.c), .extract.tile) 
  snowfall::sfStop()
  ## bind:
  out <- dplyr::rbind_all(ov.lst)
  out <- plyr::join(x@data, as.data.frame(out), type="left", by="row.index")
  return(out)
}

.extract.tile <- function(i){
  s <- stack(cov.lst[cov.c[[i]]])
  row.index <- ov.c$row.index[ov.c$TILENAME==names(cov.c)[i]]
  pnts <- x[row.index,]
  prj <- CRS(proj4string(equi7t3[[strsplit(names(cov.c)[i], "_")[[1]][1]]]))
  pnts <- spTransform(pnts, prj)
  out <- data.frame(extract(s, pnts))
  ## shorten the name to allow for binding
  names(out) <- sapply(names(out), function(i){strsplit(i, "_")[[1]][1]})
  out$Equi7t3 <- names(cov.c)[i]
  out$row.index <- row.index
  xy <- data.frame(pnts@coords)
  names(xy) <- c("X","Y")
  out <- cbind(out, xy)
  return(out)
}