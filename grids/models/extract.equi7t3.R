## Overlay and extract values of covs on multiple tiles (EQUI7)
## Tom.Hengl@isric.org

extract.equi7t3 <- function(x, y, path=".", equi7t3, ID="SOURCEID", cpus, silent=TRUE, fromTif=FALSE, cov.lst=NULL){
  x$row.index <- 1:nrow(x)
  ## get continents: 
  ov.c <- lapply(equi7t3, function(t){over(spTransform(x, CRS(proj4string(t))), t)})
  ov.t <- lapply(ov.c, function(t){which(!is.na(t$TILE))})
  ## for each point get the tile name:
  ov.c <- lapply(1:length(equi7t3), function(i){cbind(ov.c[[i]][ov.t[[i]],c("SHORTNAME","TILE")], row.index=ov.t[[i]])})
  ov.c <- do.call(rbind, ov.c)
  ov.c$TILENAME <- paste0(sapply(as.character(ov.c$SHORTNAME), function(i){strsplit(i, " ")[[1]][2]}), "_", ov.c$TILE)
  if(fromTif==TRUE){
    tiles.lst <- basename(list.dirs(path=path)[-1])    
  } else {
    tiles.lst <- basename(dirname(list.files(path=path, pattern=glob2rx("*.rds$"), recursive=TRUE)))
  }
  ## 2356 tiles
  ov.c <- ov.c[ov.c$TILENAME %in% tiles.lst, ]
  ## remove duplicates (overlap in EQUI7 system)
  ov.c <- ov.c[!duplicated(ov.c$row.index),]
  tiles <- levels(as.factor(ov.c$TILENAME)) 
  if(fromTif==TRUE){
    ## parallelize:
    cov.c <- parallel::mclapply(tiles, function(i){grep(i, cov.lst)}, mc.cores=cpus)
    names(cov.c) <- tiles
    ## remove empty tiles
    cov.c <- cov.c[lapply(cov.c,length)>0]
  } else {
    cov.c <- as.list(tiles)
    names(cov.c) <- tiles
  }
  ## extract using snowfall
  snowfall::sfInit(parallel=TRUE, cpus=cpus)
  snowfall::sfExport("x", "path", "ov.c", "cov.c", "cov.lst", ".extract.tile", "equi7t3", "fromTif")
  snowfall::sfLibrary(package="raster", character.only=TRUE)
  snowfall::sfLibrary(package="rgdal", character.only=TRUE)
  ov.lst <- snowfall::sfClusterApplyLB(1:length(cov.c), function(i){try(.extract.tile(i, x=x, path=path, ov.c=ov.c, cov.c=cov.c, equi7t3=equi7t3, fromTif=fromTif))}) 
  snowfall::sfStop()
  ## bind:
  out <- dplyr::rbind_all(ov.lst)
  out <- plyr::join(x@data, as.data.frame(out), type="left", by="row.index")
  return(out)
}

.extract.tile <- function(i, x, path, ov.c, cov.c, equi7t3, fromTif=FALSE){
  row.index <- ov.c$row.index[ov.c$TILENAME==names(cov.c)[i]]
  pnts <- x[row.index,]
  prj <- CRS(proj4string(equi7t3[[strsplit(names(cov.c)[i], "_")[[1]][1]]]))
  pnts <- spTransform(pnts, prj)
  if(fromTif==TRUE){
    s <- stack(cov.lst[cov.c[[i]]])
    out <- data.frame(extract(s, pnts))
    ## shorten the name to allow for binding
    names(out) <- sapply(names(out), function(i){strsplit(i, "_")[[1]][1]})
  } else {
    m <- readRDS(paste0(path, "/", names(cov.c)[i], "/", names(cov.c)[i], ".rds"))
    out <- sp::over(y=m, x=pnts)
    out$band1 <- NULL
  }
  out$Equi7t3 <- names(cov.c)[i]
  out$row.index <- row.index
  xy <- data.frame(pnts@coords)
  names(xy) <- c("X","Y")
  out <- cbind(out, xy)
  return(out)
}
