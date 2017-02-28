## Overlay and extract values using a tiling system:

extract.tiled <- function(x, tile.pol, path="/data/tt/SoilGrids250m/predicted250m", ID="ID", cpus=56){
  x$row.index <- 1:nrow(x)
  ov.c <- over(spTransform(x, CRS(proj4string(tile.pol))), tile.pol)
  ov.t <- which(!is.na(ov.c[,ID]))
  ## for each point get the tile name:
  ov.c <- data.frame(ID=ov.c[ov.t,ID], row.index=ov.t)
  tiles.lst <- basename(dirname(list.files(path=path, pattern=glob2rx("*.rds$"), recursive=TRUE)))
  ov.c <- ov.c[ov.c[,ID] %in% sapply(tiles.lst, function(i){strsplit(i, "T")[[1]][2]}),]
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
  out <- dplyr::bind_rows(ov.lst)
  out <- plyr::join(x@data, as.data.frame(out), type="left", by="row.index")
  return(out)
}

.extract.tile <- function(i, x, ID, path, ov.c, cov.c, tile.pol){
  row.index <- ov.c$row.index[ov.c[,ID]==names(cov.c)[i]]
  pnts <- x[row.index,]
  pnts <- spTransform(pnts, CRS(proj4string(tile.pol)))
  m <- readRDS(paste0(path, "/T", names(cov.c)[i], "/T", names(cov.c)[i], ".rds"))
  out <- sp::over(y=m, x=pnts)
  out$band1 <- NULL
  out$row.index <- row.index
  xy <- data.frame(pnts@coords)
  names(xy) <- c("X","Y")
  out <- cbind(out, xy)
  return(out)
}