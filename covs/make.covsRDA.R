## create RDA files per tile:

make.covsRDA <- function(in.path, i){
  out.rda <- paste0(in.path, "/", i, "/", i, ".rda")
  if(!file.exists(out.rda)){
    cov.lst <- list.files(path=paste(in.path, i, sep="/"), glob2rx("*_*_*_*.tif$"), full.names=TRUE)
    mask <- grep("LMK", cov.lst)
    m <- readGDAL(cov.lst[mask], silent=TRUE)
    m <- as(m, "SpatialPixelsDataFrame")
    m$LATWGS84 <- spTransform(m, CRS("+proj=longlat +datum=WGS84"))@coords[,2]
    for(k in 1:length(cov.lst)){
      if(!k==mask){
        d <- readGDAL(cov.lst[k], silent=TRUE)$band1[m@grid.index]
        ## They should all have complete values, but we run one more check just to be certain:
        dn <- which(is.na(d))
        if(length(dn)>0){ d[dn] <- quantile(d, probs=0.5, na.rm=TRUE) }
        m@data[,strsplit(basename(cov.lst[k]), "_")[[1]][1]] <- d 
      }
    }
    gc()
    save(m, file=out.rda, compress=TRUE)
  }
}