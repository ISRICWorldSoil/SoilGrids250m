## function to create RDA files per tile
## Tom.Hengl@isric.org

make.covsRDA <- function(in.path, i){
  out.rda <- paste0(in.path, "/", i, "/", i, ".rds")
  if(!file.exists(out.rda)){
    cov.lst <- list.files(path=paste(in.path, i, sep="/"), glob2rx("*_*_*_*.tif$"), full.names=TRUE)
    mask <- grep("LMK", cov.lst)
    m <- readGDAL(cov.lst[mask], silent=TRUE)
    m <- as(m, "SpatialPixelsDataFrame")
    m$LATWGS84 <- spTransform(m, CRS("+proj=longlat +datum=WGS84"))@coords[,2]
    for(k in 1:length(cov.lst)){
      if(!k==mask){
        covn = strsplit(basename(cov.lst[k]), "_")[[1]][1]
        d <- readGDAL(cov.lst[k], silent=TRUE)$band1[m@grid.index]
        dn <- which(is.na(d)|d<= -32767|d>= 65535)
        if(length(dn)>0){ 
          d[dn] <- quantile(d, probs=0.5, na.rm=TRUE)
          ## replace all missing values with median value
        }
        m@data[,covn] <- d 
      }
    }
    gc()
    saveRDS(m, file=out.rda, compress=TRUE)
  }
}