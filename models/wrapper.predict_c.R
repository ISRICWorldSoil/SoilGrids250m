## predict soil types using a model "mg" and write GeoTifs out
## by: Tom.Hengl@isric.org

wrapper.predict_c <- function(i, gm, varn, in.path, out.path){
  cov.lst <- list.files(path=paste(in.path, i, sep="/"), glob2rx("*_*_*_*.tif$"), full.names=TRUE)
  tax <- gm$lev
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
  out <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out)){  
    ## most probable class:
    try( c <- predict(gm, m, na.action = na.pass) )
    if(!class(.Last.value)[1]=="try-error"){
      gc()
      m@data[,paste0(varn, ".i")] <- as.integer(c)
      writeGDAL(m[paste0(varn, ".i")], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
    }
    ## now all classes:
    try( probs <- predict(gm, m, type="probs", na.action = na.pass) )
    if(!class(.Last.value)[1]=="try-error"){
      m@data <- data.frame(probs)
      gc()
      for(j in 1:ncol(m)){
        m$v <- round(m@data[,j]*100)
        out <- paste0(out.path, "/", i, "/", varn, "_", tax[j], "_", i, ".tif")
        writeGDAL(m["v"], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      }
    }
  }
}