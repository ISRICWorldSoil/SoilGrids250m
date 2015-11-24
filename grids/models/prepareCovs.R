## function to convert mosaics to EQUI7 tiles:

prepareCovsSoilGrids250m <- function(nms=c("DEM","MOD13Q1","MCD43A4","MOD11A2","MOD10A2","PREm","EcoTapestry","floods","NEO","Other"), s.zone, s.lst, close.gap=TRUE){
  
  if(missing(s.zone)){ s.zone <- 1:7 }
  
  ## some close gaps runs in SAGA GIS result in infinite loops, hence limit to 120secs
  #evalWithTimeout( {
  if(any(nms %in% "DEM")){    
    for(k in 1:length(DEM.out.lst)){
      for(i in s.zone){ ## skip antartica??
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        if(k==1){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="near", ot="Int16", nodata="-32767", zmin=-350, close.gap=close.gap)})
        }
        if(k==2|k==4|k==5){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="near", ot="Int16", nodata="-32767", Z_SCALE=100, close.gap=close.gap)})
        }
        if(k==3){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="near", ot="Int16", nodata="-32767", Z_SCALE=10000, close.gap=close.gap)})
        }
        if(k==6){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="near", ot="Int16", nodata="-32767", Z_SCALE=10, close.gap=close.gap)})
        }
        if(k==7|k==8){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="near", ot="Int16", nodata="-32767", Z_SCALE=1000, close.gap=close.gap)})
        }
        if(k==9){
          sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=DEM.in.lst[[k]][i], nm=DEM.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32767", Z_SCALE=100, zmin=0, close.gap=close.gap)})
        }
      }
    }
  }
  
  if(any(nms %in% "MOD13Q1")){
    for(k in 1:length(EVIM.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=EVIM.lst[k], nm=EVIM.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=-2000, close.gap=close.gap)})
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=EVIS.lst[k], nm=EVIS.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
      }
    }
  }
  
  if(any(nms %in% "MCD43A4")){  
    for(k in 1:length(NBR4.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=NBR4.lst[k], nm=NBR4.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=NBR7.lst[k], nm=NBR7.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
      }
    }
  }
  
  if(any(nms %in% "MOD11A2")){  
    for(k in 1:length(LSTD.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LSTD.lst[k], nm=LSTD.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LSTN.lst[k], nm=LSTN.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
        }
      }
      for(k in 1:length(LSSD.out.lst)){
        for(i in s.zone){
          if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LSSD.lst[k], nm=LSSD.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32767", srcnodata="65535", dstnodata="-32767", zmin=0, close.gap=close.gap)})
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LSSN.lst[k], nm=LSSN.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32767", srcnodata="-32768", zmin=0, close.gap=close.gap)})
      }
    }
  }
  
  if(any(nms %in% "MOD10A2")){
    for(k in 1:length(SNW.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=SNW.lst[k], nm=SNW.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32767", zmin=0, close.gap=FALSE)})
      }
    }
  }
  
  if(any(nms %in% "PREm")){
    for(k in 1:length(PREm.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=PREm.lst[k], nm=PREm.out.lst[k], resample.type="bilinear", ot="Int16", nodata="-32768", Z_SCALE=0.05, srcnodata=-32767, dstnodata=-32767, zmin=0, close.gap=close.gap)})
      }
    }
  }
  
  if(any(nms %in% "EcoTapestry")){
    for(k in 1:length(LIT.out.lst)){
      for(i in s.zone){ 
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LIT.lst[[k]][i], nm=LIT.out.lst[k], resample.type="near", ot="Byte", nodata="255", srcnodata=255, dstnodata=255, zmin=0, close.gap=close.gap)})
      }
    }
    for(k in 1:length(LFO.out.lst)){
      for(i in s.zone){ 
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=LFO.lst[[k]][i], nm=LFO.out.lst[k], resample.type="near", ot="Byte", nodata="255", srcnodata=255, dstnodata=255, zmin=0, close.gap=close.gap)})
      }
    }
  }
    
  if(any(nms %in% "NEO")){
    for(k in 1:length(NEO.out.lst)){
      for(i in s.zone){
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=NEO.lst[k], nm=NEO.out.lst[k], resample.type="cubicspline", ot="Int16", nodata="-32768", srcnodata=-32767, dstnodata=-32767, zmin=0, close.gap=close.gap)})
      }
    }
  }

  if(any(nms %in% "floods")){    
    for(k in 1:length(FW.out.lst)){
      for(i in s.zone){ ## skip antartica??
        if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
        sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input=FW.lst[[k]][i], nm=FW.out.lst[k], resample.type="near", ot="Byte", nodata="255", srcnodata=255, zmin=0, close.gap=close.gap)})
      }
    }
  }
    
  if(any(nms %in% "Other")){
    for(i in s.zone){ 
      if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
      sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input="/data/Groundwater/World_wtd_v2_f.tif", nm="GTDHYS3", resample.type="cubicspline", ot="Int16", nodata="-32767", srcnodata=-32768, zmin=-1, close.gap=close.gap)})
    }
    for(i in s.zone){ 
      if(missing(s.lst)){ s.lst <- 1:length(equi7t3[[i]]) }
      sapply(s.lst, function(x){tile.tif(t=equi7t3[[i]][x,], input="/data/mangroves/MANGPRf_500m.sdat", nm="MNGUSG", resample.type="bilinear", ot="Byte", nodata="255", srcnodata=255, zmin=0, close.gap=close.gap)})
    }
  }
  #} , timeout=120, onTimeout="silent")
  
}