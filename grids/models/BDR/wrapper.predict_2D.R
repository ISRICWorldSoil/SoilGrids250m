## predict 2D variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_2D <- function(i, gm_path1, gm_path2, varn, in.path, out.path, z.min, z.max, gm1.w, gm2.w){ ## , nthreads 
  out.all <- as.vector(sapply(varn, function(x){paste0(out.path, "/", i, "/", x, "_M_", i, ".tif")}))
  if(any(!file.exists(out.all))){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    gz.file <- paste0(in.path, "/", i, "/", i, ".csv.gz")
    if(file.exists(gz.file)){
      message("Predicting from csv.gz file.")
      m.grid <- h2o.importFile(path = gz.file)
    } else {
      m.grid <- as.h2o(m@data, destination_frame="m.grid")
    }
    for(x in 1:length(varn)){
      out.c <- paste0(out.path, "/", i, "/", varn[x], "_M_", i, ".tif")
      gm1 <- h2o.loadModel(paste(gm_path1[[varn[x]]]))
      v1 <- as.data.frame(h2o.predict(gm1, m.grid, na.action=na.pass))$predict
      if(missing(gm1.w)){ gm1.w = gm1@model$training_metrics@metrics$r2 }
      gm2 <- h2o.loadModel(paste(gm_path2[[varn[x]]]))
      if(missing(gm2.w)){ gm2.w = gm2@model$training_metrics@metrics$r2 }
      v2 <- as.data.frame(h2o.predict(gm2, m.grid, na.action=na.pass))$predict
      v <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
      gc()
      if(varn[x]=="BDRLOG"){ v <- v * 100 }
      #if(varn[x]=="logBDTICM"){ v <- expm1(v) }
      #if(varn[x]=="BDRICM"){ v <- ifelse(v>200, 200, v) }
      if(is.na(z.max[x])){ z.max[x] = Inf }
      m$v <- ifelse(v < z.min[x], z.min[x], ifelse(v > z.max[x], z.max[x], v))
      if(varn[x] %in% c("BDRLOG","BDRICM")){
        writeGDAL(m["v"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      } else {
        writeGDAL(m["v"], out.c, type="Int32", mvFlag=-32768, options="COMPRESS=DEFLATE") ## gsub("logBDTICM", "BDTICM", out.c)
      }
      gc()
    }
    gc()
    x = h2o.removeAll()
  }
}