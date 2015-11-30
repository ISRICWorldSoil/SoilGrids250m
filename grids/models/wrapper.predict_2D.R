## predict 2D variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_2D <- function(i, gm_path, varn, in.path, out.path, z.min, z.max){ ## , nthreads 
  out.all <- as.vector(sapply(varn, function(x){paste0(out.path, "/", i, "/", x, "_M_", i, ".tif")}))
  if(any(!file.exists(out.all))){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    gz.file <- paste0(in.path, "/", i, "/", i, ".csv.gz")
    if(file.exists(gz.file)){
      message("Predicting from csv.gz file.")
      m.grid <- h2o.importFile(localH2O, path = gz.file)
    } else {
      m.grid <- as.h2o(m@data, conn=h2o.getConnection(), destination_frame="m.grid")
    }
    for(x in 1:length(varn)){
      gm <- h2o.loadModel(gm_path[[x]], h2o.getConnection())
      out.c <- paste0(out.path, "/", i, "/", varn[x], "_M_", i, ".tif")
      v1 <- as.data.frame(h2o.predict(gm, m.grid, na.action=na.pass))$predict
      gc()
      if(varn[x]=="BDRLOG"){ v1 <- v1 * 100 }
      if(is.na(z.max[x])){ z.max[x] = Inf }
      m$v <- ifelse(v1 < z.min[x], z.min[x], ifelse(v1 > z.max[x], z.max[x], v1))
      #plot(raster(m["v"]), col=SAGA_pal[[1]], zlim=c(0,100))
      if(varn[x]=="BDRLOG"){
        writeGDAL(m["v"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      } else {
        writeGDAL(m["v"], out.c, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")  
      }
      gc()
    }
  }
}