## predict numeric variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_n <- function(i, gm_path, sd=c(2.5, 7.5, 22.5, 45, 80, 150), varn, in.path, out.path, z.min, z.max){ 
  localH2O = h2o.init(nthreads = 2) ## max_mem_size = "20g"
  m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
  m.grid <- as.h2o(m@data, conn=h2o.getConnection(), destination_frame="m.grid")
  for(x in 1:length(varn)){
    gm <- h2o.loadModel(gm_path[[x]], h2o.getConnection())
    out.c <- paste0(out.path, "/", i, "/", varn[x], "_M_sd", 1:length(sd), "_", i, ".tif")
    if(any(!file.exists(out.c))){
      for(j in 1:length(sd)){
        ## predict:
        m.depth <- as.h2o(data.frame(DEPTH=rep(sd[j], nrow(m@data))), conn=h2o.getConnection(), destination_frame="m.depth")
        newdata <- h2o.cbind(m.grid, m.depth)
        v1 <- as.data.frame(h2o.predict(gm,  newdata, na.action=na.pass))$predict
        if(varn[x]=="PHIHOX"|varn[x]=="PHIKCL"){ v1 <- v1 * 10 }
        m$v <- ifelse(v1 < z.min[x], z.min[x], v1)
        m$v <- ifelse(m$v > z.max[x], z.max[x], m$v)
        writeGDAL(m["v"], out.c[j], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
        gc()
      }
    }    
  }
  h2o.rm("m.grid", localH2O)
  remove(m.grid)
  h2o.shutdown(conn = h2o.getConnection(), prompt = FALSE)
  gc()
  gc()
}