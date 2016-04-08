## predict numeric variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_nCSV <- function(i, gm_path, sd=c(2.5, 7.5, 22.5, 45, 80, 150), varn, in.path, out.path, z.min, z.max){  
  out.all <- as.vector(sapply(1:length(sd), function(x){paste0(out.path, "/", i, "/", varn, "_M_sd", x, "_", i, ".tif")}))
  v <- list(NULL)
  if(any(!file.exists(out.all))){
    gz.file <- paste0(in.path, "/", i, "/", i, ".csv.gz")
    message("Reading gz file...")
    m.grid <- h2o.importFile(localH2O, path = gz.file)
    ## load all models:
    gm <- lapply(gm_path, function(x){h2o.loadModel(x, h2o.getConnection())})
    for(j in 1:length(sd)){
      m.depth <- as.h2o(data.frame(DEPTH=rep(sd[j], m.grid@mutable$nrows)), conn=h2o.getConnection(), destination_frame="m.depth")
      newdata <- h2o.cbind(m.grid, m.depth)
      for(x in 1:length(varn)){
        ## predict per property:
        out.c <- paste0(out.path, "/", i, "/", varn[x], "_M_sd", 1:length(sd), "_", i, ".tif")
        v1 <- as.data.frame(h2o.predict(gm[[x]], newdata, na.action=na.pass))$predict
        if(varn[x]=="PHIHOX"|varn[x]=="PHIKCL"){ v1 <- v1 * 10 }
        v[[((j-1)*length(sd)+x)]] <- ifelse(v1 < z.min[x], z.min[x], ifelse(v1 > z.max[x], z.max[x], v1))
      }
    }
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    m <- as(m, "SpatialPixels")
    gc()
    message("Writting predictions to GeoTifs...")
    sfInit(parallel=TRUE, cpus=40)
    sfExport("v", "m", "out.all")
    sfLibrary(rgdal)
    sfLibrary(sp)
    x <- sfClusterApplyLB(1:length(v), fun=function(x){ writeGDAL(SpatialPixelsDataFrame(m, data=data.frame(v[[x]])), out.all[x], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE") })
    sfStop()
  }
}