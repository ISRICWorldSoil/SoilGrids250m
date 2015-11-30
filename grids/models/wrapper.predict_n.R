## predict numeric variables using models "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_n <- function(i, gm_path1, gm_path2, sd=c(2.5, 7.5, 22.5, 45, 80, 150), varn, in.path, out.path, z.min, z.max, gm1.w, gm2.w){ ## , nthreads 
  out.all <- as.vector(sapply(varn, function(x){paste0(out.path, "/", i, "/", x, "_M_sd", 1:length(sd), "_", i, ".tif")}))
  if(any(!file.exists(out.all))){
    ## Get spatial structure:
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    m <- m[1]
    gc()
    ## load covariates:
    gz.file <- paste0(in.path, "/", i, "/", i, ".csv.gz")
    message("Predicting from csv.gz file...")
    m.grid <- h2o.importFile(localH2O, path = gz.file)
    for(j in 1:length(sd)){
      m.depth <- as.h2o(data.frame(DEPTH=rep(sd[j], nrow(m@data))), conn=h2o.getConnection(), destination_frame="m.depth")
      newdata <- h2o.cbind(m.grid, m.depth)
      out.c <- paste0(out.path, "/", i, "/", varn, "_M_sd", j, "_", i, ".tif")
      if(any(!file.exists(out.c))){
        for(x in 1:length(varn)){
          ## predict:
          gm1 <- h2o.loadModel(gm_path1[[x]], h2o.getConnection())
          if(missing(gm1.w)){ gm1.w = gm1@model$training_metrics@metrics$r2 }
          gm2 <- h2o.loadModel(gm_path2[[x]], h2o.getConnection())
          if(missing(gm2.w)){ gm2.w = gm2@model$training_metrics@metrics$r2 }
          v1 <- as.data.frame(h2o.predict(gm1, newdata, na.action=na.pass))$predict
          v2 <- as.data.frame(h2o.predict(gm2, newdata, na.action=na.pass))$predict
          v <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
          ## For pH -> convert to integers:
          if(varn[x]=="PHIHOX"|varn[x]=="PHIKCL"){ v <- v * 10 }
          m$v <- ifelse(v < z.min[x], z.min[x], ifelse(v > z.max[x], z.max[x], v))
          writeGDAL(m["v"], out.c[x], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
        }
      }
    }
    gc()
    h2o.removeAll(localH2O)
  }
}