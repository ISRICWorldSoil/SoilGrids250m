## predict numeric variables using models "gm" and write GeoTifs out - SoilGrids250m (h2o objects)
## Sub-optimal / needs revision
## by: Tom.Hengl@isric.org

wrapper.predict_n <- function(i, gm_path1, gm_path2, sd=c(0, 5, 15, 30, 60, 100, 200), varn, in.path, out.path, z.min, z.max, gm1.w=NULL, gm2.w=NULL, type="Int16", mvFlag=-32768, cpus=9){  
  out.all <- as.vector(sapply(varn, function(x){paste0(out.path, "/", i, "/", x, "_M_sd", 1:length(sd), "_", i, ".tif")}))
  ## Predict at 7 depths so that after predictions we can average to get block estimates
  #sd=c(2.5, 7.5, 22.5, 45, 80, 150)
  if(any(!file.exists(out.all))){
    ## Get spatial structure:
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    m <- m[1]
    gc()
    ## load covariates:
    gz.file <- paste0(in.path, "/", i, "/", i, ".csv.gz")
    message("Predicting from csv.gz file...")
    m.grid <- h2o.importFile(path = gz.file)
    for(j in 1:length(sd)){
      m.depth <- as.h2o(data.frame(DEPTH.f=rep(sd[j], nrow(m@data))), destination_frame="m.depth")
      newdata <- h2o.cbind(m.grid, m.depth)
      out.c <- paste0(out.path, "/", i, "/", varn, "_M_sd", j, "_", i, ".tif")
      if(any(!file.exists(out.c))){
        x = lapply(1:length(varn), function(x){predict_vars(x, varn, m, newdata, gm_path1, gm_path2, z.min, z.max, out.c, gm1.w, gm2.w, type, mvFlag)})
      }
    }
    h2o.removeAll()
    gc()
  }
}

predict_vars <- function(x, varn, m, newdata, gm_path1, gm_path2, z.min, z.max, out.c, gm1.w, gm2.w, type, mvFlag){
  ## predict all target variables for a fixed depth:
  if(!file.exists(out.c[x])){
    gm1 <- h2o::h2o.loadModel(gm_path1[[varn[x]]])
    if(is.null(gm1.w)){ gm1.w = gm1@model$training_metrics@metrics$r2 }
    gm2 <- h2o::h2o.loadModel(gm_path2[[varn[x]]])
    if(is.null(gm2.w)){ gm2.w = gm2@model$training_metrics@metrics$r2 }
    v1 <- as.data.frame(h2o::h2o.predict(gm1, newdata, na.action=na.pass))$predict
    v2 <- as.data.frame(h2o::h2o.predict(gm2, newdata, na.action=na.pass))$predict
    v <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
    ## For pH -> convert to integers:
    if(varn[x]=="PHIHOX"|varn[x]=="PHIKCL"){ v <- v * 10 }
    if(varn[x] %in% c("P", "S", "B", "Cu", "Zn")){ v <- v * 100 }
    m$v <- ifelse(v < z.min[x], z.min[x], ifelse(v > z.max[x], z.max[x], v))
    writeGDAL(m["v"], out.c[x], type=type, mvFlag=mvFlag, options="COMPRESS=DEFLATE") 
  }
}