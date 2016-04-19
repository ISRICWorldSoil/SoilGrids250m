## predict numeric variables using models "gm" and write GeoTifs out - SoilGrids250m (R objects only)
## by: Tom.Hengl@isric.org

## run per property:
wrapper.predict_np <- function(i, gm1, gm2, sd=c(0, 5, 15, 30, 60, 100, 200), varn, in.path, out.path, zmin, zmax, gm1.w=NULL, gm2.w=NULL, type="Int16", mvFlag=-32768){
  out.all <- as.vector(sapply(varn, function(x){paste0(out.path, "/", i, "/", x, "_M_sl", 1:length(sd), "_", i, ".tif")}))
  if(any(!file.exists(out.all))){
    rds.file = paste0(in.path, "/", i, "/", i, ".rds")
    if(file.exists(rds.file)&file.size(rds.file)>1e3){ 
      #message("Predicting from .rds file...")
      newdata <- readRDS(rds.file)
      ## some data.frames have a single pixel / should be skipped:
      if(nrow(newdata@data)>1){
        if(is.character(gm1)){
          ## Load RandomForest model:
          gm1 <- readRDS(gm1)
        }
        if(is.null(gm1.w)){ gm1.w = 1/gm1$prediction.error }
        if(is.character(gm2)){
          ## Load xgboost model from binary file?
          #gm2 <- xgb.load(gm2[[varn[x]]])
          gm2 <- readRDS(gm2)
        }
        if(is.null(gm2.w)){ gm2.w = 1/(min(gm2$results$RMSE, na.rm=TRUE)^2) }
        ## Predict at 7 depths so that after predictions we can average to get block estimates st=c(2.5, 7.5, 22.5, 45, 80, 150)
        for(j in 1:length(sd)){
          newdata$DEPTH.f = sd[j]
          v1 <- predict(gm1, newdata@data, na.action=na.pass)$predictions #  , num.threads=1 ## Default for number of threads
          v2 <- predict(gm2, newdata@data, na.action=na.pass)
          v <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
          ## For pH -> convert to integers:
          if(varn=="PHIHOX"|varn=="PHIKCL"){ v <- v * 10 }
          if(varn %in% c("P", "S", "B", "Cu", "Zn")){ v <- v * 100 }
          newdata$v <- ifelse(v < zmin, zmin, ifelse(v > zmax, zmax, v))
          writeGDAL(newdata["v"], out.all[j], type=type, mvFlag=mvFlag, options="COMPRESS=DEFLATE")   
        }
      }
    }
    gc()
  }
}
