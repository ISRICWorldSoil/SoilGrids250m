## predict 2D variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_2D <- function(i, gm1, gm2, varn, in.path, out.path, zmin, zmax, gm1.w=NULL, gm2.w=NULL, mask_value){  
  out.c <- paste0(out.path, "/", i, "/", varn, "_M_", i, ".tif")
  if(!file.exists(out.c)){
    rds.file = paste0(in.path, "/", i, "/", i, ".rds")
    if(file.exists(rds.file)&file.size(rds.file)>1e3){ 
      #message("Predicting from .rds file...")
      newdata <- readRDS(rds.file)
      ## some data.frames have a single pixel / should be skipped:
      if(nrow(newdata@data)>1){
        ## filter any missing layers:
        sel.na <- colSums(sapply(newdata@data, function(x){!is.na(x)}))==0
        if(any(sel.na==TRUE)){
          sel.na <- attr(sel.na, "names")[which(sel.na)]
          for(i in 1:length(sel.na)){
            newdata@data[,sel.na[i]] = mask_value[[sel.na[i]]]
          }
        }
        if(is.character(gm1)){
          ## Load RandomForest model:
          gm1 <- readRDS(gm1)
        }
        if(is.null(gm1.w)){ gm1.w = 1/gm1$prediction.error }
        if(is.character(gm2)){
          gm2 <- readRDS(gm2)
        }
        if(is.null(gm2.w)){ gm2.w = 1/(min(gm2$results$RMSE, na.rm=TRUE)^2) }
        v1 <- predict(gm1, newdata@data, na.action=na.pass)$predictions
        v2 <- predict(gm2, newdata@data, na.action=na.pass)
        v <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
        if(varn=="BDRLOG"){ v <- v * 100 }
        #if(varn=="logBDTICM"){ v <- expm1(v) }
        if(is.na(zmax)){ zmax = Inf }
        newdata$v <- ifelse(v < zmin, zmin, ifelse(v > zmax, zmax, v))
        if(varn %in% c("BDRLOG","BDRICM")){
          writeGDAL(newdata["v"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
        } else {
          writeGDAL(newdata["v"], out.c, type="Int32", mvFlag=-32768, options="COMPRESS=DEFLATE") ## gsub("logBDTICM", "BDTICM", out.c)
        }
        gc()
      }
    }
  }
}