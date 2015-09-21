## predict soil types using a model "mg" and write GeoTifs out
## by: Tom.Hengl@isric.org

wrapper.predict_c <- function(i, gm1, gm2, varn, in.path, out.path, col.legend){
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out.c)){
    load(paste0(in.path, "/", i, "/", i, ".rda"))
    ## predict probabilities:
    probs1 <- predict(gm1, m, type="probs", na.action = na.pass) ## nnet
    probs2 <- predict(gm2, m, type="prob", na.action = na.pass) ## randomForest
    lt <- list(probs1, probs2)
    ## ensemble predictions (simple average):
    probs <- Reduce("+", lt) / length(lt)
    m@data <- data.frame(probs)
    gc()
    tax <- names(m)
    for(j in 1:ncol(m)){
      out <- paste0(out.path, "/", i, "/", varn, "_", tax[j], "_", i, ".tif")
      if(!file.exists(out)){
        m$v <- round(m@data[,j]*100)
        writeGDAL(m["v"], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      }
    }
    ## most probable class:
    col.tbl <- join(col.legend, data.frame(Group=tax, int=1:length(tax))) 
    m$cl <- col.tbl[match(apply(probs,1,which.max), col.tbl$int),"Number"]   
    writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group)))
    ## TH: does not work 'colorTable=list(col.tbl$COLOR)'
  }
}