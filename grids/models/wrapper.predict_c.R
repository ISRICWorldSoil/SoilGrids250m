## predict soil types using a model "mg" and write GeoTifs out
## by: Tom.Hengl@isric.org

wrapper.predict_c <- function(i, gm1, gm2, varn, in.path, out.path, col.legend, check.names=TRUE){ #gm3,
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out.c)){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    ## Reading models from disk does not help with RAM
    #gm1 <- readRDS(gm1)
    #gm2 <- readRDS(gm2)
    ## predict probabilities:
    probs1 <- predict(gm1, m, type="probs", na.action = na.pass) ## nnet
    probs2 <- predict(gm2, m, type="prob", na.action = na.pass) ## randomForest
    ## this takes a lot of RAM - can it be reduced?
    #probs3 <- attr(predict(gm3, m, probability=TRUE, na.action = na.pass), "probabilities") ## SVM
    lt <- list(probs1[,gm1$lev], probs2[,gm1$lev])
    ## ensemble predictions (simple average):
    probs <- Reduce("+", lt) / length(lt)
    m@data <- data.frame(probs)
    cm <- colSums(m@data, na.rm=TRUE)
    if(all(!cm==0)){
      tax <- names(m)
      for(j in 1:ncol(m)){
        out <- paste0(out.path, "/", i, "/", varn, "_", tax[j], "_", i, ".tif")
        if(!file.exists(out)){
          m$v <- round(m@data[,j]*100)
          writeGDAL(m["v"], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
        }
      }
      ## most probable class:
      if(check.names==TRUE){ 
        Group = gsub("\\.", " ", gsub("\\.$", "\\)", gsub("\\.\\.", " \\(", tax))) 
      } else {
        Group = tax
      }
      col.tbl <- join(col.legend, data.frame(Group=Group, int=1:length(tax)))
      ## match most probable class (takes 1-2 mins):
      m$cl <- col.tbl[match(apply(probs,1,which.max), col.tbl$int),"Number"]   
      writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group)))  ## TH: 'colorTable=list(col.tbl$COLOR)' does not work
    }
    gc()
    gc()
  }
}