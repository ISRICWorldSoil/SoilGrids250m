## predict soil types using a model "mg" and write GeoTifs out
## by: Tom.Hengl@isric.org

wrapper.predict_c <- function(i, gm1, gm2, varn, in.path, out.path, col.legend, check.names=TRUE, soil.fix, gm1.w, gm2.w){ 
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out.c)){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    if(nrow(m@data)>1){
      mfix <- m$BICUSG5
      lfix <- levels(as.factor(paste(mfix)))
      probs <- predict_df(gm1, gm2, m@data, gm1.w, gm2.w)
      gc()
      if(any(unique(unlist(soil.fix)) %in% lfix)){
        ## correct probabilities using soil-climate matrix:
        for(k in lfix){
          sel = sapply(soil.fix, function(i){i[grep(k, i)]})
          sel <- sel[sapply(sel, function(i){length(i)>0})]
          if(length(sel)>0){ probs[mfix==k,names(sel)] <- 0 }
        }
        ## Standardize values so they sums up to 100:
        rs <- rowSums(probs, na.rm=TRUE)
        m@data <- data.frame(lapply(probs, function(i){i/rs}))
      } else {
        rs <- rowSums(probs, na.rm=TRUE)
        m@data <- probs
      }
      if(sum(rs,na.rm=TRUE)>0&length(rs)>0){
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
        col.tbl <- plyr::join(data.frame(Group=Group, int=1:length(tax)), col.legend, type="left")
        ## match most probable class (takes 1-2 mins):
        m$cl <- col.tbl[match(apply(m@data,1,which.max), col.tbl$int),"Number"]  
        writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group)))  ## TH: 'colorTable=list(col.tbl$COLOR)' does not work
      }
      gc()
    }
  }
}


predict_df <- function(gm1, gm2, newdata, gm1.w=NULL, gm2.w=NULL){
  ## read models from hard disk - saves RAM
  if(is.character(gm1)){ gm1 <- readRDS(file=gm1) }
  if(is.character(gm2)){ gm2 <- readRDS(file=gm2) }
  ## predict probabilities
  ## most computational demanding / often takes a lot of RAM
  probs1 <- predict(gm1, newdata, type="prob", na.action = na.pass) ## nnet
  #probs2 <- predict(gm2, newdata, type="prob", na.action = na.pass) ## RandomForest
  probs2 <- predict(gm2, newdata, probability=TRUE, na.action = na.pass)$predictions ## ranger
  ## ensemble predictions:
  if(is.null(gm1.w)&is.null(gm2.w)){
    ## simple average:
    lt <- list(probs1[,gm1$lev], probs2[,gm1$lev])
    probs <- data.frame(Reduce("+", lt) / length(lt))
  } else {
    ## weighted average:
    lt <- list(probs1[,gm1$lev]*gm1.w, probs2[,gm1$lev]*gm2.w)
    probs <- data.frame(Reduce("+", lt) / rep(gm1.w+gm2.w, length(lt)))
  }
  return(probs)
}
