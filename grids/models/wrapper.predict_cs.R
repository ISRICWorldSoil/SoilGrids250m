## Split randomForest models

split_rf <- function(rf, num_splits=4){
  ## Split up tree indices
  splits <- split(1:rf$num.trees, cut(1:rf$num.trees, num_splits, labels = FALSE))
  ## Split up trees
  rfs <- list(NULL)
  for(i in 1:length(splits)){
    rfs[[i]] <- rf
    rfs[[i]]$num.trees <- length(splits[[i]])
    rfs[[i]]$forest$num.trees <- length(splits[[i]])
    rfs[[i]]$forest$child.nodeIDs <- rf$forest$child.nodeIDs[splits[[i]]]
    rfs[[i]]$forest$split.varIDs <- rf$forest$split.varIDs[splits[[i]]]
    rfs[[i]]$forest$split.values <- rf$forest$split.values[splits[[i]]]
    rfs[[i]]$forest$terminal.class.counts <- rf$forest$terminal.class.counts[splits[[i]]]
  }
  return(rfs)
}

split_predict_c <- function(i, gm, in.path, out.path, split_no, varn){
  rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", split_no, ".rds")
  if(!file.exists(rds.out)){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    x = predict(gm, m@data, probability=TRUE, na.action = na.pass, num.threads=1)$predictions
    saveRDS(x, file=rds.out)
  }
}

predict_nnet <- function(i, gm, in.path, out.path, varn){
  rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds")
  if(!file.exists(rds.out)){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    x = predict(gm, m@data, type="prob", na.action = na.pass)
    saveRDS(x, file=rds.out)
  }
}

sum_predictions <- function(i, in.path, out.path, varn, gm1.w, gm2.w, col.legend, soil.fix, lev, num_splits=4, check.names=FALSE){
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out.c)){
    ## load all objects:
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    if(nrow(m@data)>1){
      mfix <- m$BICUSG5
      lfix <- levels(as.factor(paste(mfix)))
      rf.ls = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", 1:num_splits, ".rds")
      probs2 <- lapply(rf.ls, readRDS)
      probs2 <- data.frame(Reduce("+", probs2) / num_splits)
      probs1 <- readRDS(paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds"))
      ## weighted average:
      probs <- list(probs1[,lev]*gm1.w, probs2[,lev]*gm2.w)
      probs <- data.frame(Reduce("+", probs) / rep(gm1.w+gm2.w, length(probs1)))
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
      unlink(rf.ls)
      unlink(paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds"))
      gc()
    }  
  }
}
