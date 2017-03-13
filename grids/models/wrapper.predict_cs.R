## Split randomForest models and make predictions for factor and numeric type variables
## By tom.hengl@isric.org / suggestion to split RF models by Marvin N. Wright <wright at imbs.uni-luebeck.de>

## -------------------------------
## Factor-type variables
## -------------------------------

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

split_predict_c <- function(i, gm, in.path, out.path, split_no, varn, num.threads=1){
  if(dir.exists(out.path)){
    if(is.null(num_splits)){ 
      rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_rf.rds")
    } else {
      rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", split_no, ".rds")
    }
    if(any(c(!file.exists(rds.out),file.size(rds.out)==0))){
      m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
      ## round up numbers otherwise too large objects
      x = round(predict(gm, m@data, probability=TRUE, na.action = na.pass, num.threads=num.threads)$predictions*100)
      saveRDS(x, file=rds.out)
    } 
  }
}

predict_nnet <- function(i, gm, in.path, out.path, varn){
  rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds")
  if(any(c(!file.exists(rds.out),file.size(rds.out)==0))){
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    x = round(predict(gm, m@data, type="prob", na.action = na.pass)*100)
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
      if(is.null(num_splits)){
        probs2 <- readRDS(paste0(out.path, "/", i, "/", varn,"_", i, "_rf.rds"))
      } else {
        rf.ls = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", 1:num_splits, ".rds")
        probs2 <- lapply(rf.ls, readRDS)
        probs2 <- Reduce("+", probs2) / num_splits
      }
      probs1 <- readRDS(paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds"))
      ## weighted average:
      probs <- list(probs1[,lev]*gm1.w, probs2[,lev]*gm2.w)
      probs <- data.frame(Reduce("+", probs) / rep(gm1.w+gm2.w, length(probs1)))
      if(any(unique(unlist(soil.fix)) %in% lfix)){
        ## correct probabilities using soil-climate matrix:
        for(k in lfix){
          sel = sapply(soil.fix, function(i){i[which(i==k)]})
          sel <- sel[sapply(sel, function(i){length(i)>0})]
          if(length(sel)>0){ probs[mfix==k,names(sel)] <- 0 }
        }
      }
      ## Standardize values so they sums up to 100:
      rs <- rowSums(probs, na.rm=TRUE)
      m@data <- data.frame(lapply(probs, function(i){i/rs}))
      ## Write GeoTiffs:
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
        ## match most probable class
        m$cl <- col.tbl[match(apply(m@data,1,which.max), col.tbl$int),"Number"]  
        writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group)))  ## TH: 'colorTable=list(col.tbl$COLOR)' does not work
      }
      unlink(rf.ls)
      unlink(paste0(out.path, "/", i, "/", varn,"_", i, "_nnet.rds"))
      gc()
    }  
  }
}


## predict using ranger model only (model splitting not needed):
factor_predict_ranger <- function(i, gm, in.path, out.path, varn, col.legend, soil.fix, check.names=FALSE){
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(!file.exists(out.c)){
    ## load all objects:
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    #if(any(names(m) %in% "OCCGSW7.tif")){ m$OCCGSW7.tif = ifelse(m$OCCGSW7.tif>100, 0, m$OCCGSW7.tif) }
    if(nrow(m@data)>1){
      mfix <- paste(m$BICUSG5.tif)
      lfix <- levels(as.factor(mfix))
      probs = data.frame(round(predict(gm, m@data, probability=TRUE, na.action = na.pass)$predictions*100))
      if(any(unique(unlist(soil.fix)) %in% lfix)){
        ## correct probabilities using soil-climate matrix:
        for(k in lfix){
          sel = sapply(soil.fix, function(i){i[which(i==k)]})
          sel <- sel[sapply(sel, function(i){length(i)>0})]
          if(length(sel)>0){ probs[mfix==k,names(sel)] <- 0 }
        }
      }
      ## Standardize ('rescale') values so they sums up to 100:
      rs <- rowSums(probs, na.rm=TRUE)
      m@data <- data.frame(lapply(probs, function(i){i/rs}))
      ## Write GeoTiffs:
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
        ## match most probable class
        m$cl <- col.tbl[match(apply(m@data,1,which.max), col.tbl$int),"Number"]  
        writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group))) 
      }
      gc(); gc()
    }  
  }
}

## This one for testing purposes only:
predict_tile <- function(x, gm1, gm2, gm1.w, gm2.w){
  for(j in 1:length(gm2)){
    gm = gm2[[j]]
    split_predict_c(i=x, gm, in.path="/data/covs1t", out.path="/data/predicted", split_no=j, varn="TAXNWRB")
    rm(gm)
    gc()
  }
  predict_nnet(i=x, gm1, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXNWRB")
  sum_predictions(i=x, in.path="/data/covs1t", out.path="/data/predicted", varn="TAXNWRB", gm1.w=gm1.w, gm2.w=gm2.w, col.legend=col.legend, soil.fix=soil.fix, lev=lev, check.names=TRUE)
}

## simple ranger predict function:
sum_predict_ranger <- function(i, in.path, out.path, varn, num_splits){
  if(length(list.files(path = paste0(out.path, "/", i, "/"), glob2rx(paste0("^", varn, "_*_*_*_*.tif$"))))==0){
    ## load all objects:
    m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
    if(nrow(m@data)>1){
      rf.ls = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", 1:num_splits, ".rds")
      probs <- lapply(rf.ls, readRDS)
      probs <- data.frame(Reduce("+", probs) / num_splits)
      ## Standardize values so they sums up to 100:
      rs <- rowSums(probs, na.rm=TRUE)
      m@data <- data.frame(lapply(probs, function(i){i/rs}))
      ## Write GeoTiffs:
      if(sum(rs,na.rm=TRUE)>0&length(rs)>0){
        tax <- names(m)
        for(j in 1:ncol(m)){
          out <- paste0(out.path, "/", i, "/", varn, "_", tax[j], "_", i, ".tif")
          if(!file.exists(out)){
            m$v <- round(m@data[,j]*100)
            writeGDAL(m["v"], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
          }
        }
      }
      unlink(rf.ls)
    }
  }
}

## special function to derive most probable class:
most_probable_fix <- function(i, in.path, out.path, varn, col.legend, check.names=FALSE, lev){
  out.c <- paste0(out.path, "/", i, "/", varn, "_", i, ".tif")
  if(file.exists(out.c)){
    ## find tiles with missing pixels
    lm <- readGDAL(paste0(in.path, "/", i, "/LMK_", i, ".tif"))
    lm$cl <- readGDAL(out.c)$band1
    if(sum(!is.na(lm$band1)-!is.na(lm$cl))>0){
      ## load all probs:
      m <- raster::stack(list.files(path=paste0(out.path, "/", i), pattern = glob2rx(paste0(varn,  "_*_", i,".tif")), full.names = TRUE))
      names(m) <- sapply(names(m), function(x){strsplit(x, "_")[[1]][2]})
      m <- as(m, "SpatialGridDataFrame")
      ## most probable class:
      if(check.names==TRUE){ 
        Group = gsub("\\.", " ", gsub("\\.$", "\\)", gsub("\\.\\.", " \\(", lev))) 
      } else {
        Group = lev
      }
      col.tbl <- plyr::join(data.frame(Group=Group, int=1:length(lev)), col.legend, type="left")
      ## match most probable class
      m$cl <- col.tbl[match(apply(m@data,1,which.max), col.tbl$int),"Number"]
      unlink(out.c)
      writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(col.tbl$Group)))
      gc(); gc()
    }
  }
}


## -------------------------------
## Numeric variables
## -------------------------------

## 7 standard dephts
split_predict_n <- function(i, gm, in.path, out.path, split_no, varn, sd=c(0, 5, 15, 30, 60, 100, 200), method, multiplier=1, depths=TRUE, DEPTH.col="DEPTH.f", rds.file, SG.col=NULL, SG.col.name){
  if(method=="ranger"){
    if(is.null(split_no)){
      rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_rf.rds")
    } else {
      rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", split_no, ".rds")
    }
  }
  if(method=="xgboost"){
    rds.out = paste0(out.path, "/", i, "/", varn,"_", i, "_xgb.rds")
  }
  if(any(c(!file.exists(rds.out),file.size(rds.out)==0))){
    if(missing(rds.file)){ rds.file = paste0(in.path, "/", i, "/", i, ".rds") }
    if(file.exists(rds.file)&file.size(rds.file)>1e3){ 
      gc(); gc()
      #message("Predicting from .rds file...")
      m <- readRDS(rds.file)
      #if(any(names(m) %in% "OCCGSW7.tif")){ m$OCCGSW7.tif = ifelse(m$OCCGSW7.tif>100, 0, m$OCCGSW7.tif) }
      if(depths==FALSE){
        x <- matrix(data=NA, nrow=nrow(m), ncol=1)
        if(method=="ranger"){
          x[,1] <- round(predict(gm, m@data, na.action=na.pass)$predictions * multiplier)
        } else {
          x[,1] <- round(predict(gm, m@data, na.action=na.pass) * multiplier)
        }
      } else {
        x <- matrix(data=NA, nrow=nrow(m), ncol=length(sd))
        if(missing(SG.col.name)){ SG.col.name=paste0("sg_", varn) }
        for(l in 1:length(sd)){
          m@data[,DEPTH.col] = sd[l]
          if(all(!is.null(SG.col))){
            m@data[,SG.col.name] = m@data[,SG.col[l]]
          }
          if(method=="ranger"){
            v = predict(gm, m@data, na.action=na.pass)$predictions * multiplier
          } else {
            v = predict(gm, m@data, na.action=na.pass) * multiplier
          }
          x[,l] <- round(v)
        }
      }
      saveRDS(x, file=rds.out)
      gc(); gc()
    }
  }
}

## Sum up predictions
sum_predict_ensemble <- function(i, in.path, out.path, varn, num_splits, zmin, zmax, gm1.w, gm2.w, type="Int16", mvFlag=-32768, depths=TRUE, rds.file){
  if(depths==FALSE){
    out.tif = paste0(out.path, "/", i, "/", varn, "_M_", i, ".tif")
    test = !file.exists(out.tif)
  } else {
    test = !length(list.files(path = paste0(out.path, "/", i, "/"), glob2rx(paste0("^",varn,"_M_sl*_*.tif$"))))==7
  }
  if(test){
    if(missing(rds.file)){ rds.file = paste0(in.path, "/", i, "/", i, ".rds") }
    if(file.exists(rds.file)){
      m <- readRDS(rds.file)
      if(nrow(m@data)>1){
        gb = paste0(out.path, "/", i, "/", varn,"_", i, "_xgb.rds")
        if(is.null(num_splits)){
          rf.ls = paste0(out.path, "/", i, "/", varn,"_", i, "_rf.rds")
        } else {
          rf.ls = paste0(out.path, "/", i, "/", varn,"_", i, "_rf", 1:num_splits, ".rds")
        }
        if(all(file.exists(c(rf.ls,gb)))){
          ## import all predictions:
          if(is.null(num_splits)){
            v1 <- readRDS(rf.ls)
          } else {
            v1 <- lapply(rf.ls, readRDS)
            v1 <- Reduce("+", v1) / num_splits
          }
          v2 <- readRDS(gb)
          ## weighted average:
          m@data <- data.frame(Reduce("+", list(v1*gm1.w, v2*gm2.w)) / (gm1.w+gm2.w))
          if(depths==FALSE){
            ## Write GeoTiffs (2D case):
            m@data[,1] <- ifelse(m@data[,1] < zmin, zmin, ifelse(m@data[,1] > zmax, zmax, m@data[,1]))
            writeGDAL(m[1], out.tif, type=type, mvFlag=mvFlag, options="COMPRESS=DEFLATE")
          } else {
            ## Write GeoTiffs (per depth):
            for(l in 1:ncol(m)){
              out.tif = paste0(out.path, "/", i, "/", varn, "_M_sl", l, "_", i, ".tif")
              m@data[,l] <- ifelse(m@data[,l] < zmin, zmin, ifelse(m@data[,l] > zmax, zmax, m@data[,l]))
              writeGDAL(m[l], out.tif, type=type, mvFlag=mvFlag, options="COMPRESS=DEFLATE")
            }
          }
          ## cleanup:
          unlink(rf.ls) 
          unlink(gb)
          gc(); gc()
        }
      }
    }
  }
}
