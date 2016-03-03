## Function to create RDA files (stacked values of covariates) per tile
## All rows must have all complete cases for the land/soil mask
## List of covariates is at: https://goo.gl/FWD15V
## Tom.Hengl@isric.org

make.covsRDA <- function(in.path, i, csv=TRUE, fix.sums=TRUE){
  out.rda <- paste0(in.path, "/", i, "/", i, ".rds")
  if(!file.exists(out.rda)){
    cov.lst <- list.files(path=paste(in.path, i, sep="/"), glob2rx("*_*_*_*.tif$"), full.names=TRUE)
    mask <- grep("LMK", cov.lst)
    m <- readGDAL(cov.lst[mask], silent=TRUE)
    m <- as(m, "SpatialPixelsDataFrame")
    m$LATWGS84 <- spTransform(m, CRS("+proj=longlat +datum=WGS84"))@coords[,2]
    for(k in 1:length(cov.lst)){
      if(!k==mask){
        covn = strsplit(basename(cov.lst[k]), "_")[[1]][1]
        ## filter out values outside of physical range and round up:
        d <- readGDAL(cov.lst[k], silent=TRUE)$band1[m@grid.index]
        dn <- which(is.na(d)|d<= -10000|d>= 65535)
        if(length(dn)>0){ 
          ## replace all missing values with median value
          ## TH: this could lead to artifacts
          d[dn] <- quantile(d, probs=0.5, na.rm=TRUE)
        }
        m@data[,covn] <- d 
      }
    }
    #gc()
    if(fix.sums==TRUE){
      ## lithological and landform images need to sum-up to 100:
      L <- grep(pattern=glob2rx("L??USG5"), names(m))
      sL <- rowSums(m@data[,L])
      rL <- range(sL, na.rm=TRUE)
      if(rL[1]<100|rL[2]>100){
        for(j in L){ m@data[sL>0,j] <- round(m@data[sL>0,j] / sL[sL>0] * 100, 0) }
      }
    }
    ## Save output:
    saveRDS(m, file=out.rda, compress=TRUE)
    if(csv==TRUE){
      outfile <- paste0(in.path, "/", i, "/", i, ".csv")
      df <- as.data.frame(m)
      ## round up numbers to save space:
      df <- data.frame(lapply(df, function(x){signif(x, digits=4)}))
      write.csv(df, outfile)
      gc()
      gzip(outfile)
    }
  }
}

## gz files for h2o (ca 2 hrs to make):
make.csv.gz <- function(i, in.path){
  infile <- paste0(in.path, "/", i, "/", i, ".rds")
  outfile <- paste0(in.path, "/", i, "/", i, ".csv")
  if(!file.exists(paste0(outfile, ".gz"))){
    m <- readRDS(infile)
    df <- as.data.frame(m)
    ## round up numbers to save space:
    df <- data.frame(lapply(df, function(x){signif(x, digits=4)}))
    write.csv(df, outfile)
    gc()
    gzip(outfile)
  }
}