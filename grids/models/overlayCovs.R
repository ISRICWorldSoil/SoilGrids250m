## Prepare points for overlays:
## tom.hengl@isric.org

## check all RDA files (complete consistent)
check.RDA <- function(i){
  r <- raster(i)
  na.count <- sum(getValues(is.na(r)|r==0))
  if(na.count==ncell(r)){
    return(0)
  } else {
    if(na.count==0){
      return(100)
    } else {
      return(signif((1-na.count/ncell(r))*100, 3))
    } 
  }
}
sfInit(parallel=TRUE, cpus=40)
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("check.LMK", "msk.lst")
selL <- sfLapply(msk.lst, check.LMK)
sfStop()
## 243 empty tiles on 1th November 2015
selD <- data.frame(name=msk.lst[which(selL==0)])

## list all available tifs:
cov.lst <- as.vector(unlist(lapply(y, function(i){list.files(path="/data/covs", pattern=glob2rx(paste0(i, "_*_*_*.tif$")), recursive=TRUE, full.names=TRUE)})))
save(cov.lst, file="cov.lst.rda")

## remove all directories with empty landmask (CAREFULL!)
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
pr.dirs.c <- list.dirs("/data/predicted")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})

