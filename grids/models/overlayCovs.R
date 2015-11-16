## Prepare points for overlays:

## list all available tifs:
cov.lst <- as.vector(unlist(lapply(y, function(i){list.files(path="/data/covs", pattern=glob2rx(paste0(i, "_*_*_*.tif$")), recursive=TRUE, full.names=TRUE)})))
save(cov.lst, file="cov.lst.rda")

## remove all directories with empty landmask (CAREFULL!)
pr.dirs <- basename(dirname(list.files(path="/data/covs", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
pr.dirs.c <- list.dirs("/data/predicted")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})

