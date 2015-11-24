## Create output dirs based on the land mask
## 2594 folders / tiles

library(RSAGA)
library(R.utils)
gdal_translate =  "/usr/local/bin/gdal_translate"

zipfile <- "/data/models/LMK_tiled_250m.zip"
zp.lst <- unzip(zipfile, list=TRUE)
dir.lst <- sapply(zp.lst$Name, function(x){strsplit(paste(strsplit(x, "_")[[1]][2:4], collapse="_"), ".tif")[[1]][1]})
#sapply(1:nrow(zp.lst), function(x){unzip(zipfile, files=zp.lst$Name[x], overwrite=TRUE, exdir=paste0("./", dir.lst[x]))})
x <- sapply(1:nrow(zp.lst), function(x){unzip(zipfile, files=zp.lst$Name[x], overwrite=FALSE, exdir=paste0("/data/covs/", dir.lst[x]))})
## prepare mask files in SAGA GIS:
msk.lst <- list.files(path="/data/covs/", pattern=glob2rx("LMK_*_*_*.tif$"), full.names=TRUE, recursive=TRUE)
## Prepare masks as SAGA GIS grids
x <- sapply(msk.lst, function(x){if(!file.exists(set.file.extension(x, ".sdat"))){ system(paste0(gdal_translate, ' ', x, ' ', set.file.extension(x, ".sdat"), ' -of \"SAGA\" -ot \"Byte\"'))}})
## clean-up:
#rm.lst <- list.files(path="/data/covs", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(rm.lst[grep(rm.lst, pattern="LMK_")])