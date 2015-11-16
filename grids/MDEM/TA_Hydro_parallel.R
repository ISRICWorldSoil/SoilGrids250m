## Generation of DEM parameters:
library(RSAGA)
## http://www.saga-gis.org/saga_module_doc/2.2.0/
library(R.utils)
library(snowfall)
#library(gdalUtils)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
myenv =  rsaga.env(path="/usr/bin", parallel=TRUE, cmd = "/usr/local/bin/saga_cmd") ## modules="/usr/local/lib/saga", 
load("/data/models/equi7t3.rda")
## regions:
regs <- c("AF","AN","AS","EU","NA","OC","SA")

## Positive and negative openess:
saga_Openess <- function(inputFile){
  filen <- strsplit(inputFile, ".gz")[[1]][1]
  if(!file.exists(set.file.extension(filen, ".sdat"))){
    gunzip(inputFile, overwrite=TRUE, remove=FALSE)
    system(paste(gdal_translate, filen, set.file.extension(filen, ".sdat"), '-ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
  }
  if(!file.exists(set.file.extension(gsub("DEM", "POS", filen), ".zip"))){
    system(paste0("/usr/local/bin/saga_cmd -c=7 -f=q ta_lighting 5 -DEM ", set.file.extension(filen, ".sgrd"), " -POS ", set.file.extension(gsub("DEM", "POS", filen), ".sgrd"), " -NEG ", set.file.extension(gsub("DEM", "NEG", filen), ".sgrd"), " -RADIUS 20000 -METHOD 0"))
    POS.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "POS", filen))
    try( zip(set.file.extension(gsub("DEM", "POS", filen), ".zip"), files=POS.lst, zip="/usr/bin/zip") )
    NEG.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "NEG", filen))
    try( zip(set.file.extension(gsub("DEM", "NEG", filen), ".zip"), files=NEG.lst, zip="/usr/bin/zip") )
  }
  unlink(filen)
}

## run in parallel:
sfInit(parallel=TRUE, cpus=38)
sfLibrary(RSAGA)
sfLibrary(R.utils)
sfLibrary(utils)
sfExport("regs", "saga_Openess", "gdal_translate")
t <- sfLapply(regs, function(i){saga_Openess(inputFile=paste0("MDEM_",i,"_250m.tif.gz"))})
sfStop()

## OR, run in loop:
for(i in regs){
  saga_Openess(inputFile=paste0("MDEM_",i,"_250m.tif.gz"))
}

## Convert to LatLon system:
tile2.ll <- function(t, tvar){
  s_srs = proj4string(equi7t3[[t]])
  nfile <- paste0(tvar, "_", names(equi7t3)[t], "_250m.sdat")
  ofile <- gsub("250m.sdat", "ll.sdat", nfile)
  if(!file.exists(ofile)){ 
    system(paste0(gdalwarp, ' ', nfile, ' ', ofile, ' -of \"SAGA\" -s_srs \"', s_srs, '\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.008333333 0.008333333 -te -180 -90 180 90'))
  }
}

sfInit(parallel=TRUE, cpus=6)
sfLibrary(gdalUtils)
sfLibrary(sp)
sfExport("equi7t3", "tile2.ll", "gdalwarp")
#x <- sfLapply(1:length(equi7t3), tile2.ll, tvar="MTWI")
x <- sfLapply(1:length(equi7t3), tile2.ll, tvar="MSLP")
sfStop()
#system(paste0("7za a -t7z MTWI_EU_ll.7z MTWI_EU_ll.tif -mmt=40"))

## Resample to 1km:
#tmp.lst <- list.files(pattern=glob2rx("MTWI_*_ll.sgrd$"), full.names=TRUE)
tmp.lst <- list.files(pattern=glob2rx("MSLP_*_ll.sgrd$"), full.names=TRUE)
## TH: This operation takes a bit more time but produces best merge:
system(paste0('/usr/local/bin/saga_cmd -c=40 grid_tools 3 -GRIDS \"', paste(tmp.lst, collapse=";"), '\" -TYPE=7 -INTERPOL=0 -OVERLAP=4 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=-179.9958333335 -TARGET_USER_XMAX=179.9958189335 -TARGET_USER_YMIN=-89.9958261335 -TARGET_USER_YMAX=89.9958333335 -TARGET_OUT_GRID=\"SLPSRM3a.sgrd\" -TARGET_USER_SIZE=0.008333333'))
unlink("MSLP_*_ll.*")
#system(paste0("7za a -t7z TWISRM3a.7z TWISRM3a.* -mmt=40"))
system(paste0(gdalwarp, ' TWISRM3a.sdat TWISRM3a.tif -r \"near\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' SLPSRM3a.sdat SLPSRM3a.tif -r \"near\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))




