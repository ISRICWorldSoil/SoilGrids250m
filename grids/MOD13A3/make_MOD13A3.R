## Preprocessing of the global NDVI time series at 1 km resolution
## tom.hengl@isric.org

setw("/data/MOD13A3/")
library(rgdal)
library(RSAGA)
library(snowfall)
load(".RData")

## Sync all hdf files (AFTER DOWNLOADING WITH FILEZILLA)
## 859GB IN SIZE
y.span = 2000:2016
for(i in y.span){
  dr.lst <- normalizePath(list.dirs(path=paste0("/mnt/cartman/MODIS/6/MOD13A3/", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "/")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://anonymous@ladsweb.nascom.nasa.gov/allData/6/MOD13A3/', i ,'/', x))
  }
}

## list files
hdf.lst <- list.files(path="/mnt/cartman/MODIS/6/MOD13A3", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
length(hdf.lst)
## 56649 files
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=9, byrow=T))
hdfy.lst <- hdfy.lst[,7:9]
names(hdfy.lst) <- c("Year", "Day", "FileName")
hdfy.lst$YearDay <- as.factor(paste(hdfy.lst$Year, hdfy.lst$Day, sep="_"))
lvs <- levels(as.factor(paste(hdfy.lst$YearDay)))
## 201 periods
str(hdfy.lst)
summary(hdfy.lst$YearDay, maxsum=length(levels(hdfy.lst$YearDay)))
save.image()

resample_PRM = function(INPUT_FILENAME, OUTPUT_FILENAME){
  prm = set.file.extension(OUTPUT_FILENAME, ".prm")
  filename = file(prm, open="wt")
  write(paste0('INPUT_FILENAME = \"', INPUT_FILENAME, '\"'), filename)
  write(paste0('SPECTRAL_SUBSET = ( 1 )'), filename)
  write(paste0('SPATIAL_SUBSET_TYPE = INPUT_LAT_LONG'), filename)
  write(paste0('SPATIAL_SUBSET_UL_CORNER = ( 79.999999993 -179.9 )'), filename)
  write(paste0('SPATIAL_SUBSET_LR_CORNER = ( -59.999999995 179.9 )'), filename)
  write(paste0('OUTPUT_FILENAME = \"', OUTPUT_FILENAME, '\"'), filename)
  write(paste0('RESAMPLING_TYPE = BILINEAR'), filename)
  write(paste0('OUTPUT_PROJECTION_TYPE = GEO'), filename)
  write(paste0('OUTPUT_PROJECTION_PARAMETERS = ( 
    0.0 0.0 0.0
    0.0 0.0 0.0
    0.0 0.0 0.0
    0.0 0.0 0.0
    0.0 0.0 0.0 )'), filename)
  write(paste0('DATUM = NoDatum'), filename)
  close(filename)
}

## run per period (using MODIS MRT):
MODIS_mosaic = function(i){
  out.file = paste0("/data/MOD13A3/NDVI_", i, "_1km.hdf")
  out.tif =  paste0('/data/MOD13A3/MOD13A3_', i, '_ll.tif')
  if(!file.exists(out.file)){
    tmp.lst <- hdf.lst[which(hdfy.lst$YearDay==i)]
    out.lst <- paste0("/data/MOD13A3/NDVI_", i, "_list.prm")
    cat(tmp.lst, sep="\n", file=out.lst)
    try( system(paste0('/home/tom/SW/MRT/bin/mrtmosaic -i ', out.lst, ' -s \"1 0 0 0 0 0 0 0 0 0 0 0\" -o ', out.file)) )
    if(!file.exists(out.tif)){
      resample_PRM(out.file, out.tif)
      try( system(paste0('/home/tom/SW/MRT/bin/resample -p ', set.file.extension(out.tif, ".prm"))) )
    }
  }
}

# sfInit(parallel=TRUE, cpus=6)
# sfExport("MODIS_mosaic", "resample_PRM", "hdfy.lst", "hdf.lst")
# sfLibrary(RSAGA)
# sfLibrary(rgdal)
# t <- sfLapply(levels(hdfy.lst$YearDay), MODIS_mosaic)
# sfStop()
t = lapply(levels(hdfy.lst$YearDay), MODIS_mosaic)

## Compress size:
NDVI.lst = list.files("/data/MOD13A3/", pattern = glob2rx("MOD13A3_*_*_ll.1_km_monthly_NDVI.tif"), full.names = TRUE)
## 344
library(snowfall)
sfInit(parallel=TRUE, cpus=8)
sfExport("NDVI.lst")
t <- sfLapply(NDVI.lst, function(x){ if(!file.exists(gsub(".1_km", "_1km", x))) { system(paste0('gdal_translate ', x, ' ', gsub(".1_km", "_1km", x), ' -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"')) } })
sfStop()
unlink(NDVI.lst)
