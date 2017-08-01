## Preprocessing MERIT DEM (http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)
## tom.hengl@isric.org

library(rgdal)
library(parallel)
setwd("/mnt/cartman/MERIT")

## Download all tiles into the dir:
system('wget -N -nv -r -np -nH --accept "dem_tif_*.tar.gz" --reject="index.html" --cut-dirs=4 --level=5 --http-user=*** --http-password=*** "http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/"')
gz.lst = list.files(pattern=glob2rx("dem_tif_*.tar.gz"))
## Untar (in parallel not a good idea unfortunately):
#mclapply(gz.lst, FUN=function(i){untar(i)}, mc.cores=6)
lapply(gz.lst, function(i){untar(i)})

## Make a mosaic:
dem.lst <- list.files(path="/mnt/cartman/MERIT", pattern=glob2rx("*_dem.tif$"), full.names=TRUE, recursive=TRUE)
## 11150 tiles
## check if the untar went ok:
dem.size.lst = sapply(dem.lst, file.size)
str(dem.lst[dem.size.lst<1.4e8])
unlink(dem.lst[dem.size.lst<1.4e8])
#mclapply(gz.lst, FUN=function(i){untar(i, extras="--skip-old-files")}, mc.cores=6)
#lapply(gz.lst, function(i){untar(i, extras="--skip-old-files")})

cat(dem.lst, sep="\n", file="MERIT_tiles.txt")
system('gdalbuildvrt -input_file_list MERIT_tiles.txt MERIT_100m.vrt')
#system('gdalinfo MERIT_100m.vrt')
## Size is 432000, 174000
## Pixel Size = (0.000833333353512,-0.000833333353512)
#system(paste0('gdalwarp MERIT_100m.vrt MERIT_dem_1km_v28_July_2017.tif -ot \"Int16\" -co \"BIGTIFF=YES\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\" -tr ', 1/120, ' ', 1/120))
system('gdalwarp MERIT_100m.vrt  -ot \"Int16\" -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\"')
## takes 6-7 hours... 168GB file?!
system('gdal_translate MERIT_dem_100m_v28_July_2017.tif /data/Landsat/100m/MERIT_dem_100m_v28_July_2017.tif -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -co \"NUM_THREADS=6\"')
system(paste0('gdaladdo /data/Landsat/100m/MERIT_dem_100m_v28_July_2017.tif 2 4 8 16 32 64 128'))
## File now 19GB
## Add metadata:
md.Fields = c("SERIES_NAME", "ATTRIBUTE_UNITS_OF_MEASURE", "CITATION_URL", "CITATION_ORIGINATOR",	"CITATION_ADDRESS",	"PUBLICATION_DATE", "PROJECT_URL", "DATA_LICENSE")
md.Values = c("MERIT DEM: Multi-Error-Removed Improved-Terrain DEM", "meter", "http://dx.doi.org/10.1002/2017GL072874", "Institute of Industrial Sciences, The University of Tokyo", "yamadai [at] rainbow.iis.u-tokyo.ac.jp", "15 May, 2017", "http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/", "https://creativecommons.org/licenses/by-nc/4.0/")
m = paste('-mo ', '\"', md.Fields, "=", md.Values, '\"', sep="", collapse = " ")
command = paste0('gdal_edit.py ', m,' /data/Landsat/100m/MERIT_dem_100m_v28_July_2017.tif')
system (command, intern=TRUE)
system('gdalinfo /data/Landsat/100m/MERIT_dem_100m_v28_July_2017.tif')
