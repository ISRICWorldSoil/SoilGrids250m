## Processing of the GlobalForestChange global data (http://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.2.html)
## tom.hengl@isric.org

setwd("/mnt/cartman/GlobalForestChange2000-2014")
library(RCurl)
library(rgdal)
library(raster)
library(snowfall)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "gdal_translate"
  gdalwarp =  "gdalwarp"
  gdalbuildvrt = "gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
}

lst = c("loss.txt","lossyear.txt","treecover2000.txt","gain.txt","first.txt","datamask.txt","last.txt") 
## 504 tiles per layer
vars.lst = sapply(lst, function(x){strsplit(x, ".txt")[[1]][1]})

## downloading in parallel can lead to problems, hence double check that local file size local and file size on server are the same
download.wget = function(x){
  if(file.exists(basename(paste(x)))){
    curl_cmd = paste('curl -I HEAD -i', x)
    b = system(paste(curl_cmd, '|grep Content-Length |cut -d : -f 2'), intern = TRUE)
    if(!file.info(basename(paste(x)))$size==as.numeric(b)){
      download.file(method = "wget", x, basename(paste(x))) 
    }
  } else {
    download.file(method = "wget", x, basename(paste(x)))
  }
}

for(j in lst){
  filenames = read.csv(j)
  sapply(paste(filenames[,1]), function(x){ download.wget(x) })
}

## Check that all bands are available:
tif.lst <- list.files(path="./", pattern=glob2rx("*.tif$"), full.names=TRUE)
## 7 x 504 = 3528

sfInit(parallel=TRUE, cpus=10)
sfExport("tif.lst")
sfLibrary(sp)
sfLibrary(rgdal)
gd.lst <- sfLapply(tif.lst, function(i){ try( GDALinfo(i, silent=TRUE)) } )
sfStop()

rm.lst = sapply(gd.lst, function(x){is.null(x)|!attr(x, "projection")=="+proj=longlat +datum=WGS84 +no_defs "|!x[6]==0.00025})
str(tif.lst[which(rm.lst)]) ## 1563 with problems -- end of line problem?
gd.lst[which(rm.lst)[1]]
#x = readGDAL(tif.lst[which(rm.lst)[1]], band = 1)$band1
#system(paste0(gdal_translate, ' ', tif.lst[which(rm.lst)[2]], ' ./cleaned/', gsub(".tif", "_red.tif", tif.lst[which(rm.lst)[2]]), ' -ot Byte -b 1 -co "COMPRESS=DEFLATE"'))
## ERRORS 'GeoTIFF tags apparently corrupt, they are being ignored.' / 'TIFFFetchNormalTag:ASCII value for tag "GeoASCIIParams" contains null byte in value; value incorrectly truncated during reading due to implementation limitations'
tifF.lst = tif.lst[!rm.lst]

## Virtual mosaics:
unlink(list.files(pattern=glob2rx("*.vrt")))
for(j in vars.lst){
  if(!file.exists(paste0(j, '.vrt'))){
    v.lst = tifF.lst[grep(pattern=glob2rx(paste0("*_",j,"_*tif$")), tifF.lst)]
    if(length(v.lst)>0){
      cat(v.lst, sep="\n", file=paste0(j, "_local.txt"))
      system(paste0(gdalbuildvrt, ' -input_file_list ', j, '_local.txt ', j, '.vrt')) 
    }
  }
}
x = lapply(vars.lst, function(j){ print(GDALinfo(paste0(j, ".vrt"))) })
sink(file="GDALinfo_vrts.txt", type="output")
print(x)
sink()

## Per landsat band TAKES >6 hrs!!
vlist = c("loss.vrt", "lossyear.vrt", "treecover2000.vrt", "gain.vrt", "datamask.vrt", rep("first.vrt", 4), rep("last.vrt", 4))
bnames = c("loss", "lossyear", "treecover2000", "gain", "datamask", rep(c("red", "NIR", "SWIR1", "SWIR2"), 2))
bnumbers = c(rep(1, 5), rep(1:4, 2))
years = c(rep("", 5), rep(2000, 4), rep(2014, 4))
sapply(list(vlist, bnames, bnumbers, years), length)

sfInit(parallel=TRUE, cpus=9)
sfExport("bnames", "years", "bnumbers", "vlist", "gdal_translate")
sfLibrary(raster)
sfLibrary(rgdal)
x <- sfClusterApplyLB(1:length(bnames), function(k){ if(!file.exists(paste0('/data/Landsat/100m/Landsat', years[k], '_', bnames[k], '.tif'))){ system(paste0(gdal_translate, ' ', vlist[k], ' /data/Landsat/100m/Landsat', years[k], '_', bnames[k], '.tif -ot Byte -r \"average\" -tr 0.0008333333 0.0008333333 -co \"BIGTIFF=YES\" -b ', bnumbers[k], ' -co "COMPRESS=DEFLATE"')) } } )
sfStop()
#ERROR 1: /mnt/cartman/GlobalForestChange2000-2014/.//Hansen_GFC2015_first_30S_120E.tif, band 3: IReadBlock failed at X offset 0, Y offset 4107

## Make mosaics 100 m USA48:

r <- raster("/projects/USA48/elev48i0100a.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]

usa48.lst = c(rep(NA, 5), paste0(c("REDL00", "NIRL00", "SW1L00", "SW2L00", "REDL14", "NIRL14", "SW1L14", "SW2L14"), ".tif"))
sfInit(parallel=TRUE, cpus=8)
sfExport("bnames", "years", "bnumbers", "usa48.lst", "gdalwarp", "r","cellsize")
sfLibrary(raster)
sfLibrary(rgdal)
x <- sfClusterApplyLB(6:13, function(k){ system(paste0(gdalwarp, ' /data/Landsat/100m/Landsat', years[k], '_', bnames[k], '.tif /projects/USA48/', usa48.lst[k], ' -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } )
sfStop()

## Global mosaics 300 m:

