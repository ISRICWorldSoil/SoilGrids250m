## Rasterize GAUL (http://www.fao.org/geonetwork/srv/en/metadata.show?currTab=simple&id=12691) and mask out water bodies
## The GAUL dataset is distributed to the United Nations and other authorized international and national institutions/agencies. FAO grants a license to use, download and print the materials contained in the GAUL dataset solely for non-commercial purposes and in accordance with the conditions specified in the data license. The full GAUL Data License document is available for downloading.
## tom.hengl@isric.org

setwd("/data/countries")
library(raster)
library(rgdal)
library(foreign)
library(plotKML)
load(".RData")

## GAUL countries:
r = raster("/mnt/cartman/stacked250m/LCEE10.tif")
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
## Country names:
system("ogrinfo /mnt/cartman/FAO/GAUL/g2015_2014_1/g2015_2014_1.shp")
gaul.tbl = read.dbf("/mnt/cartman/FAO/GAUL/g2015_2014_1/g2015_2014_1.dbf")
## 3422 rows
str(gaul.tbl)
str(levels(gaul.tbl$ADM1_NAME))
## 3218
## There are many duplicated names unfortunately
summary(duplicated(gaul.tbl$ADM1_NAME))
gaul.tbl$ADM0_ADM1 = as.factor(make.unique(paste(gaul.tbl$ADM0_NAME, gaul.tbl$ADM1_NAME, sep="_")))
summary(duplicated(gaul.tbl$ADM0_ADM1))
## 276 'countries'
#summary(as.factor(gaul.tbl$ADM1_CODE))
gaul.tbl$Value = as.integer(gaul.tbl$ADM0_ADM1)
summary(gaul.tbl$Value)
write.dbf(gaul.tbl, "/mnt/cartman/FAO/GAUL/g2015_2014_1/g2015_2014_1.dbf")
hb.leg = gaul.tbl[,c("Value","ADM0_NAME","ADM1_NAME")]
## https://gist.github.com/primaryobjects/b56189e52a3f0e3cfdbf
friendlyUrl <- function(text, sep = '-', max = 80) {
  # Replace non-alphanumeric characters.
  url <- gsub('[^A-Za-z0-9]', sep, text)
  # Remove double separators (do this twice, in case of 4 or 3 repeats).
  doubleSep <- paste(sep, sep, sep = '')
  url <- gsub(doubleSep, sep, url)
  url <- gsub(doubleSep, sep, url)
  # Convert to lowercase and trim to max length.
  url <- substr(tolower(url), 1, max)
  # Trim leading and trailing separators.
  gsub('^-+|-$', '', url)
}
## URL friendly names:
hb.leg$ADM0_NAME_filename = friendlyUrl(iconv(hb.leg$ADM0_NAME, to="ASCII", sub=""))
hb.leg$ADM1_NAME_filename = friendlyUrl(paste(iconv(hb.leg$ADM0_NAME, to="ASCII", sub=""), iconv(hb.leg$ADM1_NAME, to="ASCII", sub=""), sep="__"))
hb.leg$ADM0_NAME_Value = as.integer(hb.leg$ADM0_NAME)
View(hb.leg)
write.csv(hb.leg, "g2015_2014_1_legend.csv")

## takes 1hr!
system(paste0('saga_cmd -c=8 grid_gridding 0 -INPUT \"/mnt/cartman/FAO/GAUL/g2015_2014_1/g2015_2014_1.shp\" -FIELD \"Value\" -GRID \"g2015_2014_1_250m.sgrd\" -GRID_TYPE 2 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("g2015_2014_1_250m.tif")
system(paste0('gdal_translate g2015_2014_1_250m.sdat g2015_2014_1_250m.tif -ot \"Int16\" -co \"COMPRESS=DEFLATE\"'))
unlink("g2015_2014_1_250m.sdat")
raster("g2015_2014_1_250m.tif")

## Filter out water bodies ----
obj <- GDALinfo("g2015_2014_1_250m.tif")
tiles <- GSIF::getSpatialTiles(obj, block.x=4, return.SpatialPolygons = FALSE)
tiles.pol <- GSIF::getSpatialTiles(obj, block.x=4, return.SpatialPolygons = TRUE)
tile.pol = SpatialPolygonsDataFrame(tiles.pol, tiles)
## Function to paralelize:
fun_mask <- function(i, tiles, dir="./tiled/", lc="/mnt/cartman/stacked250m/LCEE10.tif", ga="g2015_2014_1_250m.tif", hb.leg){
  out.tif = paste0(dir, "T", i, ".tif")
  if(!file.exists(out.tif)){
    x = readGDAL(lc, offset=unlist(tiles[i,c("offset.y","offset.x")]), region.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    x$gaul = readGDAL(ga, offset=unlist(tiles[i,c("offset.y","offset.x")]), region.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tiles[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    x$mask = ifelse(x$band1==0|x$band1==210|is.na(x$band1), NA, x$gaul)
    x$mask2 = plyr::join(data.frame(Value=x$mask), hb.leg[,c("Value","ADM0_NAME_Value")], match="first")$ADM0_NAME_Value
    if(!all(is.na(x$mask))){
      writeGDAL(x["mask"], type="Int16", mvFlag=-32768, out.tif, options=c("COMPRESS=DEFLATE"))
      writeGDAL(x["mask2"], type="Int16", mvFlag=-32768, gsub("T", "M", out.tif), options=c("COMPRESS=DEFLATE"))
    }
  }
}

#dir.create("./tiled")
#unlink(list.files("./tiled", pattern=".tif", full.names=TRUE))
fun_mask(i=801, tiles=tiles, hb.leg=hb.leg)
library(parallel)
x0 = mclapply(1:nrow(tiles), FUN=fun_mask, tiles=tiles, hb.leg=hb.leg)
## Mosaick back results of computing:
t.lst <- list.files(path="./tiled", pattern=glob2rx("^T*.tif$"), full.names=TRUE, recursive=TRUE)
## 1651
cat(t.lst, sep="\n", file="mask_tiles.txt")
system('gdalbuildvrt -input_file_list mask_tiles.txt mask.vrt')
system(paste0('gdalwarp mask.vrt GAUL_ADMIN1_landmask_250m.tif -ot \"Int16\" -dstnodata -32768 -co \"BIGTIFF=YES\" -r \"near\" -overwrite -co \"COMPRESS=DEFLATE\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner))
## Countries:
m.lst <- list.files(path="./tiled", pattern=glob2rx("^M*.tif$"), full.names=TRUE, recursive=TRUE)
cat(m.lst, sep="\n", file="mask2_tiles.txt")
system('gdalbuildvrt -input_file_list mask2_tiles.txt mask2.vrt')
system(paste0('gdalwarp mask2.vrt GAUL_ADMIN0_landmask_250m.tif -ot \"Int16\" -dstnodata -32768 -co \"BIGTIFF=YES\" -r \"near\" -overwrite -co \"COMPRESS=DEFLATE\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner))
