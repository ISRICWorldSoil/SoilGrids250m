## ESA land cover data (https://www.esa-landcover-cci.org/?q=node/175)
## tom.hengl@isric.org

setwd("/data/ESA/LandCover")
load(".RData")
library(rgdal)
library(raster)
#library(animation)
library(plotKML)
library(GSIF)
library(XML)
library(plyr)
library(maps)
library(maptools)

obj <- GDALinfo("ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif")
tile.lst <- getSpatialTiles(obj, block.x=5, return.SpatialPolygons=TRUE)
## Generating 2592 tiles...
tile.tbl <- getSpatialTiles(obj, block.x=5, return.SpatialPolygons=FALSE)
tile.tbl$xc = (tile.tbl$xu+tile.tbl$xl)/2
tile.tbl$yc = (tile.tbl$yu+tile.tbl$yl)/2
tile.tbl$ID = as.character(1:nrow(tile.tbl))
## plot tiles:
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
png(file = "World_5degree_tiles.png", width = 3600, height = 1900, type="cairo")
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="red")
lines(tile.lst, col="black")
text(x=tile.tbl$xc, y=tile.tbl$yc, labels=paste0("T", tile.tbl$ID))
dev.off()
## export tiles:
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
unlink("tiles_ll_500km.shp")
writeOGR(tile.pol, "tiles_ll_500km.shp", "tiles_ll_500km", "ESRI Shapefile")

## Finer resolution country borders:
country.f = readOGR("/data/GAUL/g2015_2014_0/g2015_2014_0.shp", "g2015_2014_0")

## legend:
leg = read.csv("ESACCI-LC-Legend.csv", header=TRUE, sep=";", row.names=NULL)
leg$COLOR <- rgb(red=leg$R/255, green=leg$G/255, blue=leg$B/255)

## plot tile using the original legend:
plot_tile <- function(m, year, leg, width, height, asp){
  out.file = tempfile(fileext = ".png")
  m$plot = as.factor(m@data[,1])
  m.leg = plyr::join(data.frame(NB_LAB=levels(m$plot)), leg)$COLOR
  png(file = out.file, width = width, height = height, type="cairo")
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  image(raster(m["plot"]), asp=asp, col=m.leg)
  plotrix::boxed.labels(x=m@grid@cellcentre.offset[1]+m@grid@cellsize[1]*m@grid@cells.dim[1]/2, y=m@grid@cellcentre.offset[2]+m@grid@cellsize[2]*(m@grid@cells.dim[2]-30), labels=year, cex=3)
  dev.off()
  return(out.file)
}

## get the correct aspect ratio (https://stackoverflow.com/questions/31745894/get-aspect-ratio-for-lat-long-plots)
library(ggplot2)
map_aspect = function(x, y) {
  x.center <- sum(range(x)) / 2
  y.center <- sum(range(y)) / 2
  x.dist <- ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist <- ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  y.dist / x.dist
}

## plot a tile as animation
tile.gif = function(fname="ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif", i, tile.tbl, leg, out.file, years=1992:2015){
  if(missing(out.file)){ out.file = paste0("./tiled/T", tile.tbl$ID[i], "_landcover_1992_2015-v2.0.7.gif") }
  if(!file.exists(out.file)){
    m = readGDAL(fname=fname, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    sd = sapply(1:ncol(m), function(x){sd(m@data[,x], na.rm=TRUE)})
    if(!all(sd==0)){
      asp = map_aspect(tile.tbl$xc[i], tile.tbl$yc[i])
      out.lst = sapply(1:ncol(m), function(j){plot_tile(m[j], year=years[j], leg, width=1800, height=1800*asp, asp)})
      ## Create animation:
      system(paste0('convert -delay 200 ', paste(out.lst, collapse=" "), ' ', out.file))
      unlink(out.lst)
      gc()
    }
 }
}

library(snowfall)
sfInit(parallel=TRUE, cpus=24)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(raster)
sfLibrary(ggplot2)
sfLibrary(plotrix)
sfExport("tile.gif", "plot_tile", "tile.tbl", "leg", "map_aspect")
out <- sfClusterApplyLB(1:nrow(tile.tbl), function(i){try( tile.gif(i=i, tile.tbl=tile.tbl, leg=leg) )})
sfStop()

## Tiling system MODIS projection ----
objS <- GDALinfo("ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7_sin.tif")
tileS.lst <- getSpatialTiles(objS, block.x=200000, return.SpatialPolygons=TRUE)
tileS.tbl <- getSpatialTiles(objS, block.x=200000, return.SpatialPolygons=FALSE)
tileS.tbl$ID = as.character(1:nrow(tileS.tbl))
## 20,301 tiles
tileS.pol = SpatialPolygonsDataFrame(tileS.lst, tileS.tbl)
writeOGR(tileS.pol, "tiles_sin.shp", "tiles_sin", "ESRI Shapefile")
## remove tiles without values
system(paste('gdal_translate ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7_sin.tif LCC_300m_sin.sdat -of \"SAGA\" -ot \"Byte\"'))
#GDALinfo("LCC_300m_sin.sdat")
system(paste0('saga_cmd -c=48 shapes_grid 2 -GRIDS=\"LCC_300m_sin.sgrd\" -POLYGONS=\"tiles_sin.shp\" -PARALLELIZED=1 -RESULT=\"ov_ADMIN_tiles_sin.shp\"'))
ov_ADMIN_sin = readOGR("ov_ADMIN_tiles_sin.shp", "ov_ADMIN_tiles_sin")
summary(selS.t <- (!ov_ADMIN_sin$LCC_300m_si.5==0)&(!ov_ADMIN_sin$LCC_300m_si.5==210))
## 6400 tiles with values
ov_ADMIN_sin = ov_ADMIN_sin[selS.t,]
tS.sel = as.character(ov_ADMIN_sin$ID)
newS.dirs <- paste0("/data/tt/LDN/Stiled/T", tS.sel)
x <- lapply(newS.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
saveRDS(ov_ADMIN_sin, "Sinusoidal_tiles_200km.rds")
#writeOGR(ov_ADMIN_sin, "Sinusoidal_tiles_200km.shp", "Sinusoidal_tiles_200km", "ESRI Shapefile")
unlink("Sinusoidal_tiles_200km.gpkg")
writeOGR(ov_ADMIN_sin, "Sinusoidal_tiles_200km.gpkg", "Sinusoidal_tiles_200km", "GPKG")
