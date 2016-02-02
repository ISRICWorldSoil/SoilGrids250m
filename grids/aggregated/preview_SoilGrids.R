## Preview SoilGrids predictions (PNGs)
## Tom.Hengl@isric.org

library(rgdal)
library(GSIF)
library(utils)
library(R.utils)
library(snowfall)
library(raster)
library(RSAGA)
library(grDevices)
library(plotKML)
library(maps)
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdal_translate = "/usr/local/bin/gdal_translate"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
soil.legends = soil.legends

## List of property maps:
tif.lst <- list.files(path="/data/GEOG", pattern="1km_ll", full.names=TRUE, recursive=TRUE)
## 224
#tif.lst[grep("WRB",tif.lst)]

aggr_SG <- function(i, r, tr=0.05, ti="1km", tn="5km"){
  out = gsub(paste0("_", ti), paste0("_", tn), basename(i))
  if(missing(r)){
    if(any(basename(i) %in% c("TAXOUSDA_M_250m_ll.tif", "TAXNWRB_250m_ll.tif"))){
      r = 'near'
    } else {
      r = 'average'
    }
  }
  if(!file.exists(out)){ 
    system(paste0(gdalwarp, ' ', i, ' ', set.file.extension(gsub(paste0("_", ti), paste0("_", tn), basename(i)), ".tif"), ' -r \"', r, '\" -tr ', tr, ' ', tr, ' -co \"COMPRESS=DEFLATE\"')) 
  }
}

sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "gdalwarp", "aggr_SG")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(tif.lst, aggr_SG)
sfStop()


## create PNGs (5km res)
out.lst <- list.files(pattern="5km_ll", full.names=TRUE, recursive=TRUE)

plot_SG <- function(i, res=150, zlim=c(0,40), width=7200, height=2987, ylim=c(-60, 85), replace=TRUE){
  out.file = paste0(normalizeFilename(basename(i)), ".png")
  if(replace==TRUE){
    r <- raster(i)
    isF = length(grep(pattern="WRB", basename(i)))>0 | length(grep(pattern="USDA", basename(i)))>0
    if(isF){
      col_scale = R_pal[["bpy_colors"]]
      breaks = c(seq(zlim[1], zlim[2], length=20), 100)
    } else {
      if(length(grep(pattern="ORC", basename(i)))>0){
        col_scale = SAGA_pal[[1]]
        breaks = c(0, 1, 2, 4, 6, 8, 11, 14, 18, 21, 25, 32, 40, 54, 70, 100, 140, 180, 240, 300, 600)
      }
      for(j in c("PHIHOX","PHIKCL","BLD","CEC","SNDPPT","SLTPPT","CLYPPT")){
        if(length(grep(pattern=j, basename(i)))>0){
          col_scale = soil.legends[[j]]$COLOR
          breaks = c(soil.legends[[j]]$MIN[1], soil.legends[[j]]$MAX)
        }
      }
      if(length(grep(pattern="CRFVOL", basename(i)))>0){
        col_scale = SAGA_pal[[2]]
        breaks = c(0, 1, 2, 4, 6, 8, 10, 14, 19, 25, 32, 40, 54, 70, 88, 120, 190, 270, 430, 600, 1000)/10
      }
      if(length(grep(pattern="BDT", basename(i)))>0){
        col_scale = SAGA_pal[[1]]
        breaks = c(0, 1, 2, 4, 6, 8, 14, 20, 34, 46, 64, 85, 105, 140, 180, 220, 310, 450, 650, 950, 1800)*100
      }
      if(length(grep(pattern="BDR", basename(i)))>0){
        col_scale = rev(R_pal[["heat_colors"]])
        breaks = seq(0,260,length=21)
      }
      if(length(grep(pattern="BDRLOG", basename(i)))>0){
        col_scale = R_pal[["heat_colors"]]
        breaks = seq(0,100,length=21)
      }
    }
    ri = cut(r, breaks=breaks, include.lowest=TRUE, right=FALSE)
    png(file = out.file, res = res, width = width, height = height, type="cairo")
    #dev.new(width = 20, height = 8.3)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    image(ri, col=col_scale, ylim=ylim)
    lines(country, col="black")
    legend("bottomleft", legend=signif(rev(breaks[-1]), 3), fill=rev(col_scale), horiz=FALSE, pt.cex=2)
    dev.off()
  }
}

sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "gdalwarp", "plot_SG", "country", "soil.legends")
sfLibrary(rgdal)
sfLibrary(raster)
sfLibrary(plotKML)
sfLibrary(grDevices)
out <- sfClusterApplyLB(out.lst, try( plot_SG ) )
sfStop()
