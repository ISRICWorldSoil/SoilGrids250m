## Spatial predictions and cross-validation errors:

library(rgdal)
library(raster)
library(plotKML)
library(maps)
library(maptools)
library(psych)
library(scales)
library(rgeos)
country <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country$names, ":"), function(x) x[1])
country = as(map2SpatialPolygons(country, IDs=IDs), "SpatialLines")
b_poly <- as(extent(c(-180,180,-65,75)), "SpatialPolygons")
country = gIntersection(country, b_poly, byid = T)
proj4string(country) = "+proj=longlat +datum=WGS84"

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalbuildvrt <- paste0(gdal.dir, "/gdalbuildvrt.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
} else {
  gdalwarp = "/usr/local/bin/gdalwarp"
  gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
  gdal_translate =  "/usr/local/bin/gdal_translate"
}
a.dir <- "/data/shang009/big/soildepth2"# dir of the project
setwd(a.dir)


tvar <- c("BDTICM", "BDRICM", "BDRLOG")
for(i in 2:3)
{
  src_d <- paste0(a.dir, paste0("/", tvar[i], "_M_1km_ll.tif"))
  out_d <- paste0(a.dir, paste0("/", tvar[i], "_M_10km_ll.tif"))
  system(paste("gdalwarp  -tr 0.08333333 0.08333333 -r average -overwrite",
               src_d, out_d))   
}

grd.lst <- c("BDTICM_M_10km_ll.tif", "BDRLOG_M_10km_ll.tif", "BDRICM_M_10km_ll.tif")
pred <- stack( grd.lst)
pred <- as(as(pred, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
#rn = quantile(pred$BDTICM_M_1km_ll_10km, c(.01, .99), na.rm=TRUE)
rn = c(100, 15000)
rx = rev(as.character(round(c(round(rn[1], 0), NA, round(mean(rn), 0), NA, round(rn[2], 0)), 2)))
pred$BDTICM_M_10km_ll <- ifelse(pred$BDTICM_M_10km_ll<rn[1], rn[1], ifelse(pred$BDTICM_M_10km_ll>rn[2], rn[2], pred$BDTICM_M_10km_ll))
pred$BDRICM_M_10km_ll <- ifelse(pred$BDRICM_M_10km_ll>220, 220, pred$BDRICM_M_10km_ll)
jpg(file="Fig_final_predictions.jpg", res=150, width=1200, height=1200*2*(65+75)/180/2)
#png(file="Fig_final_predictions.png", res=150, width=1200, height=1200*3*(65+75)/180/2)
#par(mfrow=c(3,1))
par(mfrow=c(2,1))
par(mai=c(0,0,0,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
image(log1p(raster(pred["BDTICM_M_10km_ll"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="", axes=FALSE, xlab="", ylab="", ylim=c(-65,75)) # , cex.lab=.7, cex.axis=.7
lines(country)
legend("left", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6)
image(raster(pred["BDRLOG_M_10km_ll"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,100), main="", axes=FALSE, xlab="", ylab="", ylim=c(-65,75))
lines(country)
legend("left", rev(c("0","","0.5","","1.0")), fill=rev(SAGA_pal[["SG_COLORS_YELLOW_RED"]][c(1,5,10,15,20)]), horiz=FALSE, cex=.6)
image(raster(pred["BDRICM_M_1km_ll_5kmf"]), col=rev(R_pal[["bpy_colors"]]), zlim=c(0,220), main="", axes=FALSE, xlab="", ylab="", ylim=c(-65,75)) 
lines(country)
legend("left", rev(c("0","","100","","200")), fill=R_pal[["bpy_colors"]][c(1,5,10,15,20)], horiz=FALSE, cex=.8)
dev.off()

## BDTICM alone:
par(mai=c(0,0,0,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
image(log1p(raster(pred["BDTICM_M_1km_ll_5kmf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="", axes=FALSE, xlab="", ylab="", ylim=c(-65,75), asp=1) # , cex.lab=.7, cex.axis=.7
lines(country)
legend("left", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8, bg="white")

## Cross-validation plots:



