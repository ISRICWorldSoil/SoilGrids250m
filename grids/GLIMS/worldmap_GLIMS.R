## World map of Glaciers http://www.glims.org/download/
## Cite as: Raup, B.H.; A. Racoviteanu; S.J.S. Khalsa; C. Helm; R. Armstrong; Y. Arnaud (2007). "The GLIMS Geospatial Glacier Database: a New Tool for Studying Glacier Change". Global and Planetary Change 56:101--110. (doi:10.1016/j.gloplacha.2006.07.018)
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
library(plotKML)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
system("/usr/local/bin/gdal-config --version")
gdal_rasterize <- "/usr/local/bin/gdal_rasterize"

system("7za e GLIMS.7z")
system("7za e GLIMS_raster.7z")
## import attribute table:
tbl <- read.dbf("glims_polygons.dbf")
str(tbl) ## 307159 units
summary(tbl$dbf$line_type)
tbl$dbf$glims_polygons_ID_500m <- 1:nrow(tbl$dbf)

## Overlay with rasterized polygons (500m)
r <- raster("glims_polygons_ID_500m.sdat")
r
## randomly sample points (we need about 500 points):
glaciers.pnt <- sampleRandom(r, size=200000, sp=TRUE)
plot(glaciers.pnt)
glaciers.pnt@data <- join(glaciers.pnt@data, tbl$dbf, type="left")
str(glaciers.pnt@data)
## Get actual elevation and slope:
dem <- stack(c("/data/MDEM/DEMSRP3a.tif","/data/MDEM/SLPSRM3a.tif"))
ov <- extract(dem, glaciers.pnt)
glaciers.pnt@data <- cbind(glaciers.pnt@data, as.data.frame(ov))
save(glaciers.pnt, file="glaciers.pnt.rda")
## glaciers.pnt[c("line_type","geog_area","min_elev","DEMSRP3a","SLPSRM3a")]
writeOGR(glaciers.pnt, "glaciers.pnt.shp", "glaciers.pnt", "ESRI Shapefile")
#plotKML(glaciers.pnt["geog_area"])
