## Average soil and sedimentary-deposit thickness http://dx.doi.org/10.3334/ORNLDAAC/1304
## Cite as: Pelletier, J.D., P.D. Broxton, P. Hazenberg, X. Zeng, P.A. Troch, G. Niu, Z.C. Williams, M.A. Brunke, and D. Gochis. 2016. Global 1-km Gridded Thickness of Soil, Regolith, and Sedimentary Deposit Layers. ORNL DAAC, Oak Ridge, Tennessee, USA.
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)

gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
system("/usr/local/bin/gdal-config --version")

system("7za e soil_sedimentary_thickness_1km.7z")
GDALinfo("average_soil_and_sedimentary-deposit_thickness.sdat")
