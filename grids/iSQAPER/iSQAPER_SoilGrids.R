## Soil data (SoilGrids) for the iSQAPER project
## Tom.Hengl@isric.org

library(rgdal)
library(GSIF)
library(utils)
library(R.utils)
library(snowfall)
library(raster)
library(RSAGA)
library(plotKML)
library(R.utils)
plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "/usr/local/bin/gdalwarp"
gdal_translate = "/usr/local/bin/gdal_translate"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
CN.prj = "+proj=lcc +lat_1=25 +lat_2=47 +lat_0=10 +lon_0=110 +x_0=500000 +y_0=500000 +ellps=WGS84 +units=m +no_defs"
EU.prj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

## china bounding box:
#cn <- readOGR("china_lat_longtitde.shp", "china_lat_longtitde")
cn <- readOGR("China_boundary.shp", "China_boundary")
cn.lcc <- spTransform(cn, CRS(CN.prj))
unlink("China_boundary_LCC.shp")
writeOGR(cn.lcc, "China_boundary_LCC.shp", "China_boundary_LCC", "ESRI Shapefile") 
## eu bounding box:
eu.laea <- readOGR("European_boundary.shp", "European_boundary")
proj4string(eu.laea)
#plot(eu.laea)
eu.ll <- readOGR("European_boundary_ll.shp", "European_boundary_ll")
#eu <- spTransform(eu.ll, CRS(EU.prj))
#writeOGR(eu, "EU_boundary_LAEA.shp", "EU_boundary_LAEA", "ESRI Shapefile") 
## Canada:
ca <- readOGR("lcsd000a15a_e.shp", "lcsd000a15a_e")
CA.prj <- proj4string(ca)

## List of property maps:
tif.lst <- list.files(path="/data/GEOG", pattern="1km_ll", full.names=TRUE, recursive=TRUE)
tif2.lst <- list.files(path="/data/GEOG", pattern="250m_ll", full.names=TRUE, recursive=TRUE)
## 220 bands!
LSTD.lst <- list.files(path="/data/MOD11A2", pattern=glob2rx("M_liste???_?.vrt"), full.names=TRUE, recursive=TRUE)
PREm.lst <- list.files(path="/data/PREm", pattern=glob2rx("PREm_*_1km_sum.sdat"), full.names=TRUE, recursive=TRUE)
sel.lst <- as.list(c(tif.lst[grep("WRB",tif.lst)], tif.lst[grep("sd",tif.lst)], LSTD.lst, PREm.lst))
## 191 layer
s1km <- lapply(sel.lst, FUN=raster)
lnames = basename(unlist(sel.lst))

## Hungary:
hun <-  spTransform(eu.ll[grep("HU", eu.ll$NUTS_ID),], CRS("+proj=longlat +ellps=WGS84"))
spplot(hun["NUTS_ID"])
## rasterize:
hun.grid <- vect2rast(hun["NUTS_ID"], cell.size=0.008333333) ## cell.size=0.002083333
hun.grid <- as(hun.grid, "SpatialPixelsDataFrame")
spplot(hun.grid)
hun.xy <- as(hun.grid, "SpatialPointsDataFrame") ## 158,919 points
str(hun.xy)
sfInit(parallel=TRUE,cpus=40)
sfLibrary(raster)
sfLibrary(rgdal)
hun.sg <- sfClusterApplyLB(s1km, extract, y=hun.xy)
sfStop()
hun.sg <- as.data.frame(hun.sg)
names(hun.sg) <- lnames
hun.grid@data <- cbind(hun.grid@data, hun.sg)
proj4string(hun.grid) = "+proj=longlat +datum=WGS84"
spplot(hun.grid["PHIHOX_M_sd4_1km_ll.tif"])

## Liaoning province:
liaoning <- spTransform(cn[cn$BOU2_4M_=="3",], CRS("+proj=longlat +ellps=WGS84"))
spplot(liaoning[1])
liaoning <- spTransform(cn[cn$BOU2_4M_=="6",], CRS("+proj=longlat +ellps=WGS84"))
spplot(liaoning[1])
liaoning.grid <- vect2rast(liaoning["BOU2_4M_ID"], cell.size=0.008333333) ## cell.size=0.002083333
liaoning.grid <- as(liaoning.grid, "SpatialPixelsDataFrame")
spplot(liaoning.grid)
liaoning.xy <- as(liaoning.grid, "SpatialPointsDataFrame") ## 221,579 points
str(liaoning.xy)
sfInit(parallel=TRUE,cpus=40)
sfLibrary(raster)
sfLibrary(rgdal)
liaoning.sg <- sfClusterApplyLB(s1km, extract, y=liaoning.xy)
sfStop()
liaoning.sg <- as.data.frame(liaoning.sg)
names(liaoning.sg) <- lnames
liaoning.grid@data <- cbind(liaoning.grid@data, liaoning.sg)
proj4string(liaoning.grid) = "+proj=longlat +datum=WGS84"
spplot(liaoning.grid["PHIHOX_M_sd4_1km_ll.tif"])

## bind EU and China:
m.grid <- rbind(data.frame(hun.grid[lnames]), data.frame(liaoning.grid[lnames]))
write.csv(m.grid, "Hungary_Liaoning_1km.csv")
gzip("Hungary_Liaoning_1km.csv")

## Fuzzy k-means clustering:
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)
m.hex <- h2o.importFile(localH2O, path = "Hungary_Liaoning_1km.csv.gz")
km20 <- h2o.kmeans(training_frame = m.hex, k = 20, x = lnames, keep_cross_validation_predictions = TRUE)
# Dropping constant columns: TAXNWRB_Haplic.Acrisols..Humic._1km_ll.tif, TAXNWRB_Haplic.Ferralsols..Xanthic._1km_ll.tif, TAXNWRB_Haplic.Fluvisols..Arenic._1km_ll.tif, TAXNWRB_Histic.Albeluvisols_1km_ll.tif, TAXNWRB_Lixic.Plinthosols_1km_ll.tif.
m.km20 <- as.data.frame(h2o.predict(km20, m.hex, na.action=na.pass))
## plot final results:
hun.grid$km20 <- as.factor(m.km20$predict[1:length(hun.grid)])
spplot(hun.grid["km20"])
plotKML(hun.grid["km20"])
liaoning.grid$km20 <- as.factor(m.km20$predict[(length(hun.grid)+1):length(m.km20$predict)])
spplot(liaoning.grid["km20"])
plotKML(liaoning.grid["km20"])
save(hun.grid, file="hun.grid.rda")
save(liaoning.grid, file="liaoning.grid.rda")
h2o.shutdown()

## resample to 1 km resolution (China):
sfInit(parallel=TRUE, cpus=40)
sfExport("tif.lst", "gdalwarp", "cn.lcc", "CN.prj")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(sel.lst, function(i){ if(!file.exists(gsub("_1km", "_1km_CN", basename(i)))){ system(paste0(gdalwarp, ' ', i, ' ', set.file.extension(gsub("_1km", "_1km_CN", basename(i)), ".tif"), ' -tr 1000 1000 -t_srs \"', CN.prj, '\" -te ', paste(as.vector(cn.lcc@bbox), collapse=" "),' -co \"COMPRESS=DEFLATE\"')) } })
sfStop()

## resample to 1 km resolution (Europe):
sfInit(parallel=TRUE, cpus=40)
sfExport("tif.lst", "gdalwarp", "eu.laea")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(sel.lst, function(i){ if(!file.exists(gsub("_1km", "_1km_EU", basename(i)))){ system(paste0(gdalwarp, ' ', i, ' ', set.file.extension(gsub("_1km", "_1km_EU", basename(i)), ".tif"), ' -t_srs \"', proj4string(eu.laea), '\" -tr 1000 1000 -te ', paste(as.vector(eu.laea@bbox), collapse=" "),' -co \"COMPRESS=DEFLATE\"')) } })
sfStop()

## Stack together LARGE RASTERS and run fuzzy k-means using all pixels:
eu1km.lst <- list.files(pattern="_EU_")
cn1km.lst <- list.files(pattern="_CN_")
#eu1km <- raster::stack(eu1km.lst)
eu1km_mask <- as.data.frame(as(raster(eu1km.lst[1]), "SpatialGridDataFrame"))
sel.eu <- as.numeric(row.names(eu1km_mask))
## 17M pixels
sfInit(parallel=TRUE, cpus=48)
sfExport("eu1km.lst", "sel.eu")
sfLibrary(rgdal)
sfLibrary(raster)
eu1km.df <- sfClusterApplyLB(eu1km.lst, function(i){ as.data.frame(raster(i))[sel.eu,] }) 
sfStop()
eu1km.df <- do.call(cbind, eu1km.df)
colnames(eu1km.df) = basename(eu1km.lst)
rownames(eu1km.df) = row.names(eu1km_mask)
eu1km.df <- cbind(eu1km_mask[,2:3], eu1km.df)
## 11GB object
names(eu1km.df) <- gsub("_EU_", "", names(eu1km.df))
eu1km.df$AREA = "EU"
save(eu1km.df, file="eu1km.df.rda")

## Same thing with China:
cn1km_mask <- as.data.frame(as(raster(cn1km.lst[1]), "SpatialGridDataFrame"))
sel.cn <- as.numeric(row.names(cn1km_mask))
## 17M pixels
sfInit(parallel=TRUE, cpus=48)
sfExport("cn1km.lst", "sel.cn")
sfLibrary(rgdal)
sfLibrary(raster)
cn1km.df <- sfClusterApplyLB(cn1km.lst, function(i){ as.data.frame(raster(i))[sel.cn,] }) 
sfStop()
cn1km.df <- do.call(cbind, cn1km.df)
colnames(cn1km.df) = basename(cn1km.lst)
rownames(cn1km.df) = row.names(cn1km_mask)
cn1km.df <- cbind(cn1km_mask[,2:3], cn1km.df)
## 10GB object
names(cn1km.df) <- gsub("_CN_", "", names(cn1km.df))
cn1km.df$AREA = "CN"
save(cn1km.df, file="cn1km.df.rda")

## Write to a MEGA CSV file:
EU_CN_1km <- rbind(eu1km.df, cn1km.df)
EU_CN_1km.s <- EU_CN_1km[sample.int(nrow(EU_CN_1km), round(0.02*nrow(EU_CN_1km))),]
l2names <- names(EU_CN_1km)[grep(".tif", names(EU_CN_1km))]
## TAKES >20mins!
save(EU_CN_1km.s, file="EU_CN_1km.s.rda", compress="xz")
#write.csv(EU_CN_1km, "EU_CN_1km.csv")
#gzip("EU_CN_1km.csv")
EU_CN_1km <- EU_CN_1km[,c("s1","s2","AREA")]
gc()

## Fuzzy k-means clustering:
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1)
m1km.hex <- h2o.importFile(path = "EU_CN_1km.csv.gz")
km45 <- h2o.kmeans(training_frame = m1km.hex[sample.int(nrow(EU_CN_1km), 5e6),], k = 45, x = l2names, keep_cross_validation_predictions = TRUE)
c.km45 <- as.data.frame(h2o.centers(km45))
save(c.km45, file="c.km45.rda")
m.km45 <- as.data.frame(h2o.predict(km45, m1km.hex, na.action=na.pass))
## plot final results:
EU_CN_1km$km45 <- m.km45$predict # as.factor(m.km45$predict)
#save(EU_CN_1km, file="km45_1km.rda")
eu1km_mask$km45 <- EU_CN_1km$km45[1:length(sel.eu)]
cn1km_mask$km45 <- EU_CN_1km$km45[(length(sel.eu)+1):nrow(EU_CN_1km)]
gridded(eu1km_mask) <- ~ s1+s2
writeGDAL(eu1km_mask["km45"], "eu1km_km45.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
gridded(cn1km_mask) <- ~ s1+s2
writeGDAL(cn1km_mask["km45"], "cn1km_km45.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
h2o.shutdown()

## resample to 1 km resolution (Canada):
sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "gdalwarp", "ca", "CA.prj")
sfLibrary(rgdal)
out <- sfClusterApplyLB(tif.lst, function(i){ if(!file.exists(gsub("1km_ll", "1km_CA", basename(i)))){ system(paste0(gdalwarp, ' ', i, ' ', gsub("1km_ll", "1km_CA", basename(i)), ' -t_srs \"', CA.prj, '\" -te ', paste(as.vector(ca@bbox), collapse=" "),' -co \"COMPRESS=DEFLATE\"')) } })
sfStop()

## LandMapper landforms:
library(R.utils)
gunzip("/data/EcoTapestry/EF_LF_Desc_250m.tif.gz")
GDALinfo("/data/tmp/can_m38_orig")
system(paste0(gdal_translate, " /data/tmp/can_m38_orig can_m38_250m.tif"))
r <- raster("can_m38_250m.tif")
r
system(paste0(gdalwarp, ' /data/EcoTapestry/EF_LF_Desc_250m.tif -r \"near\" can_EF_LF_250m.tif -co \"COMPRESS=DEFLATE\" -t_srs \"', proj4string(r), '\"  -tr 250 250 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
gzip("can_m38_250m.tif", remove=FALSE)
gzip("can_EF_LF_250m.tif", remove=FALSE)
