## Mosaick MODIS Flood map (http://oas.gsfc.nasa.gov/floodmap/)
## T. Hengl (tom.hengl@isric.org)

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp <- "/usr/local/bin/gdalwarp"
gdalbuildvrt <- "/usr/local/bin/gdalbuildvrt"
tiles <- read.csv("tiles.csv")
#m.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
m.lst <- c("DecJan","FebMar","AprMay","JunJul","AugSep","OctNov")
load("/data/EcoTapestry/equi7t3.rda")
## define functions:
meanf <- function(x){calc(x, mean, na.rm=TRUE)}

## Mosaics (per continent):  
for(i in 6:9){  
  #outn <- paste0('MFW_2015_', m.lst[i], '_250m.tif')
  #if(!file.exists(outn)){
  tmp.lst <- list.files(path="./tiledMFW/", pattern=glob2rx(paste0("MFW_2015_", m.lst[i], "_*_250m.tif$")), full.names=TRUE)
  if(length(tmp.lst)>195){
    unlink("my_liste.txt")
    cat(tmp.lst, sep="\n", file="my_liste.txt")
    gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0("tmp_", i, ".vrt"))
    #system(paste0(gdalwarp, ' tmp.vrt ', outn, ' -r \"bilinear\" -te -180 -90 180 90 -tr 0.002083333 0.002083333 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
    for(j in 1:length(equi7t3)){
      rout <- paste0('MFW_2015_', m.lst[i], '_', names(equi7t3)[j], "_250m.tif")
      if(!file.exists(rout)){
        t_srs = proj4string(equi7t3[[j]])
        te <- as.vector(bbox(equi7t3[[j]]))
        system(paste0(gdalwarp, ' tmp_', i, '.vrt ', rout, ' -r \"bilinear\" -t_srs \"', t_srs, '\" -tr 250 250 -te ', paste(te, collapse=" "), ' -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
      }
    }
  }
  #}
}

## Boundary fixes (big artifacts in "NA", "EU" and "AS"):
for(j in 1:length(equi7t3)){
  for(i in c(4,5)){ ## for(i in 1:length(m.lst)){ 
    outn <- paste0("FW", i, "MOD5_", names(equi7t3)[j], ".tif")
    if(!file.exists(outn)){
      m1 <- paste0("MFW_2015_", substr(m.lst[i], 1, 3), "_", names(equi7t3)[j], "_250m.tif")
      m2 <- paste0("MFW_2015_", substr(m.lst[i], 4, 6), "_", names(equi7t3)[j], "_250m.tif")
      m0 <- paste0("MFW_2015_", names(equi7t3)[j], "_mask.tif")
      mc <- paste0("MFW_2015_Sep_", names(equi7t3)[j], "_250m.tif")
      if(j %in% c(3:5)){
        for(k in c(m1,m2,m0,mc)){
          system(paste(gdal_translate, k, gsub(".tif", ".sdat", k), '-a_nodata 255 -ot \"Byte\" -of \"SAGA\"'))
        }
        system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', paste(gsub(".tif", ".sgrd", c(m1,m0,mc)), collapse=";", sep=""), '\" -FORMULA=\"ifelse(g2=1,g3,g1)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=7 -RESULT=\"tmp1.sgrd\"'))
        system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', paste(gsub(".tif", ".sgrd", c(m2,m0,mc)), collapse=";", sep=""), '\" -FORMULA=\"ifelse(g2=1,g3,g1)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=7 -RESULT=\"tmp2.sgrd\"'))
        s <- raster::stack(c("tmp1.sdat","tmp2.sdat"))
      } else {
        ## run in parallel:
        s <- raster::stack(c(m1,m2))
      }
      beginCluster()
      r <- clusterR(s, fun=meanf, filename=outn, datatype="INT1U", options=c("COMPRESS=DEFLATE"), NAflag=255)
      endCluster()
    }
    unlink("*.prj")
    unlink("*.sgrd")
    unlink("*.mgrd")
    unlink("*.sdat")
  }  
}

save.image("/data/floods/.RData")      
## end of script;