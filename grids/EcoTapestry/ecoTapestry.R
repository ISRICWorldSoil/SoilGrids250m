## inputs: global EcoPhysiography map [http://blogs.esri.com/esri/esri-insider/2014/12/09/the-first-detailed-ecological-land-unitsmap-in-the-world/];
## by: Tom.Hengl@isric.org

library(RSAGA)
library(rgdal)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
system('/usr/local/bin/saga_cmd --version')
load("equi7t3.rda")

## Geology:
for(j in 1:length(equi7t3)){
  t_srs <- proj4string(equi7t3[[j]])
  te <- as.vector(bbox(equi7t3[[j]]))
  system(paste0(gdalwarp, ' EF_Lit_Des_250m.tif ', paste0("LITUSG", "_", names(equi7t3)[j], "_250m.sdat"), ' -of \"SAGA\" -r \"near\" -te ', paste(te, collapse=" "), ' -t_srs \"', t_srs, '\" -tr 250 250 -srcnodata 255 -dstnodata 255 -ot \"Byte\"'))
}

## Boundary fixes:
for(j in c(1,3,5,7)){
  if(!file.exists(paste0('LITUSG_', names(equi7t3)[j], '_250mf.sgrd'))){
    t_srs <- proj4string(equi7t3[[j]])
    te <- as.vector(bbox(equi7t3[[j]]))
    system(paste0(ogr2ogr, ' -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"', t_srs, '\" -clipdst ', paste(te, collapse=" "), ' ', names(equi7t3)[j], '_fixes.shp Lithology_fixes.shp'))
    system(paste0('/usr/local/bin/saga_cmd grid_gridding 0 -INPUT ', names(equi7t3)[j], '_fixes.shp -OUTPUT 0 -TARGET_DEFINITION 0 -TARGET_USER_XMIN ', te[1]+125,' -TARGET_USER_XMAX ', te[3]-125, ' -TARGET_USER_YMIN ', te[2]+125, ' -TARGET_USER_YMAX ', te[4]-125, ' -TARGET_USER_SIZE 250 -TARGET_OUT_GRID ', names(equi7t3)[j], '_Lithology_fixes.sgrd'))
    tmp <- c(paste0(names(equi7t3)[j], '_Lithology_fixes.sgrd'), paste0('LITUSG_', names(equi7t3)[j], '_250m.sgrd'))
    ## remove values using the mask:
    system(paste0('/usr/local/bin/saga_cmd -c=10 grid_calculus 1 -GRIDS=\"', paste(set.file.extension(unlist(tmp), ".sgrd"), collapse=";", sep=""), '\" -FORMULA=\"ifelse(g1=1,0,g2)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=1 -RESULT=LITUSG_', names(equi7t3)[j], '_250mf.sgrd'))
  }
}

## Convert to percentages so we can fill in the gaps:
for(j in 1:7){
  t_srs <- proj4string(equi7t3[[j]])
  for(k in c(1:11,13:16)){
    no <- ifelse(k<10, paste0("0",k), paste0(k))
    outn <- paste0('./equi7/L', no,'USG_', names(equi7t3)[j], '_250m.tif')
    inn <- paste0('LITUSG_', names(equi7t3)[j], '_250mf.sgrd')
    if(!file.exists(inn)){
      inn <- paste0('LITUSG_', names(equi7t3)[j], '_250m.sgrd') 
    }
    if(!file.exists(outn)){
      system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', paste(inn, collapse=";", sep=""), '\" -FORMULA=\"ifelse(g1=0,255,ifelse(g1=',k, ',100,0))\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=1 -RESULT=', set.file.extension(outn, ".sgrd")))
      system(paste0(gdal_translate, ' ', set.file.extension(outn, ".sdat"), ' ', outn, ' -ot \"Byte\" -a_nodata 255 -a_srs \"', t_srs,'\" -co \"COMPRESS=DEFLATE\"'))
      unlink(set.file.extension(outn, ".sgrd")) 
      unlink(set.file.extension(outn, ".prj")) 
      unlink(set.file.extension(outn, ".mgrd")) 
      unlink(set.file.extension(outn, ".sdat")) 
    }
  }
}

## Landforms
GDALinfo("EF_LF_Desc_250m.tif")

## 1km:
system(paste0(gdalwarp, ' EF_LF_Desc_250m.tif ELFUSG3a.tif -r \"near\" -te -180 -90 180 90 -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.008333333 0.008333333'))

## reproject to Equi7:
for(j in 1:length(equi7t3)){
  t_srs <- proj4string(equi7t3[[j]])
  te <- as.vector(bbox(equi7t3[[j]]))
  system(paste0(gdalwarp, ' EF_LF_Desc_250m.tif ', paste0("LFOUSG", "_", names(equi7t3)[j], "_250m.sdat"), ' -of \"SAGA\" -r \"near\" -te ', paste(te, collapse=" "), ' -t_srs \"', t_srs, '\" -tr 250 250 -srcnodata 255 -dstnodata 255 -ot \"Byte\"'))
}

## Convert to percentages:
for(j in 1:7){
  t_srs <- proj4string(equi7t3[[j]])
  for(k in 1:7){
    no <- ifelse(k<10, paste0("0",k), paste0(k))
    outn <- paste0('./equi7/F', no,'USG_', names(equi7t3)[j], '_250m.tif')
    inn <- paste0('LFOUSG_', names(equi7t3)[j], '_250m.sgrd')
    if(!file.exists(outn)){
      system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', paste(inn, collapse=";", sep=""), '\" -FORMULA=\"ifelse(g1=',k, ',100,0)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=1 -RESULT=', set.file.extension(outn, ".sgrd")))
      system(paste0(gdal_translate, ' ', set.file.extension(outn, ".sdat"), ' ', outn, ' -ot \"Byte\" -a_nodata 255 -a_srs \"', t_srs,'\" -co \"COMPRESS=DEFLATE\"'))
      unlink(set.file.extension(outn, ".sgrd")) 
      unlink(set.file.extension(outn, ".prj")) 
      unlink(set.file.extension(outn, ".mgrd")) 
      unlink(set.file.extension(outn, ".sdat")) 
    }
  }
}

## compress
#for(j in 1:length(equi7t3)){ gzip(paste0("LFOUSG", "_", names(equi7t3)[j], "_250m.tif")) }
  
