## Function - Tile and close gaps - SoilGrids250m
## Tom.Hengl@isric.org
## Uses only a single core i.e. one core one tile

## generic function to tile:
tile.sdat <- function(input, t, tr=250, t_srs = proj4string(t), resample.type, Z_SCALE, srcnodata, dstnodata, TYPE){
  tmp <- set.file.extension(tempfile(), ".sdat")
  te <- as.vector(bbox(t))
  ## Byte/Int16/UInt16
  system(paste0(gdalwarp, ' ', input, ' ', tmp, ' -t_srs \"', t_srs, '\" -tr ', tr, ' ', tr, ' -r \"', resample.type, '\" -srcnodata \"', srcnodata ,'\" -dstnodata \"', dstnodata, '\" -te ', paste(te, collapse=" "), ' -of \"SAGA\" -wm 4000')) 
  if(!Z_SCALE==1){
     ## http://saga-gis.org/saga_module_doc/2.2.0/grid_calculus_1.html
     system(paste0(saga_cmd, ' -c=1 grid_calculus 1 -GRIDS ', set.file.extension(tmp, ".sgrd"), ' -FORMULA g1*', Z_SCALE, ' -INTERPOLATION 0 -TYPE ', TYPE, ' -RESULT ', set.file.extension(tmp, ".sgrd")))
  }
  return(tmp)
}

## generic function to fill all gaps:
close.gaps <- function(inputTile, maskTile, outTile, ot, nodata, a_srs, zmin, method){
  tmp <- set.file.extension(tempfile(), ".sgrd")
  r <- raster(inputTile)
  zmin.count <- getValues(r<zmin)
  s.zmin.count <- sum(zmin.count,na.rm=TRUE)
  na.count <- getValues(is.na(r))
  na.count.mask <- sum(getValues(is.na(raster(set.file.extension(maskTile, ".sdat"))))&(zmin.count|na.count), na.rm=TRUE) 
  ## Filter only blocks which miss less than 95% of missing values? 
  if( na.count.mask>0 & !s.zmin.count==ncell(r) & na.count.mask<0.95*ncell(r) ){
    if(s.zmin.count>0){
      if(zmin>=0){
        #system(paste0(saga_cmd, ' -c=1 grid_calculus 1 -GRIDS ', set.file.extension(inputTile, ".sgrd"), ' -FORMULA=\"ifelse(g1+', abs(zmin),'<0,-99999,g1)\" -INTERPOLATION 0 -RESULT ', set.file.extension(inputTile, ".sgrd")))
      #} else {
        system(paste0(saga_cmd, ' -c=1 grid_calculus 1 -GRIDS ', set.file.extension(inputTile, ".sgrd"), ' -FORMULA=\"ifelse(g1<', zmin, ',-99999,g1)\" -INTERPOLATION 0 -RESULT ', set.file.extension(inputTile, ".sgrd")))
      }
    }
    if(method=="stepwise"){
      ## http://saga-gis.org/saga_module_doc/2.2.3/grid_tools_29.html
      system(paste0(saga_cmd, ' -c=1 grid_tools 29 -INPUT ', set.file.extension(inputTile, ".sgrd"), ' -MASK ', maskTile, ' -RESULT ', tmp, ' -INTERPOLATION=1 -PYRAMIDS=1'))
    } else {
      ## http://saga-gis.org/saga_module_doc/2.2.0/grid_tools_7.html 
      system(paste0(saga_cmd, ' -c=1 grid_tools 7 -INPUT ', set.file.extension(inputTile, ".sgrd"), ' -MASK ', maskTile, ' -RESULT ', tmp))
    }
    system(paste0(gdal_translate, ' ', set.file.extension(tmp, ".sdat"), ' ', outTile, ' -ot \"', ot, '\" -a_nodata \"', nodata, '\" -a_srs \"', a_srs, '\" -co \"COMPRESS=DEFLATE\"'))
    unlink(set.file.extension(tmp, ".sdat"))    
  } else {
    system(paste0(gdal_translate, ' ', inputTile, ' ', outTile, ' -ot \"', ot, '\" -a_nodata \"', nodata, '\" -a_srs \"', a_srs, '\" -co \"COMPRESS=DEFLATE\"'))
  }
  unlink(set.file.extension(inputTile, ".sdat"))
}

tile.tif <- function(t, input, nm, resample.type="near", ot, nodata, Z_SCALE=1, srcnodata=-32767, dstnodata=srcnodata, TYPE=7, zmin=srcnodata, close.gap=TRUE, method="stepwise"){
   ns <- strsplit(paste(t$SHORTNAME), " ")[[1]][2]
   maskTile <- paste0("/data/covs/", ns, "_", t$TILE, "/LMK_", ns, "_", t$TILE, ".sgrd")
   outTile <- set.file.extension(gsub("LMK", nm, maskTile), ".tif")
   if(!file.exists(outTile)&dir.exists(dirname(outTile))){
     outn <- tile.sdat(input, t=t, t_srs=proj4string(t), resample.type=resample.type, Z_SCALE=Z_SCALE, srcnodata=srcnodata, dstnodata=dstnodata, TYPE=TYPE)
     if(close.gap==TRUE){
       close.gaps(inputTile=outn, maskTile=maskTile, outTile=outTile, ot=ot, nodata=nodata, a_srs=proj4string(t), zmin=zmin, method=method)
     } else {
       system(paste0(gdal_translate, ' ', outn, ' ', outTile, ' -ot \"', ot, '\" -a_nodata \"', nodata, '\" -a_srs \"', proj4string(t), '\" -co \"COMPRESS=DEFLATE\"'))
     }
   }
}