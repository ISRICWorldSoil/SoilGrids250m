
rp@data[,paste0("LSTD_", i)] <- raster::extract(raster(paste0('/data/MOD11A2/M_liste', mS.lst[i], '_D.vrt')), rp)

rp.df$LSTD <- c((rp$LSTD_1+rp$LSTD_2)/2, (rp$LSTD_3+rp$LSTD_4)/2, (rp$LSTD_5+rp$LSTD_6)/2, (rp$LSTD_7+rp$LSTD_8)/2, (rp$LSTD_9+rp$LSTD_10)/2, (rp$LSTD_11+rp$LSTD_12)/2)

#rp.df$sday <- as.vector(sapply(1:6, function(x){rep(x, length(rp))}))
#rp.df <- rp.df[rp.df$SNW>210,]
#m.SNW <- step(lm(log1p(SNW)~cos(sday*50*pi/180)+LSTD+abs(y), rp.df))
abline(coef = coef(m.SNW)[c(1,2)], col="red", lwd=3)

#raster("Lat_500m.sdat")
system(paste0(gdalwarp, ' Lat_500m.sdat Lat_10km.tif -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
#system(paste0(gdalwarp, ' LSTD_M_annual_500m.sdat LSTD_M_annual_10km.tif -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' /data/MOD11A2/M_listeJan_N.vrt LSTD_M_Jan_10km.tif -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' /data/MOD11A2/M_listeAug_N.vrt LSTD_M_Aug_10km.tif -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
#system("7z a Lat_500m.7z Lat_500m.*")

#unzip("LSTD_500m.zip")
#system("7za e LSTD_500m.zip")

for(i in 1:length(m.lst)){
  system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"SNW_Ms_', m.lst[i], '_500m.sgrd\" -FORMULA=\"ifelse(g1<0,-32767,g1)\" -INTERPOLATION=0 -USE_NODATA=0 -TYPE=7 -RESULT=SNW_Ms_', m.lst[i], '_500m.sgrd'))
}


system(paste0("7za a -t7z SNW_Ms_500m.7z ", paste0("SNW_Ms_", m.lst, "_500m.*", collapse=" "), " -mmt=40"))

## Convert to GeoTiffs
for(i in 1:length(Ml)){
  inn <- set.file.extension(Ml[i], ".sgrd")
  outn <- gsub("_M_", "_Ms", inn)
  if(!file.exists(set.file.extension(outn, ".7z"))){
    system(paste0(gdal_translate, ' ', Ml[i], ' ', set.file.extension(Ml[i], ".sdat"), ' -ot \"Int16\" -a_nodata \"-32767\" -of \"SAGA\"'))
    ## http://saga-gis.org/saga_module_doc/2.2.0/grid_calculus_1.html
    system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', inn, '; SNW_mask_', mon[j], '.sgrd\" -FORMULA=\"ifelse(g1<0,g2,g1)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=4 -RESULT=', outn))
    unlink(set.file.extension(Ml[i], ".sdat"))
    unlink(set.file.extension(Ml[i], ".sgrd"))
    unlink(set.file.extension(Ml[i], ".prj"))
    unlink(set.file.extension(Ml[i], ".mgrd"))
    unlink(Ml[i])
  }
}


## fill in missing values using the model
for(j in 1:length(m.lst)){
  outn <- paste0("SNW_Ms_", m.lst[j],"_500m.sgrd")
  if(!file.exists(set.file.extension(outn, ".tif"))){
  tmp <- set.file.extension(tempfile(), ".sdat")
  system(paste0(gdalwarp, ' /data/MOD11A2/LSTD_M_', m.lst[j], '_1km.tif ', tmp, ' -ot \"Int16\" -r \"bilinear\" -srcnodata \"-32768\" -dstnodata \"-32767\" -of \"SAGA\" -tr 0.004166667 0.004166667 -te -180 -90 180 90 -multi')) ## 
  ## http://saga-gis.org/saga_module_doc/2.2.0/grid_calculus_1.html
  system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"', set.file.extension(tmp, ".sgrd"), ';Lat_500m.sgrd\" -FORMULA=\"10^(', signif(coefficients(m.SNW)[1], 4), signif(coefficients(m.SNW)[2], 4), '*g1+', signif(coefficients(m.SNW)[3], 4), '*abs(g2))-1\" -INTERPOLATION=0 -USE_NODATA=0 -TYPE=7 -RESULT=SNW_mask_', m.lst[j],'.sgrd'))  
  if(j==6|j==1){
  if(j==6){ tlat = 68 } else { tlat = 76 }
  system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"SNW_M_', m.lst[j], '_500m.sgrd;SNW_mask_', m.lst[j], '.sgrd;Lat_500m.sgrd\" -FORMULA=\"ifelse(g3>',tlat,',g2,ifelse(g1<10,g2,g1))\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=4 -RESULT=', outn))
  } else {
  system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"SNW_M_', m.lst[j], '_500m.sgrd;SNW_mask_', m.lst[j], '.sgrd\" -FORMULA=\"ifelse(g1<10,g2,g1)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=4 -RESULT=', outn))
  }
  system(paste0(gdal_translate, ' ',set.file.extension(outn, ".sdat"), ' ', set.file.extension(outn, ".tif"), ' -ot \"Int16\" -a_nodata \"-32767\" -co \"COMPRESS=DEFLATE\"'))
  }
  }
  
  system(paste0(gdalwarp, ' SNW_M_SepOct_500m.sdat SNW_M_SepOct_10km.tif -ot \"Int16\" -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\" -multi'))
  system(paste0(gdalwarp, ' SNW_mask_JanFeb.sdat SNW_mask_JanFeb_10km.tif -ot \"Int16\" -srcnodata \"-32767\" -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\" -multi'))
  system(paste0(gdalwarp, ' /data/RTMP/RtmpqrJmRC/file61b36a6eb83d.sdat LSTD_NovDec_10km.tif -ot \"Int16\" -srcnodata \"-32767\" -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\" -multi'))
  system(paste0(gdalwarp, ' Lat_500m.sdat Lat_10km.tif -tr 0.1 0.1 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\" -multi'))
  system(paste0(gdalwarp, ' Lat_500m.sdat Lat_10km.sdat -tr 0.1 0.1 -te -180 -90 180 90 -of \"SAGA\"'))
  #system("7z a SNW_mask_NovDec.7z SNW_mask_NovDec.*")
  system("7z a SNW_M_500m.7z SNW_M_*_500m.*")
  ## clean up:
  unlink(paste0("SNW_Ms_", m.lst,"_500m.sdat"))
  unlink(paste0("SNW_Ms_", m.lst,"_500m.mgrd"))
  unlink(paste0("SNW_Ms_", m.lst,"_500m.sgrd"))
  unlink(paste0("SNW_Ms_", m.lst,"_500m.prj"))
  
## function(i){try( fix.SNW(i, tiles=tiles, infile=paste0("SNW_M_", m.lst[j], "_500m.sdat"), g1=tmp, g2="Lat_500m.sdat", ng1="LSTD", ng2="y", tlat=tlat) ) } )

