## Mosaick function for SoilGrids in ll projection system

make_mosaick_ll <- function(varn, i, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/GEOG", ot="Int16", dstnodata=-32768, dominant=FALSE, resample="near", metadata=NULL){
  if(is.null(i)){
    out.tif <- paste0(out.path, "/", varn, "_250m_ll.tif")
  } else {
    out.tif <- paste0(out.path, "/", varn, "_", i, "_250m_ll.tif")
  }
  if(!file.exists(out.tif)){
    if(is.null(i)){
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_*.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_*.tif$")), full.names=TRUE, recursive=TRUE)
    }
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -ot \"', paste(ot), '\" -dstnodata \"',  paste(dstnodata), '\" -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000'))
    system(paste0('gdaladdo ', out.tif, ' 2 4 8 16 32 64 128'))
    if(!is.null(metadata)){ 
      m = paste('-mo ', '\"', names(metadata), "=", as.vector(metadata), '\"', sep="", collapse = " ")
      command = paste0('gdal_edit.py ', m,' ', out.tif)
      system (command, intern=TRUE)
    }
    ## 1 km resolution:
    if(dominant==TRUE){
      system(paste0('gdal_translate -of GTiff -r \"near\" -tr ', 1/120, ' ', 1/120, ' ', vrt.tmp, ' ', gsub("/data/GEOG", "/data/GEOG/SoilGrids1km", gsub("_250m.tif", "_1km.tif", out.tif)), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
    } else {
      system(paste0('gdal_translate -of GTiff -r \"average\" -tr ', 1/120, ' ', 1/120, ' ', vrt.tmp, ' ', gsub("/data/GEOG", "/data/GEOG/SoilGrids1km", gsub("_250m.tif", "_1km.tif", out.tif)), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
    }
    unlink(vrt.tmp)
    unlink(out.tmp)
  }
}
