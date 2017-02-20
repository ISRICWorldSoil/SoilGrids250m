## Mosaick function for SoilGrids in ll projection system

make_mosaick_ll <- function(varn, i, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/GEOG", ot="Int16", dstnodata=-32768, dominant=FALSE, resample="near"){
  out.tif <- paste0(out.path, "/", varn, "_", i, "_250m_ll.tif")
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -ot \"', paste(ot), '\" -dstnodata \"',  paste(dstnodata), '\" -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000'))
    system(paste0('gdaladdo ', out.tif, ' 2 4 8 16 32 64 128'))
    if(dominant==TRUE){
      system(paste0('gdal_translate -of GTiff -r \"near\" -tr 1000 1000 ', vrt.tmp, ' ', gsub("_250m.tif", "_1km.tif", out.tif), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
    } else {
      system(paste0('gdal_translate -of GTiff -r \"average\" -tr 1000 1000 ', vrt.tmp, ' ', gsub("_250m.tif", "_1km.tif", out.tif), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
    }
    unlink(vrt.tmp)
    unlink(out.tmp)
  }
}
