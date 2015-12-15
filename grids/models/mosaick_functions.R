## Functions for mosaicking SoilGrids
## Tom.Hengl@isric.org

## Create mosaicks:
mosaick.equi7t3 <- function(i, j, varn, in.path, r, te, tr, ot, dstnodata, out.path='/data/GEOG/'){
  if(i=="dominant"){
    out.tif <- paste0(out.path, j, '/', varn, '_', j, '_250m_r.tif')
  } else {
    out.tif <- paste0(out.path, j, '/', varn, '_', i, '_', j, '_250m_r.tif')
  }
  if(!file.exists(out.tif)){
    if(i=="dominant"){
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    }
    if(length(tmp.lst)>0){
      out.tmp <- tempfile(fileext = ".txt")
      vrt.tmp <- tempfile(fileext = ".vrt")
      cat(tmp.lst, sep="\n", file=out.tmp)
      system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -te ', paste(te, collapse=" "),' -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 4000')) ## -co \"COMPRESS=DEFLATE\"
    }
  }
}

## Merge everything into a single mosaick
make_mosaick <- function(i, varn, ext, resample1="bilinear", resample2="average", r250m=TRUE, tr=0.002083333, r="bilinear", in.path="/data/predicted", ot="Byte", dstnodata=255, equi7t3, out.path='/data/GEOG/'){
  if(i=="dominant"){
    out.tif <- paste0(out.path, varn, "_250m_ll.tif")
  } else {
    out.tif <- paste0(out.path, varn, "_", i, "_250m_ll.tif")
  }
  if(!file.exists(out.tif)){
    ## per continent:
    x <- sapply(1:length(equi7t3), function(x){mosaick.equi7t3(j=names(equi7t3)[x], i=i, varn=varn, te=ext[[x]], tr=tr, r=r, in.path=in.path, ot=ot, dstnodata=dstnodata)})
    if(i=="dominant"){
      in.tif <- list.files(path=out.path, pattern=glob2rx(paste0(varn, "_*_250m_r.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      in.tif <- list.files(path=out.path, pattern=glob2rx(paste0(varn, "_", i, "_*_250m_r.tif$")), full.names=TRUE, recursive=TRUE)
    }
    if(length(in.tif)>0){
      out.tmp <- tempfile(fileext = ".txt")
      vrt.tmp <- tempfile(fileext = ".vrt")
      cat(in.tif[c(2,3,5,7,6,1,4)], sep="\n", file=out.tmp)
      system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
      ## relatively fast
      if(r250m == TRUE){
        system(paste0(gdal_translate, ' -of GTiff ', vrt.tmp, ' ', out.tif, ' -ot \"', ot, '\" -a_nodata \"', dstnodata, '\" -r \"', resample1, '\" -co \"COMPRESS=DEFLATE\" -co \"TILED=YES\" -co \"BLOCKXSIZE=512\" -co \"BLOCKYSIZE=512\" -co \"BIGTIFF=YES\"'))
      }
      system(paste0(gdal_translate, ' -of GTiff -r \"', resample2,'\" -tr 0.008333333 0.008333333 ', vrt.tmp, ' ', gsub("250m_ll.tif", "1km_ll.tif", out.tif), ' -ot \"', ot, '\" -a_nodata \"', dstnodata, '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
      unlink(out.tmp)
      unlink(vrt.tmp)
      unlink(in.tif)      
    }
  }
}