## Functions for mosaicking SoilGrids using the EQUI7 grid system (https://github.com/TUW-GEO/Equi7Grid)
## Tom.Hengl@isric.org

load("equi7t1.rda")
load("equi7t3.rda")
tile.names <- names(equi7t1)
## Continents:
ext <- as.list(1:7)
ext[[1]] <- c(-32.42, -42.90, 64.92, 41.08) ## "AF"
ext[[2]] <- c(-66.4, -56.17, 55.44, -44.30) ## "AN"
ext[[3]] <- c(40.35, -4.67, 180, 87.37) ## "AS"
ext[[4]] <- c(-31.4, 32.2, 60.6, 82.40) ## "EU"
ext[[5]] <- c(-180, -9.71, -10.3, 83.3) ## "NA"
ext[[6]] <- c(92.68, -53.25, 180, 26.38) ## "OC"
ext[[7]] <- c(-122.85, -56.29, -16.04, 20.23) ## "SA"
names(ext) <- names(equi7t3)

## Create mosaicks:
mosaick.equi7t3 <- function(i, j, varn, in.path, r, te, tr, ot, dstnodata, out.path='/data/GEOG', compress, vrt.tmp=NULL){
  if(i=="dominant"){
    out.tif <- paste0(out.path, "/", j, '/', varn, '_', j, '_250m_r.tif')
  } else {
    out.tif <- paste0(out.path, "/", j, '/', varn, '_', i, '_', j, '_250m_r.tif')
  }
  if(!file.exists(out.tif)){
    if(is.null(vrt.tmp)){
      if(i=="dominant"){
        tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
      } else {
        tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
      }
      if(length(tmp.lst)>0){
        out.tmp <- tempfile(fileext = ".txt")
        vrt.tmp <- tempfile(fileext = ".vrt")
        cat(tmp.lst, sep="\n", file=out.tmp)
        system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
      } else {
        stop("Empty list")
      }
    }
    ## Two extra tiles for >180 degrees:
    if(j=="AS"){
      system(paste0('gdalwarp ', vrt.tmp, ' ', gsub("/AS/", "/chukotka/", out.tif), ' -t_srs \"+proj=longlat +datum=WGS84\" -overwrite -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -te -180 54 -168.3 83.3 -tr ', tr, ' ', tr, ' -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000')) ## chukotka
    }
    if(j=="OC"){
      system(paste0('gdalwarp ', vrt.tmp, ' ', gsub("/OC/", "/pacific/", out.tif), ' -t_srs \"+proj=longlat +datum=WGS84\" -overwrite -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -te -180 -62 -120 15 -tr ', tr, ' ', tr, ' -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000')) ## Pacific islands
    }
    if(compress==TRUE){
      system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -overwrite -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -te ', paste(te, collapse=" "),' -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"')) ## <-- compression takes MORE time. Maybe not necessary ?
    } else {
      system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -overwrite -r \"', r,'\" -ot \"', ot, '\" -dstnodata \"',  dstnodata, '\" -te ', paste(te, collapse=" "),' -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 2000'))
    }
  }
}

## Merge everything into a single mosaick
make_mosaick <- function(i, varn, ext, resample1="bilinear", resample2="average", r250m=TRUE, tr=0.002083333, r="bilinear", in.path="/data/predicted", ot="Byte", dstnodata=255, tile.names, out.path='/data/GEOG', compress=TRUE, build.pyramids=TRUE, vrt.tmp=NULL){
  if(i=="dominant"){
    out.tif <- paste0(out.path, "/", varn, "_250m_ll.tif")
    r = "near"
  } else {
    out.tif <- paste0(out.path, "/", varn, "_", i, "_250m_ll.tif")
  }
  if(!file.exists(out.tif)){
    ## build mosaics per continent:
    if(is.null(vrt.tmp)){ 
      x <- sapply(1:length(ext), function(x){ mosaick.equi7t3(j=tile.names[x], i=i, varn=varn, te=ext[[x]], tr=tr, r=r, in.path=in.path, ot=ot, dstnodata=dstnodata, compress=compress, out.path=out.path) })
    } else {
      x <- sapply(1:length(ext), function(x){ mosaick.equi7t3(j=tile.names[x], i=i, varn=varn, te=ext[[x]], tr=tr, r=r, in.path=in.path, ot=ot, dstnodata=dstnodata, compress=compress, out.path=out.path, vrt.tmp=vrt.tmp[x]) })
    }
    if(i=="dominant"){
      in.tif <- list.files(path=out.path, pattern=glob2rx(paste0(varn, "_*_250m_r.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      in.tif <- list.files(path=out.path, pattern=glob2rx(paste0(varn, "_", i, "_*_250m_r.tif$")), full.names=TRUE, recursive=TRUE)
    }
    if(length(in.tif)>0){
      outG.tmp <- tempfile(fileext = ".txt")
      vrtG.tmp <- tempfile(fileext = ".vrt")
      ## sort based on priority? 
      sort.lst = lapply(c("AN","chukotka","pacific","AS","OC","SA","AF","EU","NA"), function(x){grep(x, in.tif)})
      sort.lst = unlist(sort.lst[sapply(sort.lst, function(x){length(x)>0})])
      cat(in.tif[sort.lst], sep="\n", file=outG.tmp)
      system(paste0('gdalbuildvrt -input_file_list ', outG.tmp, ' ', vrtG.tmp))
      if(r250m == TRUE){
        if(build.pyramids==TRUE){
          system(paste0('gdalwarp ', vrtG.tmp, ' ', out.tif, ' -ot \"', ot, '\" -dstnodata \"', dstnodata, '\" -overwrite -r \"', resample1, '\" -co \"COMPRESS=DEFLATE\" -co \"TILED=YES\" -co \"BLOCKXSIZE=512\" -co \"BLOCKYSIZE=512\" -wm 2000 -co \"BIGTIFF=YES\"'))
          system(paste0('gdaladdo ', out.tif, ' 2 4 8 16 32 64 128'))
        } else {
          system(paste0('gdalwarp ', vrtG.tmp, ' ', out.tif, ' -ot \"', ot, '\" -dstnodata \"', dstnodata, '\" -overwrite -r \"', resample1, '\" -co \"COMPRESS=DEFLATE\" -wm 2000 -co \"BIGTIFF=YES\"'))
        }
      }
      ## gdal_translate relatively faster?
      system(paste0('gdal_translate -of GTiff -r \"', resample2,'\" -tr 0.008333333 0.008333333 ', vrtG.tmp, ' ', gsub("250m_ll.tif", "1km_ll.tif", out.tif), ' -ot \"', ot, '\" -a_nodata \"', dstnodata, '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
      unlink(outG.tmp)
      unlink(vrtG.tmp)
      unlink(in.tif)      
    }
  }
}