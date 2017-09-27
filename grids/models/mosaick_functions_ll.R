## Mosaick function for SoilGrids in ll projection system
## by: tom.hengl@isric.org

make_mosaick_ll <- function(varn, i, in.path="/data/tt/SoilGrids250m/predicted250m", out.path="/data/GEOG", ot="Int16", dstnodata=-32768, dominant=FALSE, resample="near", metadata=NULL, aggregate=FALSE, te, tr, only.metadata=TRUE){
  if(missing(i)){
    out.tif <- paste0(out.path, "/", varn, "_250m_ll.tif")
  } else {
    out.tif <- paste0(out.path, "/", ifelse(varn=="BLD.f", "BLDFIE", ifelse(varn=="CECSUM", "CECSOL", varn)), "_", i, "_250m_ll.tif")
  }
  if(!file.exists(out.tif)){
    if(missing(i)){
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_T*.tif$")), full.names=TRUE, recursive=TRUE)
      if(!varn=="ACDWRB_M_ss"){
        n.s = length(strsplit(varn, "_")[[1]])+2
        tmp.lst <- tmp.lst[sapply(tmp.lst, function(k){length(strsplit(basename(k), "_")[[1]])})<n.s]
      }
    } else {
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_T*.tif$")), full.names=TRUE, recursive=TRUE)
    }
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -ot \"', paste(ot), '\" -dstnodata \"',  paste(dstnodata), '\" -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000 -tr ', tr, ' ', tr, ' -te ', te))
    system(paste0('gdaladdo ', out.tif, ' 2 4 8 16 32 64 128'))
    if(!is.null(metadata)){ 
      m = paste('-mo ', '\"', names(metadata), "=", as.vector(metadata), '\"', sep="", collapse = " ")
      command = paste0('gdal_edit.py ', m,' ', out.tif)
      system (command, intern=TRUE)
    }
    ## 1 km resolution:
    if(aggregate==TRUE){
      if(dominant==TRUE){
        system(paste0('gdal_translate -of GTiff -r \"near\" -tr ', 1/120, ' ', 1/120, ' ', vrt.tmp, ' ', gsub("/data/GEOG", "/data/GEOG/SoilGrids1km", gsub("_250m_ll.tif", "_1km_ll.tif", out.tif)), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
      } else {
        system(paste0('gdal_translate -of GTiff -r \"average\" -tr ', 1/120, ' ', 1/120, ' ', vrt.tmp, ' ', gsub("/data/GEOG", "/data/GEOG/SoilGrids1km", gsub("_250m_ll.tif", "_1km_ll.tif", out.tif)), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
      }
    }
    unlink(vrt.tmp)
    unlink(out.tmp)
  }
  if(!is.null(metadata)&only.metadata==TRUE){ 
    m = paste('-mo ', '\"', names(metadata), "=", as.vector(metadata), '\"', sep="", collapse = " ")
    command = paste0('gdal_edit.py ', m,' ', out.tif)
    system (command, intern=TRUE)
  }
}
