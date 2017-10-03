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

## Convert to sinusoidal projection ----
latlon2sin = function(input.file, output.file, mod.grid, tmp.dir="/data/tmp/", proj, pixsize, cleanup.files=TRUE, te, resample="near"){
  ## reproject grid in tiles:
  out.files = paste0(tmp.dir, "T", mod.grid$ID, "_", set.file.extension(basename(input.file), ".tif"))
  te.lst = apply(mod.grid@data[,1:4], 1, function(x){paste(x, collapse=" ")})
  if(missing(proj)){ proj = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" }
  sfInit(parallel=TRUE, cpus=48)
  sfExport("mod.grid", "te.lst", "proj", "out.files")
  #sfLibrary(rgdal)
  x <- sfClusterApplyLB(1:length(out.files), function(i){ invisible( system(paste0('gdalwarp ', input.file, ' ', out.files[i], ' -r \"', resample, '\" -t_srs \"', proj, '\" -tr ', pixsize, ' ', pixsize, ' -te ', te.lst[i]), show.output.on.console = FALSE, intern = TRUE) ) }) ## -co \"COMPRESS=DEFLATE\"
  sfStop()
  ## mosaic:
  tmp.lst = list.files(path=tmp.dir, pattern=basename(input.file), full.names=TRUE)
  out.tmp <- tempfile(fileext = ".txt")
  vrt.tmp <- tempfile(fileext = ".vrt")
  cat(tmp.lst, sep="\n", file=out.tmp)
  system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
  if(missing(te)){
    system(paste0('gdalwarp ', vrt.tmp, ' ', output.file, ' -ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\"'))
  } else {
    system(paste0('gdalwarp ', vrt.tmp, ' ', output.file, ' -ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\" -te ', te))
  }
  if(cleanup.files==TRUE){ unlink(out.files) }
}

## resample maps to coarser resolution ----
aggr_SG <- function(i, r, tr=1/120, tr.metric=1000, out.dir="/data/aggregated/1km/", ti="250m", tn="1km"){
  if(missing(r)){
    if(any(basename(i) %in% c("TAXOUSDA_250m_ll.tif", "TAXNWRB_250m_ll.tif", "TAXNWRB_300m_sin.tif", "TAXOUSDA_300m_sin.tif", paste0("TEXMHT_M_sl",1:7,"_250m_ll.tif")))){
      r = 'near'
    } else {
      r = 'average'
    }
  }
  if(any(basename(i) %in% c("OCSTHA_M_30cm_300m_sin.tif", "OCSTHA_M_100cm_300m_sin.tif", "OCSTHA_M_200cm_300m_sin.tif", "TAXNWRB_300m_sin.tif", "TAXOUSDA_300m_sin.tif"))){
    tr = tr.metric
    ti = "300m"
  }
  out.tif = paste0(out.dir, set.file.extension(gsub(paste0("_", ti), paste0("_", tn), basename(i)), ".tif"))
  if(!file.exists(out.tif)){
    system(paste0('gdalwarp ', i, ' ', out.tif, ' -r \"', r, '\" -tr ', tr, ' ', tr, ' -co \"COMPRESS=DEFLATE\"'))
  }
}
