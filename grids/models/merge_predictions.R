## Create mosaicks from predictions (ca. 2500 EQUI7T3 tiles) - SoilGrids250m
## Tom.Hengl@isric.org

gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
gdal_translate =  "/usr/local/bin/gdal_translate"
gdal_merge.py = "/usr/local/bin/gdal_merge.py"
system("/usr/local/bin/gdal-config --version")
load("../equi7t3.rda")

## Create dirs:
x <- lapply(paste0("/data/GEOG/", names(equi7t3)), dir.create, recursive=TRUE)

## clean-up:
#del.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*_ll.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## Create mosaicks:
mosiack.equi7t3 <- function(i, j, varn, r="near", te){
  out.tif <- paste0('/data/GEOG/', j, '/', varn, '_', i, '_', j, '_250m_ll.tif')
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path="/data/predicted", pattern=glob2rx(paste0(varn, "_", i, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"', r,'\" -ot \"Byte\" -dstnodata \"255\" -te ', paste(te, collapse=" "),' -tr 0.002083333 0.002083333 -co \"COMPRESS=DEFLATE\"'))
  }
}

## Continents:
ext <- as.list(1:7)
#names(ext) <- names(equi7t3)
ext[[1]] <- c(-32.42, -42.90, 64.92, 41.08) ## "AF"
ext[[2]] <- c(-66.4, -56.17, 55.44, -44.30) ## "AN"
ext[[3]] <- c(40.35, -4.67, 180, 87.37) ## "AS"
ext[[4]] <- c(-31.4, 32.2, 60.6, 82.40) ## "EU"
ext[[5]] <- c(-180, -9.71, -10.3, 83.3) ## "NA"
ext[[6]] <- c(92.68, -53.25, 180, 26.38) ## "OC"
ext[[7]] <- c(-122.85, -56.29, -16.04, 20.23) ## "SA"

#test it
x = sapply(1:length(equi7t3), function(x){mosiack.equi7t3(j=names(equi7t3)[x], i="Fibrists", varn="TAXOUSDA", te=ext[[x]])})

## Reprojection TAKES CA 3-5 hrs (compression is most time-consuming?)
#levs <- m_TAXOUSDA$lev
levs <- gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", m_TAXNWRB$lev)))
for(x in 1:length(names(equi7t3))){
  sfInit(parallel=TRUE, cpus=40)
  sfExport("mosiack.equi7t3", "equi7t3", "x", "gdalbuildvrt", "gdalwarp", "ext", "levs")
  out <- sfClusterApplyLB(1:length(levs), fun=function(i){mosiack.equi7t3(j=names(equi7t3)[x], i=levs[i], varn="TAXNWRB", te=ext[[x]])}) ## "TAXOUSDA"
  sfStop()
}

## Merge everything into a single mosaick (relatively fast)
for(i in levs){
  out.tif <- paste0("/data/GEOG/TAXOUSDA_", i, "_250m_ll.tif")
  in.tif <- list.files(path="/data/GEOG", pattern=i, full.names=TRUE, recursive=TRUE)
  out.tmp <- tempfile(fileext = ".txt")
  cat(in.tif[c(2,3,5,7,6,1,4)], sep="\n", file=out.tmp)
  system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' tmp.vrt'))
  system(paste0(gdal_translate, " -of GTiff  tmp.vrt ", out.tif, " -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\""))
  system(paste0(gdal_translate, " -of GTiff -r \"average\" -tr 0.008333333 0.008333333 -co \"COMPRESS=DEFLATE\" tmp.vrt ", gsub("250m_ll.tif", "1km_ll.tif", out.tif), " -ot \"Byte\" -a_nodata \"255\""))
  #gzip(out.tif)
  #system(paste0(gdal_merge.py, " -o ", out.tif, " ", paste(in.tif, collapse = " "))) ## -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\"
}


