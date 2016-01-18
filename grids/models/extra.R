
sfInit(parallel=TRUE, cpus=40)
sfExport("make.csv.gz", "pr.dirs")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(R.utils)
x <- sfLapply(pr.dirs, fun=function(i){ try(make.csv.gz(i, in.path="/data/covs") ) })
sfStop()


#gzip(out.tif)
#system(paste0(gdal_merge.py, " -o ", out.tif, " ", paste(in.tif, collapse = " ", "-ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\""))) ## -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\"))

## Parellel mosaicks per continent:
# for(x in 1:length(names(equi7t3))){
#   sfInit(parallel=TRUE, cpus=40)
#   sfExport("mosiack.equi7t3", "equi7t3", "x", "gdalbuildvrt", "gdalwarp", "ext", "levs")
#   out <- sfClusterApplyLB(1:length(levs), fun=function(i){mosiack.equi7t3(j=names(equi7t3)[x], i=levs[i], varn="TAXNWRB", te=ext[[x]])}) ## "TAXOUSDA"
#   sfStop()
# }