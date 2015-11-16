## Generation of DEM parameters:
library(RSAGA)
myenv <- rsaga.env(path="C:/Progra~2/SAGA-GIS")
## http://www.saga-gis.org/saga_module_doc/2.1.4/
library(R.utils)
library(snowfall)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
}

saga_TA <- function(inputFile, smaskFile, fill.gaps=TRUE){
  ## DEM file:
  filen <- strsplit(inputFile, ".gz")[[1]][1]
  if(!file.exists(set.file.extension(filen, ".sdat"))){
    gunzip(inputFile, overwrite=TRUE, remove=FALSE)
    system(paste(gdal_translate, filen, set.file.extension(filen, ".sdat"), '-ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
    unlink(filen)
    if(fill.gaps==TRUE){
      ## mask file:
      smkfilen <- strsplit(smaskFile, ".gz")[[1]][1]
      if(!file.exists(set.file.extension(smkfilen, ".sdat"))){
        gunzip(smaskFile, overwrite=TRUE, remove=FALSE)
        system(paste(gdal_translate, smkfilen, set.file.extension(smkfilen, ".sdat"), '-ot \"Int16\" -of \"SAGA\" -a_nodata 0'))
      }
      unlink(smkfilen)
      ## Fill in missing DEM pixels:
      suppressWarnings( rsaga.geoprocessor(lib="grid_tools", module=25, param=list(GRID=set.file.extension(filen, ".sgrd"), MASK=set.file.extension(smkfilen, ".sgrd"), CLOSED=set.file.extension(filen, ".sgrd")), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    }
  }
  ## Slope, curvature:
  if(!file.exists(set.file.extension(gsub("DEM", "SLP", filen), ".zip"))){
    suppressWarnings( rsaga.geoprocessor(lib="ta_morphometry", module=0, param=list(ELEVATION=set.file.extension(filen, ".sgrd"), SLOPE=set.file.extension(gsub("DEM", "SLP", filen), ".sgrd"), C_PROF=set.file.extension(gsub("DEM", "CRV", filen), ".sgrd")), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    #suppressWarnings( rsaga.slope(in.dem=set.file.extension(filen, ".sgrd"), out.slope=set.file.extension(gsub("DEM", "SLP", filen), ".sgrd"), check.module.exists = FALSE, show.output.on.console = FALSE, env=myenv) )
    slp.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "SLP", filen))
    try( zip(set.file.extension(gsub("DEM", "SLP", filen), ".zip"), files=slp.lst) )
    unlink(slp.lst)
    #suppressWarnings( rsaga.profile.curvature(in.dem=set.file.extension(filen, ".sgrd"), out.vcurv=set.file.extension(gsub("DEM", "CRV", filen), ".sgrd"), check.module.exists = FALSE, show.output.on.console = FALSE, env=myenv) )
    crv.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "CRV", filen))
    try( zip(set.file.extension(gsub("DEM", "CRV", filen), ".zip"), files=crv.lst) )
    unlink(crv.lst)
  }
  ## SAGA TWI
  if(!file.exists(set.file.extension(gsub("DEM", "TWI", filen), ".zip"))){
    ## Resample to 500 m:
    suppressWarnings( rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=set.file.extension(filen, ".sgrd"), TARGET_OUT_GRID=set.file.extension(gsub("250m", "500m", filen), ".sgrd"), SCALE_UP_METHOD=6, TARGET_DEFINITION=0, TARGET_USER_SIZE=500), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    suppressWarnings( rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM=set.file.extension(gsub("250m", "500m", filen), ".sgrd"), TWI=set.file.extension(gsub("DEM", "TWI", filen), ".sgrd")), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    twi.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "TWI", filen))
    try( zip(set.file.extension(gsub("DEM", "TWI", filen), ".zip"), files=twi.lst) )
    unlink(twi.lst)
  }
  ## MrVBF:
  if(!file.exists(set.file.extension(gsub("DEM", "VBF", filen), ".zip"))){
    suppressWarnings( rsaga.geoprocessor(lib="ta_morphometry", module=8, param=list(DEM=set.file.extension(filen, ".sgrd"), MRVBF=set.file.extension(gsub("DEM", "VBF", filen), ".sgrd"), T_SLOPE=10, P_SLOPE=3), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    vbf.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "VBF", filen))
    try( zip(set.file.extension(gsub("DEM", "VBF", filen), ".zip"), files=vbf.lst) )
    unlink(vbf.lst)
  }
  ## Valley depth:
  if(!file.exists(set.file.extension(gsub("DEM", "VDP", filen), ".zip"))){
    suppressWarnings( rsaga.geoprocessor(lib="ta_channels", module=7, param=list(ELEVATION=set.file.extension(filen, ".sgrd"), VALLEY_DEPTH=set.file.extension(gsub("DEM", "VDP", filen), ".sgrd")), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    vdp.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "VDP", filen))
    try( zip(set.file.extension(gsub("DEM", "VDP", filen), ".zip"), files=vdp.lst) )
    unlink(vdp.lst)
  }
  ## Deviation from Mean Value:
  if(!file.exists(set.file.extension(gsub("DEM", "DVM", filen), ".zip"))){
    suppressWarnings( rsaga.geoprocessor(lib="statistics_grid", module=1, param=list(GRID=set.file.extension(filen, ".sgrd"), DEVMEAN=set.file.extension(gsub("DEM", "DVM", filen), ".sgrd"), RADIUS=11), check.module.exists = FALSE, show.output.on.console = FALSE, warn=FALSE, env=myenv) )
    dvm.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "DVM", filen))
    try( zip(set.file.extension(gsub("DEM", "DVM", filen), ".zip"), files=dvm.lst) )
    unlink(dvm.lst)
  }
  ## clean up:
  dem.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=filen)
  unlink(dem.lst)
  smk.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=smkfilen)
  unlink(smk.lst)
}

#saga_TA(inputFile="DEM_AF_036_090.tif.gz", smaskFile="SMK_AF_036_090.tif.gz")
#saga_TA(inputFile="MDEM_AN_250m.tif.gz", fill.gaps=FALSE)
#sfInit(parallel=TRUE, type="SOCK", socketHosts=c("sheep1@sheep1","sheep1@sheep1","sheep2@sheep2","sheep2@sheep2","sheep3@sheep3","sheep3@sheep3","sheep3@sheep3"))
#inputFiles <- c("HU_DEM.tif.gz", "ET_DEMSRE6a_100m.tif.gz", "eu_DEM_km.tif.gz")
#sfInit(parallel=TRUE, cpus=3)
#sfLibrary(RSAGA)
#sfLibrary(R.utils)
#sfExport("inputFiles", "saga_TA", "gdal_translate")
#t <- sfLapply(inputFiles, saga_TA)
#sfStop()


