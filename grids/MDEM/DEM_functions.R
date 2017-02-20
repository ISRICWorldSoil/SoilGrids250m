## Functions to derive DEM derivatives using SAGA GIS
## tom.hengl@isric.org

## Convert to LatLon system and mosaick:
equi7_latlon <- function(t, tvar, equi7t3, cellsize, ot="Int16", dstnodata="-32767"){
  s_srs = proj4string(equi7t3[[t]])
  nfile <- paste0(tvar, "_", names(equi7t3)[t], "_250m.sdat")
  ofile <- gsub("250m.sdat", "ll.sdat", nfile)
  if(!file.exists(ofile)){ 
    system(paste0(gdalwarp, ' ', nfile, ' ', ofile, ' -of \"SAGA\" -s_srs \"', s_srs, '\" -t_srs \"+proj=longlat +datum=WGS84\" -ot \"', ot, '\" -dstnodata \"', dstnodata, '\" -tr ', cellsize, ' ', cellsize, ' -te -180 -90 180 90'))
  }
}

## Derive some standard DEM variables of interest for soil mapping:
saga_DEM_derivatives <- function(INPUT, MASK=NULL, sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM","MRN","TPI"), RADIUS=c(9,13), cpus=56){
  if(pkgmaker::file_extension(INPUT)=="tif"){ 
    system(paste0(gdal_translate, ' ', INPUT, ' ', gsub(".tif", ".sdat", INPUT), ' -of \"SAGA\" -ot \"Int16\"'))
    INPUT = gsub(".tif", ".sgrd", INPUT)
  }
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Uplslope curvature:
  if(any(sel %in% "CRV")){
    if(!file.exists(gsub(".sgrd", "_downlocal.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 26 -DEM=\"', INPUT, '\" -C_DOWN_LOCAL=\"', gsub(".sgrd", "_downlocal.sgrd", INPUT), '\" -C_UP_LOCAL=\"', gsub(".sgrd", "_uplocal.sgrd", INPUT), '\" -C_UP=\"tmp.sgrd\" -C_LOCAL=\"tmp.sgrd\" -C_DOWN=\"', gsub(".sgrd", "_down.sgrd", INPUT), '\"') ) ) )
    }
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    if(!file.exists(gsub(".sgrd", "_slope.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
    }
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    if(!file.exists(gsub(".sgrd", "_vbf.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
    }
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    if(!file.exists(gsub(".sgrd", "_vdepth.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
    }
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    if(!file.exists(gsub(".sgrd", "_openp.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
    }
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    if(!file.exists(gsub(".sgrd", "_devmean.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=', RADIUS[1] ) ) )
    }
    if(!file.exists(gsub(".sgrd", "_devmean2.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean2.sgrd", INPUT), '\" -RADIUS=', RADIUS[2] ) ) )
    }
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    if(!file.exists(gsub(".sgrd", "_twi.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_hydrology 15 -DEM=\"', INPUT, '\" -SLOPE_MIN=0.04 -SLOPE_OFF=0.3 -AREA_MOD=\"', gsub(".sgrd", "_catchm.sgrd", INPUT), '\" -SLOPE_TYPE=0 -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) ) ## gsub("100", "250", INPUT)
    }
  }
  ## Melton Ruggedness Number:
  if(any(sel %in% "MRN")){
    if(!file.exists(gsub(".sgrd", "_mrn.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_hydrology 23 -DEM=\"', INPUT, '\" -AREA=\"tmp.sgrd\" -MRN=\"', gsub(".sgrd", "_mrn.sgrd", INPUT), '\" -ZMAX=\"tmp.sgrd\"' ) ) )
    }
  }
  ## TPI:
  if(any(sel %in% "TPI")){
    if(!file.exists(gsub(".sgrd", "_tpi.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 18 -DEM=\"', INPUT, '\" -STANDARD=1 -TPI=\"', gsub(".sgrd", "_tpi.sgrd", INPUT), '\" -RADIUS_MIN=0 -RADIUS_MAX=2000 -DW_WEIGHTING=3 -DW_BANDWIDTH=75' ) ) )
    }
  }
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
    slp.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "SLP", filen))
    try( zip(set.file.extension(gsub("DEM", "SLP", filen), ".zip"), files=slp.lst) )
    unlink(slp.lst)
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

## Positive and negative openess:
saga_Openess <- function(inputFile){
  filen <- strsplit(inputFile, ".gz")[[1]][1]
  if(!file.exists(set.file.extension(filen, ".sdat"))){
    gunzip(inputFile, overwrite=TRUE, remove=FALSE)
    system(paste(gdal_translate, filen, set.file.extension(filen, ".sdat"), '-ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
  }
  if(!file.exists(set.file.extension(gsub("DEM", "POS", filen), ".zip"))){
    system(paste0("/usr/local/bin/saga_cmd -c=7 -f=q ta_lighting 5 -DEM ", set.file.extension(filen, ".sgrd"), " -POS ", set.file.extension(gsub("DEM", "POS", filen), ".sgrd"), " -NEG ", set.file.extension(gsub("DEM", "NEG", filen), ".sgrd"), " -RADIUS 20000 -METHOD 0"))
    POS.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "POS", filen))
    try( zip(set.file.extension(gsub("DEM", "POS", filen), ".zip"), files=POS.lst, zip="/usr/bin/zip") )
    NEG.lst <- sapply(c(".sgrd", ".sdat", ".prj", ".mgrd"), set.file.extension, filename=gsub("DEM", "NEG", filen))
    try( zip(set.file.extension(gsub("DEM", "NEG", filen), ".zip"), files=NEG.lst, zip="/usr/bin/zip") )
  }
  unlink(filen)
}
