## predict texture classes per pixel using Silt, Sand and Clay data
## by: tom.hengl@isric.org

library(soiltexture)
library(plyr)
library(raster)

tex.c <- data.frame(class.n=TT.classes.tbl(class.sys="USDA.TT", collapse=", ")[,"abbr"], class.i=1:12)
trim <- function (x){ gsub("^\\s+|\\s+$", "", x) }

frac2TEX <- function(x){
  require(soiltexture)
  TT <- soiltexture::TT.points.in.classes(tri.data=x, class.sys="USDA.TT", PiC.type="t", tri.sum.tst=FALSE)
  ## filter transitional classes:
  no.TT <- sapply(TT, function(x){length(strsplit(x, ",")[[1]])})
  sel.TT <- which(no.TT > 1)
  TEX <- TT
  ## replace transitional classes with a random class:
  for(i in sel.TT){
    TEX[[i]] <- trim(strsplit(TT[i], ",")[[1]][ceiling(runif(1)*no.TT[i])])
  }
  return(TEX)
}

predictTEXclass <- function(i, in.path, depths=1:7){
  for(d in depths){
    outn <- paste0(in.path, "/", i, "/TEXMHT_M_sl", d, "_", i, ".tif")
    if(!file.exists(outn)){
      tex.tif.lst <- paste0(in.path, "/", i, "/", c("CLYPPT","SNDPPT","SLTPPT"), "_M_sl", d, "_", i, ".tif")
      if(all(file.exists(tex.tif.lst))){
        try( x <- rgdal::readGDAL(tex.tif.lst[1]) )
        if(!class(.Last.value)[1]=="try-error"){
          if(sum(!is.na(x$band1))>4){
            x$band2 <- rgdal::readGDAL(tex.tif.lst[2])$band1
            x$band3 <- rgdal::readGDAL(tex.tif.lst[3])$band1
            try( x0 <- as(x, "SpatialPixelsDataFrame"), silent=TRUE )
            if(!class(.Last.value)[1]=="try-error"&exists("x0")){
              names(x0@data) <- c("CLAY", "SAND", "SILT")
              sel <- complete.cases(x0@data)
              if(length(sel)>2){
                x0 <- x0[sel,]
                tex.i <- frac2TEX(x0@data)  
                tex.i <- as.vector(unlist(tex.i))
                ## convert to integers:
                x0$TEX_M <- plyr::join(data.frame(class.n=tex.i), tex.c, type="left", match="first")$class.i
                rgdal::writeGDAL(x0["TEX_M"], outn, mvFlag=255, type="Byte", options="COMPRESS=DEFLATE")
              }
            } else {
              return(i)
            }
          }
        } else {
          return(i)
        }
      }
    }
  }
}

## test it:
predictTEXclass(i="NA_060_036", in.path="/data/tt/SoilGrids250m/predicted250m", depths=1:7)

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1])

library(snowfall)
library(raster)
library(rgdal)
library(soiltexture)
library(plyr)
sfInit(parallel=TRUE, cpus=48)
sfExport("predictTEXclass", "frac2TEX", "trim", "pr.dirs", "tex.c")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(soiltexture)
sfLibrary(plyr)
out <- sfClusterApplyLB(pr.dirs, function(i){try( predictTEXclass(i, in.path="/data/tt/SoilGrids250m/predicted250m") )})
sfStop()
