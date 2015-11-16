## Average Eight_Day_Snow_Cover 500 m
## by: Tom.Hengl@isric.org
## Example tile name: MOD10A2.A2001033.h27v12.005.2006358171335.hdf

library(utils)
library(R.utils)
library(snowfall)
library(parallel)
library(raster)
library(rgdal)
library(plyr)
library(RSAGA)
library(scales)
library(GSIF)

gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
system('/usr/local/bin/saga_cmd --version')
source("/data/models/tile.R")

## MODIS tiles and days of the year
load(file="modis_grid.rda")
str(modis_grid$TILE)
## only 317 cover land
modis_grid$TILEh <- paste0(ifelse(nchar(modis_grid$h)==1, paste0("h0", modis_grid$h), paste0("h", modis_grid$h)), ifelse(nchar(modis_grid$v)==1, paste0("v0", modis_grid$v), paste0("v", modis_grid$v)))
## 8-day period
load("dirs.rda")
days <- as.numeric(format(seq(ISOdate(2015,1,1), ISOdate(2015,12,31), by="month"), "%j"))-1
mS.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
## 6-annual periods to get smoother values - less missing values:
m.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
dsel <- cut(as.numeric(dirs), breaks=c(days[c(1,3,5,7,9,11)], 365), labels = m.lst)
m.length <- c(diff(days[c(1,3,5,7,9,11)]), 61)
## table of codes:
conv <- read.csv("MOD10A2_8day_table.csv")
## years: 2001, 2002, 2010, 2011, 2012, 2013, 2014

aggr_tile <- function(lst, outm, ml){
  if(!file.exists(outm)){
    ## TH: using raster::calc is way slower than e.g. SAGA GIS, but easier to program
    s <- raster::stack(lst)
    ## replace values:
    s <- as(s, "SpatialGridDataFrame")
    xsum <- data.frame( parallel::mclapply(s@data, function(x){conv[match(x,conv$value),"sum"]}, mc.cores=40) )
    nasum <- data.frame( parallel::mclapply(s@data, function(x){!is.na(x)}, mc.cores=40) ) ## mc.cores=ncol(s)
    ## average number of hours per month with snow:
    s$TOTAL <- (rowSums(xsum, na.rm=TRUE)/rowSums(nasum))*(ml/8)*24
    s$TOTAL <- ifelse(is.nan(s$TOTAL)|s$TOTAL==Inf, NA, s$TOTAL)
    writeGDAL(s["TOTAL"], outm, type="Int16", options="COMPRESS=DEFLATE", mvFlag="-32768")
    gc()
  }
}

## Generate average per month (per year):
aggr_parallel <- function(i){
  lstD <- list.files(path="./tiled/", pattern=glob2rx(paste0("MOD10A2.*.", modis_grid$TILEh[i], ".005.*.tif$")), full.names=TRUE)
  if(length(lstD)>0){
    for(k in 1:length(m.lst)){
      xsel <- unlist(sapply(dirs[which(m.lst[k] == dsel)], function(x){grep(basename(lstD), pattern=glob2rx(paste0("*.A2???",x,".*.*.*.tif$")))})) 
      if(length(xsel)>0){
        ## About 30 snow images per month
        try( aggr_tile(lstD[xsel], outm=paste0("./MtiledX/SNW_M_", m.lst[k], "_", modis_grid$TILEh[i], ".tif"), ml=m.length[k]) )
      }
    }
  }
}

#aggr_parallel(which(modis_grid$TILEh=="h24v05"))
aggr_parallel(1:length(modis_grid$TILEh))

# sfInit(parallel=TRUE, cpus=4)
# sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs", "conv", "m.length")
# sfLibrary(raster)
# sfLibrary(rgdal)
# sfLibrary(plyr)
# sfLibrary(sp)
# t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel)
# sfStop()

## Mosaic:
M.lstD <- list.files(path="./MtiledX/", pattern="SNW_M", full.names = TRUE)
## total number of tiles = 1738 (1738 / 6 = 290)

make_mosaic <- function(j, mon){
  #outm <- paste0("SNW_M_", mon[j], "_500m.tif")
  outm <- paste0("SNW_M_", mon[j], "_500m.sdat")
  if(!file.exists(paste0(outm,'.gz'))){
    lst.s <- M.lstD[grep(pattern=mon[j], M.lstD)]
    vrt <- paste0('M_liste', mon[j], '.vrt')
    txt <- paste0('M_liste', mon[j], '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -srcnodata \"-32768\" -of \"SAGA\" -ot \"Int16\" -tr 0.004166667 0.004166667 -te -180 -90 180 90'))
    #gzip(outm)
    #unlink(vrt)
    #unlink(txt)
  }
}

sfInit(parallel=TRUE, cpus=6)
sfExport("make_mosaic", "M.lstD", "m.lst", "gdalbuildvrt", "gdalwarp")
sfLibrary(R.utils)
sfLibrary(rgdal)
t <- sfLapply(1:6, make_mosaic, mon=m.lst)
sfStop()


## Filter all missing values using the temp map
## First, fit a model to fill in all missing values:
Ml <- paste0("SNW_M_", m.lst, "_500m.sdat")
#sapply(paste0(Ml,".gz"), gunzip)
rp <- sampleRandom(stack(Ml), size=5000, sp=TRUE)
for(i in 1:length(m.lst)){
  rp@data[,paste0("LSTD_", i)] <- raster::extract(raster(paste0('/data/MOD11A2/LSTD_M_', m.lst[i], '_1km.tif')), rp)
}
rp.df <- do.call("rbind", replicate(6, data.frame(rp)[,c("x","y")], simplify = FALSE))
rp.df$SNW <- as.vector(unlist(rp@data[,1:6]))
rp.df$LSTD <- c(rp$LSTD_1, rp$LSTD_2, rp$LSTD_3, rp$LSTD_4, rp$LSTD_5, rp$LSTD_6) 
rp.df <- rp.df[rp.df$SNW>0,]
## model SNW probability as a function of mean daily temp and latitude
m.SNW <- step(lm(log10(SNW+1)~LSTD+abs(y), rp.df))
summary(m.SNW)$adj.r.squared
## 57% variation explained by TEMP
plot(m.SNW$fitted.values~m.SNW$model[[1]], asp=1, cex=.8, pch=21, bg=alpha("blue", 0.6), xlab="Measured SNW", ylab="SNW ~ TMP + Latitude")
coefficients(m.SNW)
save(m.SNW, file="m.SNW.rda")
plot(log10(SNW+1)~LSTD, rp.df, cex=.8, pch=21, bg=alpha("blue", 0.6))
points(x=m.SNW$model$LSTD, y=m.SNW$fitted.values, pch=21, cex=.8, bg=alpha("red", 0.6))
#plot(log1p(SNW)~abs(y), rp.df)
## test the model:
10^predict(m.SNW, data.frame(LSTD=310,y=0))-1

## create latitude grid:
if(!file.exists("Lat_500m.sdat")){
  system(paste0(gdalwarp, ' /data/MOD11A2/TMDMOD3.tif LSTD_M_annual_500m.sdat -ot \"Int16\" -srcnodata \"-32767\" -of \"SAGA\" -tr 0.004166667 0.004166667 -te -180 -90 180 90'))
  system(paste0('/usr/local/bin/saga_cmd -c=40 pj_proj4 17 -GRID=\"LSTD_M_annual_500m.sgrd\" -LAT=\"Lat_500m.sgrd\"'))
}

grd <- GDALinfo("Lat_500m.sdat")
tiles <- getSpatialTiles(grd, block.x=10, block.y=10)
## function to run correction block by block:
fix.SNW <- function(i, tiles, infile, g1, g2, ng1, ng2, tlat){
  fname = paste0("./filtered/", gsub("500m.sdat", paste0("500m_",i,".tif"), infile))
  if(!file.exists(fname)){
    g <- readGDAL(infile, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)
    g@data[,ng1] <- readGDAL(g1, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)$band1
    g@data[,ng2] <- readGDAL(g2, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)$band1
    g$out <- 10^predict(m.SNW, g@data)-1
    if(tlat==90){
      g$out <- ifelse(is.na(g@data[,1])|g@data[,1]<=0, g$out, g@data[,1])
    } else {
      g$out <- ifelse(g@data[,ng2]>tlat, g$out, ifelse(is.na(g@data[,1])|g@data[,1]<0, g$out, g@data[,1]))
    }
    writeGDAL(g["out"], fname=fname, drivername="GTiff", type="Int16", mvFlag=-32767, options="COMPRESS=DEFLATE")  
  }
}

## test it:
fix.SNW(i=580, tiles, infile="SNW_M_JanFeb_500m.sdat", g1=tmp, g2="Lat_500m.sdat", ng1="LSTD", ng2="y", tlat=76)

## Create filtered SNW maps:
for(j in 1:length(m.lst)){
  vrt <- paste0('Ms_liste', m.lst[j], '.vrt')
  outn <- paste0("SNW_Ms_", m.lst[j],"_500m.tif")
  if(!file.exists(vrt)){
    tmp <- set.file.extension(tempfile(), ".sdat")
    system(paste0(gdalwarp, ' /data/MOD11A2/LSTD_M_', m.lst[j], '_1km.tif ', tmp, ' -ot \"Int16\" -r \"bilinear\" -srcnodata \"-32768\" -dstnodata \"-32767\" -of \"SAGA\" -tr 0.004166667 0.004166667 -te -180 -90 180 90 -multi'))
    ## filter by tile:
    if(j==6|j==1){
      if(j==6){ tlat = 68 } else { tlat = 76 }
    } else {
      tlat = 90
    }  
    sfInit(parallel=TRUE, cpus=40)
    sfExport("fix.SNW", "tmp", "m.lst", "m.SNW", "tiles", "tlat")
    sfLibrary(rgdal)
    t <- sfLapply(1:nrow(tiles), fix.SNW, tiles=tiles, infile=paste0("SNW_M_", m.lst[j], "_500m.sdat"), g1=tmp, g2="Lat_500m.sdat", ng1="LSTD", ng2="y", tlat=tlat)
    sfStop()  
    ## build a mosaick:
    lstF.s <- list.files(path="./filtered", pattern=glob2rx(paste0("SNW_M_", m.lst[j],"_500m_*.tif$")), full.names=TRUE)
    txt <- paste0('Ms_liste', m.lst[j], '.txt')
    cat(lstF.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    #system(paste0(gdalwarp, ' ', vrt, ' ', gsub("500m", "10km", outn), ' -ot \"Int16\" -tr 0.1 0.1 -te -180 -90 180 90'))
    #system(paste0(gdalwarp, ' ', vrt, ' ', outn, ' -r \"near\" -ot \"Int16\" -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"')) ## 
  }
}

