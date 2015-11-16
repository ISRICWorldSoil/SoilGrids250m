## Filter the Groundwater map / fill-in all missing pixels
## Data from: Y. Fan et al. (2014) "Global Patterns of Groundwater Table Depth"

library(scales)
library(rgdal)
library(raster)
library(GSIF)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
source("/data/models/tile.R")
#system("7za e World_wtd_v2.7z")

## filter missing pixels by using TWI and slope maps
s = raster::stack(c("World_wtd_v2.sdat","/data/MDEM/SLPSRM3a.tif","/data/MDEM/TWISRM3a.tif"))
rp <- sampleRandom(s, size=5000, sp=TRUE)
rp.df <- data.frame(rp)
#rp.df <- rp.df[!is.na(rp.df$World_wtd_v2),]
m.WTD <- step(lm(log10(World_wtd_v2)~TWISRM3a+SLPSRM3a, rp.df))
summary(m.WTD)$adj.r.squared
## Only 15%!
plot(m.WTD$fitted.values~m.WTD$model[[1]], asp=1, cex=.8, pch=21, bg=alpha("blue", 0.6))
coefficients(m.WTD)
## clearly linear relationship
save(m.WTD, file="m.WTD.rda")

grd <- GDALinfo("World_wtd_v2.sdat")
tiles <- getSpatialTiles(grd, block.x=10, block.y=10)
## function to run correction block by block:
fix.WTD <- function(i, tiles, infile, g1, g2, g3, ng1, ng2){
  fname = paste0("./filtered/", gsub("v2.sdat", paste0("v2_",i,".tif"), infile))
  if(!file.exists(fname)){
    g <- readGDAL(infile, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)
    g@data[,ng1] <- readGDAL(g1, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)$band1
    g@data[,ng2] <- readGDAL(g2, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)$band1
    g@data[,"landmask"] <- readGDAL(g3, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)$band1
    try( { g$out <- 10^predict(m.WTD, g@data)-1 
    g$out <- ifelse(is.na(g@data[,1])|g@data[,1]<=0, ifelse(is.na(g$landmask), NA, g$out), g@data[,1]) 
    if(is.numeric(g$out)){
      writeGDAL(g["out"], fname=fname, drivername="GTiff", type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
    } )
  }
}

## test it:
fix.WTD(320, tiles=tiles, infile="World_wtd_v2.sdat", g1="/data/MDEM/TWISRM3a.tif", g2="/data/MDEM/SLPSRM3a.tif", g3="/data/MOD11A2/TMDMOD3.tif", ng1="TWISRM3a", ng2="SLPSRM3a")
plot(raster("./filtered/World_wtd_v2_320.tif"))

sfInit(parallel=TRUE, cpus=40)
sfExport("fix.WTD", "m.WTD", "tiles")
sfLibrary(rgdal)
t <- sfLapply(1:nrow(tiles), fix.WTD, tiles=tiles, infile="World_wtd_v2.sdat", g1="/data/MDEM/TWISRM3a.tif", g2="/data/MDEM/SLPSRM3a.tif", g3="/data/MOD11A2/TMDMOD3.tif", ng1="TWISRM3a", ng2="SLPSRM3a")
sfStop()  
## build a mosaick:
lstF.s <- list.files(path="./filtered", pattern=glob2rx(paste0("*.tif$")), full.names=TRUE)
txt <- paste0('list.txt')
cat(lstF.s, sep="\n", file=txt)
system(paste0(gdalbuildvrt, ' -input_file_list ', txt, ' World_wtd_v2_f.vrt'))
system(paste0(gdalwarp, ' World_wtd_v2_f.vrt World_wtd_v2_f.tif -ot \"Int16\" -tr 0.008333333 0.008333333 -te -180 -90 180 90'))
system("7za a World_wtd_v2_f.7z World_wtd_v2_f.tif")

