#compare with raster
rm(list = ls(all = TRUE))
library(rgdal)
library(gdalUtils)
library(raster)
library(RSAGA)
library(gridExtra)

w.dir <- "/home/shang009/big/soildepth/surficial"
t.dir <- paste0("/home/shang009/big/soildepth2")
a.dir <- "/home/shang009/big/"# dir of the project
source(paste0(a.dir, "/soildepth/code/head/functions.r"))
#r0<-tmp0
#r1 <-tmp1
comparedepth <- function(r0, r1, tname)
{
    r1$band1[is.na(r0$band1)] <- NA
    r1$dif<- r1$band1-r0$band1
    r2 <- cor(r1$band1,r0$band1, use = "complete.obs")
    me <- mean(r1$dif, na.rm = T)
    rmse <- signif(sqrt(mean(r1$dif ^ 2, na.rm = T)), 3)
    r0$cl <- toclass(r0$band1, dclass)
    r1$cl <- toclass(r1$band1, dclass)
    com.c <- getcor(r0$cl, r1$cl, length(dclass$class))
    plotList <- NULL
    p.at <- seq(0,15000,1000)
    plotList[[1]] <- spplot(r1, zcol =c("band1"), at = p.at, xlab = "Predication by this sutdy")
    plotList[[2]] <- spplot(r0, zcol =c("band1"),  at = p.at, xlab = "Regional map")
    plotList[[3]] <- spplot(r1, zcol =c("dif"),  xlab = "Difference")
    do.call(grid.arrange, c(plotList, ncol=2, nrow =2))
    #plot(plotList)
    dev.copy(png,paste0("../", tname, ".png"),  width = 1000, height = 700, units = "px")
    dev.off()
    return(list(r2 = r2, me = me, rmse =rmse, com.c = com.c))
}

kmldtb <- function(tmp,tname,i)
{
    kml(tmp, colour = log1p(dtb), z.lim =c(3,10),
        raster_name = paste0(tname, i, ".png"),
        colour_scale = SAGA_pal[[1]],
        folder.name =   paste0(tname, i, ".kml" ),
        file.name = paste0(tname, i, ".kml"))
    flist <-  list.files(pattern = tname)
    file.copy(from=flist, to = paste0(w.dir, "/comparekml/"), overwrite = T)
    file.remove(flist)
}



st <- NULL
#iowa:T477
i<-1
setwd(paste0(w.dir, "/iowa"))
unzip("depth_to_bedrock.zip", exdir = "tmp")
tname <- "Iowa_dtb"
src_d <- paste0(w.dir, "/iowa/tmp/depth_to_bedrock.img") 
system(paste("gdalwarp  -tr 0.008333333 0.008333333 -r average",
             "-t_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' -overwrite",
             src_d,"./tmp/tmp0.tif"))   
tmp0 <- readGDAL("./tmp/tmp0.tif")
tmp0$band1 <- tmp0$band1*30.48

src_d <- paste0(t.dir, "/BDTICM_M_1km_ll.tif")
system(paste("gdalwarp  -te ",
             tmp0@bbox[1,1],tmp0@bbox[2,1],tmp0@bbox[1,2],tmp0@bbox[2,2]," -overwrite",
             src_d," ./tmp/tmp1.tif"))   

tmp1 <- readGDAL("./tmp/tmp1.tif")
#tmp1$band1 <- expm1(tmp1$band1)
st[[1]] <- comparedepth(tmp0, tmp1, tname)
del_unzip()
tmp0$dtb <- tmp0$band1
#kmldtb(tmp0,tname,0)


#Ohio: T442, T478
i <-2
setwd(paste0(w.dir, "/ohio"))
tname <- "Ohio_dtb"
src_d <- paste0(w.dir, "/ohio/depthb")
system(paste("gdalwarp  -tr 0.008333333 0.008333333 -r average",
             "-t_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' -overwrite",
             src_d,"./tmp0.tif"))
tmp0 <- readGDAL("./tmp0.tif")
tmp0$band1 <- tmp0$band1*30.48
tmp0$band1[tmp0$band1<0] <- 0

src_d <- paste0(t.dir, "/BDTICM_M_1km_ll.tif")
system(paste("gdalwarp  -te ",
             tmp0@bbox[1,1],tmp0@bbox[2,1],tmp0@bbox[1,2],tmp0@bbox[2,2]," -overwrite",
             src_d," ./tmp1.tif"))   

tmp1 <- readGDAL("./tmp1.tif")
#tmp1$band1 <- expm1(tmp1$band1)
st[[2]] <- comparedepth(tmp0, tmp1, tname)
del_unzip()
tmp0$dtb <- tmp0$band1
#kmldtb(tmp0,tname,0)


st[[1]]
st[[2]]


