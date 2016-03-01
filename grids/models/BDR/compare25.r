#compare with raster
rm(list = ls(all = TRUE))
library(rgdal)
library(gdalUtils)
library(raster)
library(RSAGA)
library(plotKML)
library(gridExtra)
library(lattice)
library(hexbin)

data(R_pal, package = "plotKML")

w.dir <- "/home/shang009/big/soildepth/surficial"
t.dir <- "/home/shang009/big/soildepth2"
a.dir <- "/home/shang009/big/"# dir of the project
source(paste0(a.dir, "/soildepth/code/head/functions.r"))
r0<-tmp0
r1 <-tmp1
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
    p.at <- sapply(sapply(seq(0,5,0.25), function(x)10^x), log1p)
    x1 <- ceiling(r0@bbox[1,1]/100000)*100000
    x2 <- floor(r0@bbox[1,2]/100000)*100000
    y1 <- ceiling(r0@bbox[2,1]/100000)*100000
    y2 <- floor(r0@bbox[2,2]/100000)*100000
    plotList[[1]] <- spplot(r1, zcol = "band1", 
        scales=list(draw = TRUE,x=list(at=seq(x1,x2, 100000),labels=seq(x1,x2, 100000)/1000),
                    y=list(at=seq(y1,y2, 100000),labels=seq(y1,y2, 100000)/1000)), at = p.at,
        main = "Predicted", formula = as.formula(log1p(band1)~s1+s2) , colorkey = list(space = "right", height = 0.4),xlab="Easting (km)", ylab = "Northing (km)")
    plotList[[1]]$legend$right$args$key$labels$at <- sapply(sapply(seq(0,5,1), function(x)10^x), log1p)
    plotList[[1]]$legend$right$args$key$labels$labels <- c(0,10,100,1000,10000,"100000")
    plotList[[3]] <- spplot(r0, zcol = "band1", 
              scales=list(draw = TRUE,x=list(at=seq(x1,x2, 100000),labels=seq(x1,x2, 100000)/1000),
              y=list(at=seq(y1,y2, 100000),labels=seq(y1,y2, 100000)/1000)), at = p.at, 
        main = "Regional study", formula = as.formula(log1p(band1)~s1+s2),xlab="Easting (km)", ylab = "Northing (km)" )
    plotList[[3]]$legend$right$args$key$labels$at <- sapply(sapply(seq(0,5,1), function(x)10^x), log1p)
    plotList[[3]]$legend$right$args$key$labels$labels <- c(0,10,100,1000,10000,"100000")
    #plotList[[3]] <- spplot(r1, zcol =c("dif"),  xlab = "Difference")
    out <- r1$band1+1
    meas <- r0$band1+1
    plotList[[2]] <- hexbinplot(out~meas,
                 colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main= "",
                 xlab="Regional study", ylab="Predicted",
                 type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8,
                 scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)),
                 xlim=range(meas, na.rm=TRUE), ylim=range(meas, na.rm=TRUE),
                 #xlim=c(0,10000), ylim=c(0,10000),
                 asp=1, xbins=20, density=40, panel=pfun
                 , colorcut=c(0,0.002,0.01,0.03,0.07,0.15,0.25,0.5,1)
                 )    
    do.call(grid.arrange, c(plotList, ncol=3, nrow =1))
    #plot(plotList)
    dev.copy(png,paste0("../", tname, ".png"),  width = 1500, height = 700, units = "px")
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
prj <- CRS("+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
tmp0 <- reproject(tmp0, prj )
tmp1 <- reproject(tmp1, prj )
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
prj <- CRS("+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
tmp0 <- reproject(tmp0, prj )
tmp1 <- reproject(tmp1, prj )
#tmp1$band1 <- expm1(tmp1$band1)
st[[2]] <- comparedepth(tmp0, tmp1, tname)
del_unzip()
tmp0$dtb <- tmp0$band1
#kmldtb(tmp0,tname,0)


st[[1]]
st[[2]]


