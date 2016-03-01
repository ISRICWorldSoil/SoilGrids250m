library(rgdal)
library(plotKML)
library(plyr)
t.dir <- "/home/shang009/big/soildepth2"
#t.dir <- "C:/Users/shanggv/Downloads"
ohio.dir <- "/data/shang009/big/soildepth/profs/well"
#ohio.dir <- "E:\\data\\soildata\\depth\\points\\profs\\well\\"



#the left upper cornor and cellsize
setwd(t.dir)
t<-readGDAL("BDTICM_M_1km_ll.tif", offset = c(0,0), region.dim = c(10,10))
xst <- as.double(t@grid@cellcentre.offset[1])
yst <- as.double(t@grid@cellcentre.offset[2])
cellsize <- as.double(t@grid@cellsize[1])
t2<-readGDAL("DEMSRE3a.tif", offset = c(0,0), region.dim = c(10,10))
xst2 <- as.double(t2@grid@cellcentre.offset[1])
yst2 <- as.double(t2@grid@cellcentre.offset[2])



# function to get cross section
getcross <- function(xst, yst, xst2, yst2, cellsize, gkm, xy1, xy2, wells.depth,prj)
{
  offxy <- c(floor((min(xy1[1],xy2[1])-xst)/cellsize),floor((yst-max(xy1[2],xy2[2]))/cellsize))
  offxy2 <- c(floor((min(xy1[1],xy2[1])-xst2)/cellsize),floor((yst2-max(xy1[2],xy2[2]))/cellsize)-1)
  dimxy <- ceiling(abs((xy1-xy2)/cellsize))
  #make the line and coverted it into points
  t<-readGDAL("BDTICM_M_1km_ll.tif", offset = offxy[2:1], region.dim = dimxy[2:1])
  d<-readGDAL("DEMSRE3a.tif", offset = offxy2[2:1], region.dim = dimxy[2:1])
  s <- readGDAL(paste0(tname, ".tif"))
  coords <- rbind(xy1,xy2)
  sp <- SpatialPoints(coords)
  proj4string(sp) <- proj4string(t)
  spdf <- SpatialPointsDataFrame(sp,data.frame(name = c(attr(xy1,"name"),attr(xy2,"name"))))
  L <- Line(sp)
  Ls <- Lines(list(L), ID = "a")
  SL <- SpatialLines(list(Ls))
  proj4string(SL) <- proj4string(t)
  SLDF <- SpatialLinesDataFrame(SL, data.frame(z = 1, row.names = "a"))
  SLP <- vect2rast(SLDF, cell.size = cellsize*gkm)
  SLPDF <- as(SLP, "SpatialPixelsDataFrame")

  SLPDF$z <- NULL
  SLPDF$n <-1:length(SLPDF)
  spp <- as(SLPDF,"SpatialPointsDataFrame")
  #overlay
  wells.depth$n <- over(wells.depth, SLPDF)[,1]
  wp <- subset(wells.depth, !is.na(wells.depth$n))
  spo <- wells.depth@data[,c(7,6)]
  spo <- spo[!is.na(spo[,1]),]  
  names(spo)[2] <- "bdt"
  spo<-aggregate(.~n,spo,mean)
  #over() has some bugs producing NA
  spp$bdt <- over(spp,t)[,1]
  spp$bdt2 <- over(spp,s)[,1]*30.48
  spp$dem <- over(spp,d)[,1]
  spp$bdd <- spp$dem-spp$bdt/100
  spp$bdd2 <- spp$dem-spp$bdt2/100
  spo<- join(spo,spp@data[,c("n","dem")],type="left")
  spo$bdd <- spo$dem-spo$bdt/100
  spdf <- spTransform(spdf,prj)
  writeOGR(spdf,paste0(t.dir,"/csection"),paste0(tname,"point"), driver="ESRI Shapefile", overwrite_layer =T)
  SLDF <- spTransform(SLDF,prj)
  writeOGR(SLDF,paste0(t.dir,"/csection"),paste0(tname,"line"), driver="ESRI Shapefile", overwrite_layer =T)
  wp <- spTransform(wp,prj)
  writeOGR(wp,paste0(t.dir,"/csection"),paste0(tname,"wells"), driver="ESRI Shapefile", overwrite_layer =T)
  return(list(spo=spo, spp = spp@data))
}


#ohio
#set the starting point and ending point
tname <- "Ohio"
xy1 <- c(-84.606,41.584)
attr(xy1, "name")<-"Montpelier"
xy2 <- c(-82.032,39.028)
attr(xy2, "name")<-"Pomeroy"
prj <- CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
###read the points from ohio
setwd(ohio.dir)
us <- read.csv("wells_us.txt", sep="\t")
us$BDRICM <- us$D_BR*100
names(us)[2:3] <- c("LONWGS84", "LATWGS84")
wells.depth <- us[us$LONWGS84>xy1[1] & us$LONWGS84< xy2[1]& us$LATWGS84>xy2[2] & us$LATWGS84<xy1[2], ]
coordinates(wells.depth) <- ~ LONWGS84+LATWGS84
proj4string(wells.depth) <- proj4string(t)    


# get the offset and dim, then read the block
setwd(t.dir)
gkm <- 1
cs <- getcross(xst, yst,xst2, yst2, cellsize, gkm, xy1, xy2, wells.depth,prj)

# set the distance by lonlat or by the actual values
dis <- ((xy1[1]-xy2[1])^2 +(xy1[2]-xy2[2])^2)^0.5/cellsize
dis <- 357
dnum <- sum(!is.na(cs$spp$bdt))
pat <- seq(0,dnum,100*dnum/dis)
plot(cs$spp[,c(1,4)], type = "l", xlab = paste("Distance from", attr(xy1, "name") ,"(km)") , 
     ylim=c(range(c(cs$spp$dem,cs$spp$bdd,cs$spp$bdd2,cs$spo$bdd),na.rm = T)[1], 450),  
     ylab = "Elevation (m)", xaxt = "n")
axis(1, at = c(pat,dnum),
     labels = c(attr(xy1, "name"), (1:(length(pat)-1))*100,attr(xy2, "name")))
points(cs$spo[,c(1,4)])
lines(cs$spp[,c(1,5)],col="red")
lines(cs$spp[,c(1,6)],col="blue")
dev.copy(tiff,paste0(tname, "_cs.tif"),  width = 1000, height = 400,  units = "px")
dev.off()


#Australia
#set the starting point and ending point
tname <- "Australia"
xy1 <- c(142.136,-34.212)
attr(xy1, "name")<-"Mildura"
xy2 <- c(151.202,-33.896)
attr(xy2, "name")<-"Sydney"


#Iowa
#set the starting point and ending point
tname <- "Iowa"
xy1 <- c(-96.3,43.343)
attr(xy1, "name")<-"Alvord"
xy2 <- c(-91.23,40.71)
attr(xy2, "name")<-"Wever"
prj <- CRS("+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
###read the points from Iowa
setwd(ohio.dir)
us <- read.csv("wells_us.txt", sep="\t")
us$BDRICM <- us$D_BR*100
names(us)[2:3] <- c("LONWGS84", "LATWGS84")
us <- us[!is.na(us$LONWGS84) & !is.na(us$LONWGS84),]
wells.depth <- us[us$LONWGS84>xy1[1] & us$LONWGS84< xy2[1]& us$LATWGS84>xy2[2] & us$LATWGS84<xy1[2], ]
coordinates(wells.depth) <- ~ LONWGS84+LATWGS84
proj4string(wells.depth) <- proj4string(t)    


setwd(t.dir)
gkm <- 1
cs <- getcross(xst, yst, xst2, yst2,cellsize, gkm, xy1, xy2, wells.depth,prj)

# set the distance by lonlat or by the actual values
dis <- ((xy1[1]-xy2[1])^2 +(xy1[2]-xy2[2])^2)^0.5/cellsize
dis <- 515
dnum <- sum(!is.na(cs$spp$bdt))
pat <- seq(0,dnum,100*dnum/dis)
plot(cs$spp[,c(1,4)], type = "l", xlab = paste("Distance from", attr(xy1, "name") ,"(km)") , 
     ylim=c(range(c(cs$spp$dem,cs$spp$bdd,cs$spp$bdd2,cs$spo$bdd),na.rm = T)[1], 500),  
     ylab = "Elevation (m)", xaxt = "n")
axis(1, at = c(pat[1:(length(pat)-1)],dnum),
     labels = c(attr(xy1, "name"), (1:(length(pat)-2))*100, attr(xy2, "name")))

points(cs$spo[,c(1,4)])
lines(cs$spp[,c(1,5)],col="red")
lines(cs$spp[,c(1,6)],col="blue")
dev.copy(tiff,paste0(tname, "_cs.tif"),  width = 1000, height = 300,  units = "px")
dev.off()



