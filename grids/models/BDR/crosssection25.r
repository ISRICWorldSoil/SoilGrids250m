library(rgdal)
library(plotKML)
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


# function to get cross section
getcross <- function(xst, yst, cellsize, gkm, xy1, xy2, wells.depth)
{
  offxy <- c(floor((min(xy1[1],xy2[1])-xst)/cellsize),floor((yst-max(xy1[2],xy2[2]))/cellsize))
  dimxy <- ceiling(abs((xy1-xy2)/cellsize))
  #make the line and coverted it into points
  t<-readGDAL("BDTICM_M_1km_ll.tif", offset = offxy[2:1], region.dim = dimxy[2:1])
  coords <- rbind(xy1,xy2)
  sp <- SpatialPoints(coords)
  spdf <- SpatialPointsDataFrame(sp,data.frame(name = c(attr(xy1,"name"),attr(xy2,"name"))))
  L <- Line(sp)
  Ls <- Lines(list(L), ID = "a")
  SL <- SpatialLines(list(Ls))
  SLDF <- SpatialLinesDataFrame(SL, data.frame(z = 1, row.names = "a"))
  SLP <- vect2rast(SLDF, cell.size = cellsize*gkm)
  SLPDF <- as(SLP, "SpatialPixelsDataFrame")
  proj4string(SLPDF) <- proj4string(t)
  SLPDF$z <- NULL
  SLPDF$n <-1:length(SLPDF)
  spp <- as(SLPDF,"SpatialPointsDataFrame")
  #overlay
  wells.depth$n <- over(wells.depth, SLPDF)[,1]
  spo <- wells.depth@data[,c(7,6)]
  spo <- spo[!is.na(spo[,1]),]
  names(spo)[2] <- "bdt"
  spp$bdt <- over(spp,t)[,1]
  writeOGR(SLDF,t.dir,paste0(tname,"point"), driver="ESRI Shapefile", overwrite_layer =T)
  writeOGR(SLDF,t.dir,paste0(tname,"line"), driver="ESRI Shapefile", overwrite_layer =T)
  return(list(spo=spo, spp = spp@data))
}


#ohio
#set the starting point and ending point
tname <- "Ohio"
xy1 <- c(-84.606,41.584)
attr(xy1, "name")<-"Montpelier"
xy2 <- c(-82.032,39.028)
attr(xy2, "name")<-"Pomeroy"
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
cs <- getcross(xst, yst, cellsize, gkm, xy1, xy2, wells.depth)

# set the distance by lonlat or by the actual values
dis <- ((xy1[1]-xy2[1])^2 +(xy1[2]-xy2[2])^2)^0.5/cellsize
dnum <- dim(cs$spp)[1]
pat <- seq(0,dnum,100*dnum/dis)
plot(cs$spp, type = "l", xlab = paste("Distance from", attr(xy1, "name") ,"(km)") , 
     ylim= range(c(cs$spo$bdt,cs$spp$bdt),na.rm = T)[2:1], 
     ylab = "Absolute depth to bedrock (cm)", xaxt = "n")
axis(1, at = c(pat,dnum),
     labels = c(attr(xy1, "name"), (1:(length(pat)-1))*100,attr(xy2, "name")))
points(cs$spo)
dev.copy(png,paste0(tname, "_cs.png"),  width = 1000, height = 700, units = "px")
dev.off()


#Australia
#set the starting point and ending point
tname <- "Australia"
xy1 <- c(142.136,-34.212)
attr(xy1, "name")<-"Mildura"
xy2 <- c(151.202,-33.896)
attr(xy2, "name")<-"Sydney"
###read the points from Australia
setwd(ohio.dir)
us <- read.csv("wells_as2.txt", sep="\t")
us$BDRICM <- us$D_BR*100
names(us)[2:3] <- c("LONWGS84", "LATWGS84")
us <- us[!is.na(us$LONWGS84) & !is.na(us$LONWGS84),]
wells.depth <- us[us$LONWGS84>xy1[1] & us$LONWGS84< xy2[1]& us$LATWGS84>xy1[2] & us$LATWGS84<xy2[2], ]
coordinates(wells.depth) <- ~ LONWGS84+LATWGS84
proj4string(wells.depth) <- proj4string(t)    


setwd(t.dir)
gkm <- 1
cs <- getcross(xst, yst, cellsize, gkm, xy1, xy2, wells.depth)

# set the distance by lonlat or by the actual values
dis <- ((xy1[1]-xy2[1])^2 +(xy1[2]-xy2[2])^2)^0.5/cellsize
dnum <- dim(cs$spp)[1]
pat <- seq(0,dnum,100*dnum/dis)
plot(cs$spp, type = "l", xlab = paste("Distance from", attr(xy1, "name") ,"(km)") , 
     ylim= range(c(cs$spo$bdt,cs$spp$bdt),na.rm = T)[2:1], 
     ylab = "Absolute depth to bedrock (cm)", xaxt = "n")
axis(1, at = c(pat,dnum),
     labels = c(attr(xy1, "name"), (1:(length(pat)-1))*100,attr(xy2, "name")))
points(cs$spo)
dev.copy(png,paste0(tname, "_cs.png"),  width = 1000, height = 700, units = "px")
dev.off()




