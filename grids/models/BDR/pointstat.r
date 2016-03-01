library(rgdal)
library(sp)

a.dir <- "/home/shang009/big/"# dir of the project
m.dir <- paste0(a.dir, "/soildepth2")
setwd(m.dir)
ovA <- read.csv("ov.BDR_SoilGrids250m.csv")

ovB <- ovA[, c("LONWGS84","LATWGS84","BDRICM","BDTICM","BDRLOG")]
rm(ovA)

sp <- ovB[!is.na(ovB$LONWGS84)&!is.na(ovB$LATWGS84),]

coordinates(sp) <-  ~ LONWGS84+LATWGS84

c <- readOGR(dsn="Continents.shp", layer="Continents")
c@data <- as.data.frame(c$PLACENAME)
names(c) <- "p"
sp$p <- over(sp,c)[,1]
sp$p2 <- sp$p
sp$p2[sp$p=="Oceania East" |sp$p=="Oceania West"] <- "Australia"

 t1<- tapply(sp$BDRICM, sp$p2, function(x){ 
   x[x==250] <-NA
   return(cbind(min(x,na.rm=T), mean(x,na.rm=T), median(x,na.rm=T), max(x,na.rm=T),sum(!is.na(x)))) })

t1 <- do.call(rbind,t1)
x <-sp$BDRICM
t1 <- rbind(t1,cbind(min(x,na.rm=T), mean(x,na.rm=T), median(x,na.rm=T), max(x,na.rm=T),sum(!is.na(x))) )
row.names(t1) <- c("Africa","Antarctica","Asia","Oceania","Europe","North America","South America", "World")


t2<- tapply(sp$BDTICM, sp$p2, function(x){ 
  return(cbind(min(x,na.rm=T), mean(x,na.rm=T), median(x,na.rm=T), max(x,na.rm=T),sum(!is.na(x)))) })

t2 <- do.call(rbind,t2)
x <-sp$BDTICM
t2 <- rbind(t2,cbind(min(x,na.rm=T), mean(x,na.rm=T), median(x,na.rm=T), max(x,na.rm=T),sum(!is.na(x))) )
row.names(t2) <- c("Africa","Antarctica","Asia","Oceania","Europe","North America","South America", "World")

