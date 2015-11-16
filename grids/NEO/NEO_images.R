## Download NEO images (monthly) from http://neo.sci.gsfc.nasa.gov/

library(XML)
library(RCurl)
library(snowfall)
library(raster)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
gdalinfo <- paste0(gdal.dir, "/gdalinfo.exe")





## via WMS:
var.name = c("MODAL2_M_AER_OD")
wms = "http://neowms.sci.gsfc.nasa.gov/wms/wms?"
i = 1

t1 <- newXMLNode("GDAL_WMS")
t1.v <- newXMLNode("Service", attrs= c(name="WMS"), parent=t1)
t1.s <- newXMLNode("ServerUrl", wms, parent=t1.v)
t1.l <- newXMLNode("CoverageName", var.name[i], parent=t1)
xml.out <- paste(var.name[i], ".xml", sep="")
saveXML(t1, file=xml.out)
system(paste(gdalinfo, xml.out))


for(i in c("MOD05_L2", "")
system("wget ftp://neoftp.sci.gsfc.nasa.gov/rgb/")
system("wget http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_E_SKY_WV&date=2015-09-01")

system(paste('wget --user-agent=\"Googlebot/2.1 (+http://www.googlebot.com/bot.html)\" --accept \"*.TIF\" -r \"http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_M_SKY_WV\"'))

system(paste('wget --user-agent=\"Googlebot/2.1 (+http://www.googlebot.com/bot.html)\" -r \"http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=1694760&cs=rgb&format=FLOAT.TIFF&width=3600&height=1800\"'))

for(i in c(2001, 2005, 2010, 2014)){
  dr.lst <- normalizePath(list.dirs(path=paste0("X:\\MODIS\\6\\MOD05_L2\\", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "\\\\")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=1694760&cs=rgb&format=FLOAT.TIFF&width=3600&height=1800', i ,'/', x)) ## –cut-dirs=4
  }
}
