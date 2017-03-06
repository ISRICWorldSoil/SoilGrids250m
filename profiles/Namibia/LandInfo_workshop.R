## Photographs is GE:

library(plotKML)
library(rgdal)
setInternet2(TRUE)

## (2) PhotoOverlay (geotagged photo):
lst.jpg <- c("Soil_sampling_next_to_While_Lady_Lodge.JPG", 
"Testing_Land_Cover_app_next_to_While_Lady_Lodge.JPG", "River_Crossing_Lodge_terrace.JPG")

for(j in 1:lst.jpg){
 x <- getWikiMedia.ImageInfo(lst.jpg[j]) 
 f <- spPhoto(filename = x$url$url, exif.info = x$metadata)
 kml(f, filename=set.file.extension(lst.jpg[j], ".kml"))
}

## soil profiles:
NSS <- read.csv("NSS_topsoil_tex.csv")
names(NSS)
coordinates(NSS) <- ~ Longitude.E + Latitude.S
proj4string(NSS) <- "+proj=longlat +datum=WGS84"
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
# color only:
kml(NSS, folder.name="texture", subfolder.name="horizon1", shape = shape, colour = Textural.Class..Lab., points_names=NSS$Textural.Class..Lab.)                                     