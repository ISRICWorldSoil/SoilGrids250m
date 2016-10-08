## ftp.soilgrids.org logs processing by Tom.Hengl@isric.org
## google (66.249.64.81) is always replicating the content
## The FTP uses the VSFTP file format something like this
## { 'date': date, 'username': username, 'ip': ip, 'path': path, 'size':
##  size, 'speed': speed }

library(RCurl)
library(rgdal)
library(rjson)
library(plyr)

load(".RData")
system("7za e vsftpd.zip")
lst.gz = list.files(pattern=".gz")
sapply(lst.gz, function(i){system(paste("7za e", i))})
unlink(lst.gz)
## Read all data to R:
lst.log = list.files(pattern=".log")
#x = scan(lst.log[1], what = "character", sep = "\n")
logs = lapply(lst.log, function(i){scan(i, what = "character", sep = "\n")})
logs = unlist(logs)
str(logs)
## about 833,014 hits - TAKES FEW MINS TO PROCESS
sg.logs = data.frame(Time=rep(NA, length(logs)), IP=NA, COMMENT=NA)
sg.logs$Time = as.POSIXct(paste(sapply(logs, function(x){strsplit(x, " \\[")[[1]][1]})), format="%a %b %d %H:%M:%S %Y")
sg.logs$IP = sapply(logs, function(x){strsplit(strsplit(x, "::ffff:")[[1]][2], "\\\", ")[[1]][1]})
sg.logs$IP = sapply(sg.logs$IP, function(x){strsplit(x, " ")[[1]][1]})
sg.logs$IP = gsub("\"", "", sg.logs$IP)
sg.logs$COMMENT = sapply(logs, function(x){strsplit(strsplit(x, "::ffff:")[[1]][2], "\\\", ")[[1]][2]})
str(sg.logs)
## Many very frequent IPs (webcrawlers):
sum.ip = summary(as.factor(sg.logs$IP), maxsum = length(unique(sg.logs$IP)))
sum.ip = data.frame(IP=attr(sum.ip, "names"), Count=as.integer(sum.ip))

## geocode IPs:
sg.ips = data.frame(IP=unique(sg.logs$IP))
str(sg.ips)
## 2387 IPs (takes few minutes)
q = lapply(sg.ips$IP, function(x){ rjson::fromJSON(getURI(paste0("freegeoip.net/json/", x))) })
length(q)
sg.ips = cbind(sg.ips, plyr::rbind.fill(lapply(q, as.data.frame)))
str(sg.ips)
sg.ips$Count = join(sg.ips, sum.ip)$Count
saveRDS(sg.ips, paste0("soilgrids_ftp_Geolocations_", gsub("-","_", Sys.Date()), ".rds"))
summary(as.factor(sg.ips$country_name))
# China     United States             India 
# 854               363                95 
# Germany       Netherlands    United Kingdom 
# 94                68                66 
# Tunisia            Brazil            France 
# 62                51                44 
# Belgium      South Africa         Australia 
# 39                37                32
## World plot:
library(leaflet)
library(htmlwidgets)
m = leaflet(sg.ips) %>% addTiles() %>% addCircleMarkers(~longitude, ~latitude, radius = ~log1p(Count), color = c('blue'))
saveWidget(m, file="SoilGrids_ftp_downloads_worldmap.html")

## Merge coords back:
sg.logs$longitude = join(sg.logs[,c("IP","Time")], sg.ips, match="first")$longitude
sg.logs$latitude = join(sg.logs[,c("IP","Time")], sg.ips, match="first")$latitude
saveRDS(sg.logs, paste0("soilgrids_ftp_logs_", gsub("-","_", Sys.Date()), ".rds"))

## Plot in Google Earth:
sg.logs.sel = sg.logs[sg.logs$IP %in% sum.ip$IP[sum.ip$Count<200],]
## 44,000 records only
library(spacetime)
library(plotKML)
coordinates(sg.logs.sel) <- c("longitude", "latitude")
proj4string(sg.logs.sel) <- CRS("+proj=longlat +datum=WGS84")
writeOGR(sg.logs.sel, "ftp_logs_soilgrids.shp", "ftp_logs_soilgrids", "ESRI Shapefile")
sg.logs.sel$Day = as.Date(sg.logs.sel$Time)
fmd_ST <- STIDF(sp=as(sg.logs.sel, "SpatialPoints"), time=sg.logs.sel$Day-0.5, data=sg.logs.sel@data[,c("Day","IP")], endTime=sg.logs.sel$Day+0.5)
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(fmd_ST, dtime = 24*3600, colour = IP, shape = shape, labels = "", file.name="ftp_logs_soilgrids.kml", folder.name="FTP logs")
system("zip -m ftp_logs_soilgrids.kmz ftp_logs_soilgrids.kml")
system("gnome-open ftp_logs_soilgrids.kmz")
