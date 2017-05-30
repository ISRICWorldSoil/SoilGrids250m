## ftp.soilgrids.org logs processing by Tom.Hengl@isric.org
## google (66.249.64.81) is always replicating the content
## The FTP uses the VSFTP file format something like this
## { 'date': date, 'username': username, 'ip': ip, 'path': path, 'size':
##  size, 'speed': speed }

load(".RData")
library(RCurl)
library(rgdal)
#library(rjson)
library(plyr)

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
## 15,003,109 hits - TAKES FEW MINS TO PROCESS
sg.logs = data.frame(Time=rep(NA, length(logs)), IP=NA, COMMENT=NA)
sg.logs$Time = as.POSIXct(paste(sapply(logs, function(x){strsplit(x, " \\[")[[1]][1]})), format="%a %b %d %H:%M:%S %Y")
sg.logs$IP = sapply(logs, function(x){strsplit(strsplit(x, "::ffff:")[[1]][2], "\\\", ")[[1]][1]})
sg.logs$IP = sapply(sg.logs$IP, function(x){strsplit(x, " ")[[1]][1]})
sg.logs$IP = gsub("\"", "", sg.logs$IP)
sg.logs$COMMENT = sapply(logs, function(x){strsplit(strsplit(x, "::ffff:")[[1]][2], "\\\", ")[[1]][2]})
str(sg.logs)
## Many very frequent IPs (webcrawlers?):
sum.ip = summary(as.factor(sg.logs$IP), maxsum = length(unique(sg.logs$IP)))
sum.ip = data.frame(IP=attr(sum.ip, "names"), Count=as.integer(sum.ip))
str(sum.ip)

## geocode IPs:
sg.ips = data.frame(IP=unique(sg.logs$IP))
str(sg.ips)
## 12,273 IPs (takes few minutes)
write(sg.ips, "soilgrids_ips.rds")
#q = lapply(sg.ips$IP, function(x){ rjson::fromJSON(getURI(paste0("freegeoip.net/json/", x))) })
q = lapply(sg.ips$IP, function(x){ try( data.frame(t(strsplit(getURI(paste0("freegeoip.net/csv/", x)), ",")[[1]]), stringsAsFactors = FALSE) ) })
#{"ip":"192.30.253.112","country_code":"US","country_name":"United States","region_code":"CA","region_name":"California","city":"San Francisco","zip_code":"94107","time_zone":"America/Los_Angeles","latitude":37.7697,"longitude":-122.3933,"metro_code":807}
length(q)
na.lst = sapply(q, function(i){is.data.frame(i)}) 
## Second round:
for(i in which(!na.lst)){
  q[[i]] = try( data.frame(t(strsplit(getURI(paste0("freegeoip.net/csv/", sg.ips$IP[i])), ",")[[1]]), stringsAsFactors = FALSE) )
}
na.lst = sapply(q, function(i){is.data.frame(i)}) 
summary(na.lst)
sg.ips = plyr::rbind.fill(q)
names(sg.ips) = c("IP","country_code","country_name","region_code","region_name","city","zip_code","time_zone","latitude","longitude","metro_code","X1","X2")
#sg.ips = cbind(sg.ips, plyr::rbind.fill(lapply(q, as.data.frame)))
str(sg.ips)
sg.ips$Count = join(sg.ips, sum.ip)$Count
summary(as.factor(sg.ips$country_name))
# China               United States                     Germany 
# 8449                         939                         303 
# India                 Netherlands              United Kingdom 
# 216                         168                         147 
# Brazil                      Canada                      France 
# 110                         104                          97 
# Australia                     Nigeria                    Ethiopia 
# 67                          62                          58 
# Italy                       Japan                    Pakistan 
# 56                          54                          53
save.image()
## Look up domain names per IP (VERY SLOW!):
sg.ips$DNS = paste(sapply(sg.ips$IP, function(i){strsplit(system(paste('nslookup', i), intern = TRUE)[5], "\tname = ")[[1]][2]}))

## World plot:
summary(as.numeric(sg.ips$latitude))
sg.ips$latitude = as.numeric(sg.ips$latitude)
sg.ips$longitude = as.numeric(sg.ips$longitude)
#plot(sg.ips$longitude, sg.ips$latitude)
summary(sg.ips$Count)
saveRDS(sg.ips, paste0("soilgrids_ftp_Geolocations_", gsub("-","_", Sys.Date()), ".rds"))
write.csv(sg.ips, paste0("soilgrids_ftp_Geolocations_", gsub("-","_", Sys.Date()), ".csv"))
system(paste0("7za a soilgrids_ftp_Geolocations_", gsub("-","_", Sys.Date()), ".zip soilgrids_ftp_Geolocations_", gsub("-","_", Sys.Date()), ".csv"))

library(leaflet)
library(htmlwidgets)
unlink("SoilGrids_ftp_downloads_worldmap.html")
m = leaflet(sg.ips[!is.na(sg.ips$longitude),c("longitude","latitude","Count","DNS")]) %>% addTiles() %>% addCircleMarkers(lng =~longitude, lat=~latitude, radius=~log1p(Count), popup = ~DNS, color = c('red'))
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

