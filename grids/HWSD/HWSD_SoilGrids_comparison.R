## Comparison HWSD and SoilGrids250 using WoSIS points
## Two versions of polygon-based soil property maps:
## (1) Regridded Harmonized World Soil Database v1.2 (http://dx.doi.org/10.3334/ORNLDAAC/1247)
## (2) ISRIC-WISE soil property maps at 1 km (http://www.isric.org/data/isric-wise-derived-soil-property-estimates-30-30-arcsec-global-grid-wise30sec)
## Tom.Hengl@isric.org and Maria.RuiperezGonzales@wur.nl

library(sp)
library(rgdal)
library(raster)
library(aqp)
library(GSIF)
library(plyr)
library(dplyr)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
ogr2ogr <- "/usr/local/bin/ogr2ogr"
ogrinfo <- "/usr/local/bin/ogrinfo"
source("/data/models/extract.equi7t3.R")
load("/data/models/equi7t3.rda")
#system(paste(ogrinfo, '-ro WFS:\"http://wfs.isric.org/geoserver/wosis/wfs\"'))

var.lst = c("BD","OC","ph")
var_DAAC.lst = c("BULK_DEN","OC","PH_H2O")
w.shp <- c("geoserver_bulk_density_fine_earth", "geoserver_organic_carbon", "geoserver_ph_h2o")
## download WoSIS points (focus on Bulk density soil, ph in H2O and Organic carbon)
#for(j in c(w.shp)){
#  system(paste0(ogr2ogr, ' -f \"ESRI Shapefile\" ', j, '.shp WFS:\"http://wfs.isric.org/geoserver/wosis/wfs" ', j, ' -clipsrc -180 -90 180 90'))
#}
## uncompress:
for(j in w.shp){  system(paste0('7za e ', j, '.zip -y')) }
pnt.lst <- lapply(w.shp, function(x){readOGR(paste0(x,".shp"), x)})
names(pnt.lst) = var.lst
#plot(pnt.lst[[1]])
summary(pnt.lst[["BD"]]$value) ## Bulk densities < 100 kg/cubic-m?
summary(pnt.lst[["ph"]]$value)
summary(pnt.lst[["OC"]]$value)
hist(log1p(pnt.lst[["OC"]]$value), col="grey") ## Are all "0" OC measured?

## aggregate: get WoSIS average values for depths 0-30 and 30-100 cm:
wa_depths <- function(sp, depths=c(0,30,100)){
  x <- as.data.frame(sp)
  x$depth <- x$top + (x$bottom - x$top)/2
  x$thickness <- x$bottom - x$top
  #hist(x$thickness)
  x$cl <- cut(x$depth, breaks=depths)
  x <- x[!is.na(x$cl)&x$thickness>0&!is.na(x$thickness),]
  #summary(x$cl)
  ## TH: we could also use weigthed average, but a simple average is probably also fine
  ## weigthed mean for top/sub soil:
  #df <- ddply(x, .(cl,profile_id), summarize, aggregated = sum(value*thickness)/sum(thickness))
  df <- ddply(x, .(cl,profile_id), summarize, aggregated = mean(value))
  dfxy <- join(df, as.data.frame(sp)[,c("profile_id","coords.x1","coords.x2")], type="left", match = "first")
  return(dfxy)
}
sfInit(parallel=TRUE, cpus=3)
sfExport("pnt.lst", "wa_depths")
sfLibrary(dplyr)
sfLibrary(plyr)
sfLibrary(rgdal)
out.WoSIS <- sfClusterApplyLB(pnt.lst, wa_depths)
sfStop()
names(out.WoSIS) <- var.lst
hist(out.WoSIS[["BD"]]$aggregated)
summary(out.WoSIS[["OC"]]$aggregated) ## permilles
save(out.WoSIS, file="out.WoSIS.rda")

## WoSIS average values for GlobalSoilMap depths (5):
sfInit(parallel=TRUE, cpus=3)
sfExport("pnt.lst", "wa_depths")
sfLibrary(dplyr)
sfLibrary(plyr)
sfLibrary(rgdal)
out.WoSIS_SG <- sfClusterApplyLB(pnt.lst, wa_depths, depths=c(0,5,15,30,60,100))
sfStop()
names(out.WoSIS_SG) <- var.lst
levels(out.WoSIS_SG[[1]]$cl)
save(out.WoSIS_SG, file="out.WoSIS_SG.rda")

## WoSIS average values for WISE depths (D1= 0-20 cm, D2=20-40, D3=40-60, D4=40-80, D5=80-100, D6=100-150, D7=150-200 cm):
sfInit(parallel=TRUE, cpus=3)
sfExport("pnt.lst", "wa_depths")
sfLibrary(dplyr)
sfLibrary(plyr)
sfLibrary(rgdal)
out.WoSIS_WISE <- sfClusterApplyLB(pnt.lst, wa_depths, depths=c(0,20,40,60,80,100,150,200))
sfStop()
names(out.WoSIS_WISE) <- var.lst
summary(out.WoSIS_SG[[1]]$cl)
#(0,5]   (5,15]  (15,30]  (30,60] (60,100] 
#5374    14334    12566    16636    14681 
save(out.WoSIS_WISE, file="out.WoSIS_WISE.rda")

## overlay WoSIS and HWSD:
## Rasters downloaded from http://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1247 
ts_DAAC.lst <- unlist(lapply(c("S","T"), function(x){list.files(pattern=glob2rx(paste0(x, "_*.nc4$")), full.names=TRUE, recursive=TRUE)}))
r.lst <- stack(ts_DAAC.lst)
names(r.lst) <- unlist(lapply(c("s","t"), function(x){paste(x, var.lst, sep="_")})) 
#plot(r.lst, col=SAGA_pal[[1]]) ## Many missing pixels
## convert to Geotifs:
lapply(names(r.lst), function(i){writeRaster(r.lst[[i]], filename=paste0(i, ".tif"), options=c("COMPRESS=DEFLATE"), overwrite=TRUE)})
rm(r.lst)
ts.lst <- unlist(lapply(c("s","t"), function(x){list.files(pattern=glob2rx(paste0(x, "_*.tif$")), full.names=TRUE, recursive=TRUE)}))

## Overlay and aggregate values from HWSD tifs:
extract_HWSD <- function(x, sp, OC.corf=10){
  pr <- strsplit(basename(x), "_")[[1]][1]
  varn <- strsplit(strsplit(basename(x),".tif")[[1]][1], "_")[[1]][2]
  r <- raster(x)
  names(r) <- varn
  if(pr=="s"){
    sp <- sp[sp$cl=="(30,100]",]
    d = 30+(100-30)/2
  }
  if(pr=="t"){
    sp <- sp[sp$cl=="(0,30]",]
    d = 0+(30-0)/2
  }
  coordinates(sp) <- ~ coords.x1 + coords.x2
  proj4string(sp) = proj4string(r)
  ov <- extract(r, sp, sp=TRUE)
  ov$depth = d
  ov <- as.data.frame(ov)
  if(varn=="OC"){ ov[,varn] <- ov[,varn]*OC.corf } ## permilles
  names(ov)[which(names(ov)=="aggregated")] = paste("WoSIS", varn, sep="_")
  return(ov)
}
## overlay in parallel:
sfInit(parallel=TRUE, cpus=6)
sfExport("out.WoSIS", "ts.lst", "extract_HWSD")
sfLibrary(raster)
sfLibrary(rgdal)
ovHWSD.lst <- sfClusterApplyLB(1:length(ts.lst), function(i){extract_HWSD(ts.lst[i], sp=out.WoSIS[[strsplit(strsplit(basename(ts.lst[i]),".tif")[[1]][1], "_")[[1]][2]]])})
sfStop()

## bind everything together per property:
ovHWSD.lst <- lapply(1:3, function(i){rbind.fill(ovHWSD.lst[c(i,i+3)])})
str(ovHWSD.lst)
names(ovHWSD.lst) <- var.lst
hist(ovHWSD.lst[["BD"]]$BD)
hist(ovHWSD.lst[["OC"]]$OC)
save(ovHWSD.lst, file="ovHWSD.lst.rda")

## WISE rasters (BULK, ORGC and PHAQ):
ts2.lst <- list.files(pattern=glob2rx("*_D?.tif$"), full.names=TRUE, recursive=TRUE)
in2.lst <- unlist(lapply(var.lst, rep, 7))
extract_WISE <- function(x, sp, dp=c(0,20,40,60,80,100,150,200)){
  j <- as.numeric(strsplit(strsplit(basename(x), "_D")[[1]][2], ".tif")[[1]][1])
  varn <- strsplit(basename(x),"_")[[1]][1]
  r <- raster(x)
  if(varn=="ORGC"){ names(r) <- "OC" }
  if(varn=="BULK"){ names(r) <- "BD" }
  if(varn=="PHAQ"){ names(r) <- "ph" }
  sp <- sp[sp$cl==paste0("(",dp[j],",",dp[j+1],"]"),]
  d = dp[j]+(dp[j+1]-dp[j])/2
  coordinates(sp) <- ~ coords.x1 + coords.x2
  proj4string(sp) = proj4string(r)
  ov <- extract(r, sp, sp=TRUE)
  ov$depth = d
  ov <- as.data.frame(ov)
  ## filter negative values:
  ov[,names(r)] <- ifelse(ov[,names(r)]<0, NA, ov[,names(r)])
  names(ov)[which(names(ov)=="aggregated")] = paste("WoSIS", names(r), sep="_")
  return(ov)
}

sfInit(parallel=TRUE, cpus=21)
sfExport("out.WoSIS_WISE", "ts2.lst", "in2.lst", "extract_WISE")
sfLibrary(raster)
sfLibrary(rgdal)
ovWISE.lst <- sfClusterApplyLB(1:length(ts2.lst), function(i){extract_WISE(ts2.lst[i], sp=out.WoSIS_WISE[[in2.lst[i]]])})
sfStop()
ovWISE.lst <- lapply(1:3, function(i){rbind.fill(ovWISE.lst[c((i-1)*7+1:7)])})
names(ovWISE.lst) <- var.lst
save(ovWISE.lst, file="ovWISE.lst.rda")

## overlay WoSIS with SoilGrids250m (we focus on first 5 depths)
extract_SG <- function(j, sp, depths=6, classes=c("(0,5]",  "(5,15]", "(15,30]", "(30,60]", "(60,100]", "(100,200]")){
  sp.xy <- sp[!duplicated(sp$profile_id),c("profile_id","coords.x1","coords.x2")]
  coordinates(sp.xy) <- ~ coords.x1 + coords.x2
  proj4string(sp.xy) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(j=="BD"){
    ovSG <- extract.equi7t3(x=sp.xy, y=paste0("BLD_M_sd", 1:depths), equi7t3=equi7t3, path="/data/predicted", cpus=48, fromTif=TRUE, ID="profile_id", strip.names=3)
    cf = 1/1000
  }
  if(j=="OC"){
    ovSG <- extract.equi7t3(x=sp.xy, y=paste0("ORCDRC_M_sd", 1:depths), equi7t3=equi7t3, path="/data/predicted", cpus=48, fromTif=TRUE, ID="profile_id", strip.names=3)
    cf = 1
  }
  if(j=="ph"){
    ovSG <- extract.equi7t3(x=sp.xy, y=paste0("PHIHOX_M_sd", 1:depths), equi7t3=equi7t3, path="/data/predicted", cpus=48, fromTif=TRUE, ID="profile_id", strip.names=3)
    cf = 1/10
  }
  ## bind together:
  ovSG <- data.frame(profile_id=rep(ovSG$profile_id, depths), out=unlist(ovSG[,grep("_sd", names(ovSG))])*cf, cl=unlist(lapply(classes, function(g){rep(g, length(ovSG$profile_id))})))
  names(ovSG)[2] = j
  ovSG <- plyr::join(ovSG, sp[,c("profile_id","cl","aggregated","coords.x1","coords.x2")], type="left")
  names(ovSG)[which(names(ovSG)=="aggregated")] = paste("WoSIS", j, sep="_")
  return(ovSG)
}

ovSG.lst <- list(NULL)
for(i in 1:length(var.lst)){
  ovSG.lst[[i]] <- extract_SG(j=var.lst[i], sp=out.WoSIS_SG[[var.lst[i]]])
}
names(ovSG.lst) <- var.lst
save(ovSG.lst, file="ovSG.lst.rda")
#plot(ovSG.lst[["ph"]][,c("WoSIS_ph","ph")], asp=1, xlim=c(3,9.5))
#plot(ovHWSD.lst[["ph"]][,c("WoSIS_ph","ph")], asp=1, xlim=c(3,9.5))
#plot(ovWISE.lst[["ph"]][,c("WoSIS_ph","ph")], asp=1, xlim=c(3,9.5))

## plot comparisons next to each other:
sel.plt <- list(
  data.frame(HWSD=ovHWSD.lst[["BD"]]$BD, WoSIS=ovHWSD.lst[["BD"]]$WoSIS_BD), 
  data.frame(WISE=ovWISE.lst[["BD"]]$BD, WoSIS=ovWISE.lst[["BD"]]$WoSIS_BD), 
  data.frame(SoilGrids=ovSG.lst[["BD"]]$BD, WoSIS=ovSG.lst[["BD"]]$WoSIS_BD),
  data.frame(HWSD=ovHWSD.lst[["ph"]]$ph, WoSIS=ovHWSD.lst[["ph"]]$WoSIS_ph), 
  data.frame(WISE=ovWISE.lst[["ph"]]$ph, WoSIS=ovWISE.lst[["ph"]]$WoSIS_ph), 
  data.frame(SoilGrids=ovSG.lst[["ph"]]$ph, WoSIS=ovSG.lst[["ph"]]$WoSIS_ph), 
  data.frame(HWSD=ovHWSD.lst[["OC"]]$OC, WoSIS=ovHWSD.lst[["OC"]]$WoSIS_OC), 
  data.frame(WISE=ovWISE.lst[["OC"]]$OC, WoSIS=ovWISE.lst[["OC"]]$WoSIS_OC),
  data.frame(SoilGrids=ovSG.lst[["OC"]]$OC, WoSIS=ovSG.lst[["OC"]]$WoSIS_OC))
names(sel.plt) <- unlist(lapply(c("BD fine earth t / cubic-m", "pH in water", "SOC in g/kg"), rep, 3))
xlim.lst <- c(rep(list(c(0.4,2.5)), 3), rep(list(c(3.5,9.5)), 3), rep(list(c(1,600)), 3))

pfun <- function(x,y, ...){ 
  panel.hexbinplot(x,y, ...)  
  panel.abline(0,1,lty=1,lw=2,col="black") 
}

plotList <- list(NULL)
for(i in 1:length(sel.plt)){
  varn = names(sel.plt)[i]
  pred = sel.plt[[i]][,1]
  meas = sel.plt[[i]][,2]
  xlab = names(sel.plt[[i]])[2]
  ylab = names(sel.plt[[i]])[1]
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  R.squared = signif(1-var(meas - pred, na.rm=TRUE)/var(meas, na.rm=TRUE), 2)
  main = paste0(varn, " (A.V.E.: ", R.squared, ")")
  if(!varn == "SOC in g/kg"){  
    plotList[[i]] <- hexbinplot(pred ~ meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=main, type="g", xlab=xlab, ylab=ylab, lwd=1, lcex=8, inner=.2, cex.labels=.8, xlim=xlim.lst[[i]], ylim=xlim.lst[[i]], asp=1, xbins=30, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) ## 
  } else {
    d.meas <- min(meas, na.rm=TRUE)
    pred <- pred+ifelse(d.meas==0, 1, d.meas)
    meas <- meas+ifelse(d.meas==0, 1, d.meas)
    plotList[[i]] <- hexbinplot(pred ~ meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=main, type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=30, xlim=xlim.lst[[i]], ylim=xlim.lst[[i]], xlab=xlab, ylab=ylab, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))  
  }
}

do.call(grid.arrange, c(plotList[c(1:3)], ncol=3))
do.call(grid.arrange, c(plotList[c(4:6)], ncol=3))
do.call(grid.arrange, c(plotList[c(7:9)], ncol=3))
## Conclusion: HWSD/WISE estimates of BD and SOC of critically low accuracy when compared to WoSIS points (ground truth)
