## Comparison HWSD and SoilGrids250 using WoSIS points
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
system(paste(ogrinfo, '-ro WFS:\"http://wfs.isric.org/geoserver/wosis/wfs\"'))

var.lst = c("BD","OC","ph")
w.shp <- c("geoserver_bulk_density_fine_earth", "geoserver_organic_carbon", "geoserver_ph_h2o")
## download WoSIS points (focus on Bulk density soil, ph in H2O and Organic carbon)
#for(j in c(w.shp)){
#  system(paste0(ogr2ogr, ' -f \"ESRI Shapefile\" ', j, '.shp WFS:\"http://wfs.isric.org/geoserver/wosis/wfs" ', j, ' -clipsrc -180 -90 180 90'))
#}
#for(j in w.shp){  system(paste0('7za e ', j, '.zip -y')) }
pnt.lst <- lapply(w.shp, function(x){readOGR(paste0(x,".shp"), x)})
names(pnt.lst) = var.lst
#plot(pnt.lst[[1]])

## aggregate: get WoSIS average values for depths 0-30 and 30-100 cm:
wa_depths <- function(sp, depths=c(0,30,100)){
  x <- as.data.frame(sp)
  x$depth <- x$top + (x$bottom - x$top)/2
  x$thickness <- x$bottom - x$top
  #hist(x$thickness)
  x$cl <- cut(x$depth, breaks=depths)
  x <- x[!is.na(x$cl)&x$thickness>0&!is.na(x$thickness),]
  #summary(x$cl)
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
levels(out.WoSIS_SG[[1]]$cl)

## WoSIS average values for GlobalSoilMap depths:
sfInit(parallel=TRUE, cpus=3)
sfExport("pnt.lst", "wa_depths")
sfLibrary(dplyr)
sfLibrary(plyr)
sfLibrary(rgdal)
out.WoSIS_SG <- sfClusterApplyLB(pnt.lst, wa_depths, depths=c(0,5,15,30,60,100))
sfStop()
names(out.WoSIS_SG) <- var.lst
save(out.WoSIS_SG, file="out.WoSIS_SG.rda")

## overlay WoSIS HWSD:
hwsd.lst <- sapply(c("s","t"), function(x){paste0(x, "_", var.lst)})
#for(j in hwsd.lst){  try( system(paste0('7za e ', j, '.7z -y')) ) }
ts.lst <- unlist(lapply(c("s","t"), function(x){list.files(pattern=glob2rx(paste0(x, "_*.tif$")), full.names=TRUE, recursive=TRUE)}))

extract_tif <- function(x, sp){
  pr <- strsplit(basename(x), "_")[[1]][1]
  r <- raster(x)
  varn <- strsplit(strsplit(basename(x),".tif")[[1]][1], "_")[[1]][2]
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
  if(varn=="OC"){ ov[,varn] <- ov[,varn]*10 } ## permilles
  names(ov)[which(names(ov)=="aggregated")] = paste("WoSIS", varn, sep="_")
  return(ov)
}
## overlay in parallel:
sfInit(parallel=TRUE, cpus=6)
sfExport("out.WoSIS", "ts.lst", "extract_tif")
sfLibrary(raster)
sfLibrary(rgdal)
ovHWSD.lst <- sfClusterApplyLB(1:length(ts.lst), function(i){extract_tif(ts.lst[i], sp=out.WoSIS[[strsplit(strsplit(basename(ts.lst[i]),".tif")[[1]][1], "_")[[1]][2]]])})
sfStop()
## bind everything together per property:
ovHWSD.lst <- lapply(1:3, function(i){rbind.fill(ovHWSD.lst[c(i,i+3)])})
str(ovHWSD.lst)
names(ovHWSD.lst) <- var.lst
hist(ovHWSD.lst[["BD"]]$BD)
save(ovHWSD.lst, file="ovHWSD.lst.rda")

## overlay with SoilGrids250m (first 5 depths)
extract_SG <- function(j, sp, depths=5, classes=c("(0,5]",  "(5,15]", "(15,30]", "(30,60]", "(60,100]")){
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

## WISE 30arcmin:
#system('7za e wise30sec_fin1.7z -y')

## plot comparisons next to each other:
sel.plt <- list(
  data.frame(HWSD=ovHWSD.lst[["BD"]]$BD, WoSIS=ovHWSD.lst[["BD"]]$WoSIS_BD), 
  data.frame(HWSD=ovHWSD.lst[["ph"]]$ph, WoSIS=ovHWSD.lst[["ph"]]$WoSIS_ph), 
  data.frame(HWSD=ovHWSD.lst[["OC"]]$OC, WoSIS=ovHWSD.lst[["OC"]]$WoSIS_OC), 
  data.frame(SoilGrids=ovSG.lst[["BD"]]$BD, WoSIS=ovSG.lst[["BD"]]$WoSIS_BD), 
  data.frame(SoilGrids=ovSG.lst[["ph"]]$ph, WoSIS=ovSG.lst[["ph"]]$WoSIS_ph), 
  data.frame(SoilGrids=ovSG.lst[["OC"]]$OC, WoSIS=ovSG.lst[["OC"]]$WoSIS_OC))
names(sel.plt) <- rep(c("BD fine earth t / cubic-m", "pH in water", "SOC in g/kg"), 2)
xlim.lst <- list(c(0.4,2.5), c(3.5,9.5), c(1,600), c(0.4,2.5), c(3.5,9.5), c(1,600))

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
  R.squared = signif(1-var(meas - pred, na.rm=TRUE)/var(meas, na.rm=TRUE), 2)
  main = paste0(varn, " (R-squared: ", R.squared, ")")
  if(!varn == "SOC in g/kg"){  
    plotList[[i]] <- hexbinplot(pred ~ meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=main, type="g", xlab=xlab, ylab=ylab, lwd=1, lcex=8, inner=.2, cex.labels=.8, xlim=xlim.lst[[i]], ylim=xlim.lst[[i]], asp=1, xbins=30, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) ## 
  } else {
    d.meas <- min(meas, na.rm=TRUE)
    pred <- pred+ifelse(d.meas==0, 1, d.meas)
    meas <- meas+ifelse(d.meas==0, 1, d.meas)
    plotList[[i]] <- hexbinplot(pred ~ meas, colramp=colorRampPalette(R_pal[["bpy_colors"]][1:18]), main=main, type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=30, xlim=xlim.lst[[i]], ylim=xlim.lst[[i]], xlab=xlab, ylab=ylab, panel=pfun, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) ## 
  }
}

do.call(grid.arrange, c(plotList[c(1,4)], ncol=2))
do.call(grid.arrange, c(plotList[c(2,5)], ncol=2))
do.call(grid.arrange, c(plotList[c(3,6)], ncol=2))

## Comparison with ISRIC WISE