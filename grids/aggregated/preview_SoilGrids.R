## Preview SoilGrids predictions (PNGs); generate aggregated maps in parallel, generate metadata (XML files)
## Tom.Hengl@isric.org

setwd("/data/aggregated")
load(".RData")
library(rgdal)
library(GSIF)
library(utils)
library(R.utils)
library(snowfall)
library(parallel)
library(raster)             
library(RSAGA)
library(grDevices)
library(plotKML)
library(maps)
library(RColorBrewer)
library(XML)

source("/data/models/mosaick_functions_ll.R")
plotKML.env(convert="convert", show.env=FALSE)
system("gdal-config --version")
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
## rename some legends:
soil.legends[4] = "BLDFIE"
soil.legends[5] = "CECSOL"

## Reproject some maps to Sinusoidal projection system ----
r.mod = raster("/data/MCD12Q1/LandCover_2001001_L1_500m.tif")
mod.grid = readRDS("/data/models/Sinusoidal_tiles_200km.rds")

latlon2sin(input.file="/data/GEOG/OCSTHA_M_30cm_250m_ll.tif", output.file="/data/GEOG/OCSTHA_M_30cm_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
latlon2sin(input.file="/data/GEOG/OCSTHA_M_100cm_250m_ll.tif", output.file="/data/GEOG/OCSTHA_M_100cm_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
latlon2sin(input.file="/data/GEOG/OCSTHA_M_200cm_250m_ll.tif", output.file="/data/GEOG/OCSTHA_M_200cm_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
latlon2sin(input.file="/data/GEOG/TAXNWRB_250m_ll.tif", output.file="/data/GEOG/TAXNWRB_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
latlon2sin(input.file="/data/GEOG/TAXOUSDA_250m_ll.tif", output.file="/data/GEOG/TAXOUSDA_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
latlon2sin(input.file="/data/GEOG/OCSTHA_M_30cm_250m_ll.tif", output.file="/data/GEOG/OCSTHA_M_30cm_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))
## GAUL 2014 data set (https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/countries/GAUL_250m.R)
latlon2sin(input.file="/data/countries/GAUL_ADMIN1_landmask_250m.tif", output.file="/data/GEOG/GAUL_ADMIN1_landmask_300m_sin.tif", mod.grid=mod.grid, pixsize=300, te=paste(as.vector(extent(r.mod))[c(1,3,2,4)], collapse=" "))

## Total soil organic carbon stock per GAUL ----
summary_OCS_tiles <- function(i, tileS.tbl, admin="/data/GEOG/GAUL_ADMIN1_landmask_300m_sin.tif", ocs="/data/GEOG/OCSTHA_M_100cm_300m_sin.tif"){
  m = readGDAL(fname=admin, offset=unlist(tileS.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
  m@data[,2] = readGDAL(fname=ocs, offset=unlist(tileS.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
  names(m) = c("Value","OCS")
  if(sum(!is.na(m$OCS))>0){
    ## Aggregate per class combination
    SOC_agg.admin <- plyr::ddply(m@data, .(Value), summarize, Total_OCS_mT=sum(OCS*300^2/10000, na.rm=TRUE)/1e6, Area_km2=sum(!is.na(OCS))*0.09)
    return(SOC_agg.admin)
  }
}

tileS.tbl = readRDS("/data/models/tileS.tbl.rds")
ov_ADMIN_sin = readRDS("/data/models/Sinusoidal_tiles_200km.rds")
tS.sel = as.character(ov_ADMIN_sin$ID)
hb.leg = read.csv("/data/countries/g2015_2014_1_legend.csv")
## run in parallel:
sfInit(parallel=TRUE, cpus=48)
sfExport("summary_OCS_tiles", "tileS.tbl", "tS.sel")
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(dplyr)
sfLibrary(utils)
out.lst <- sfClusterApplyLB(as.numeric(tS.sel), function(x){ summary_OCS_tiles(x, tileS.tbl=tileS.tbl) })
sfStop()
Tocs = rbind.fill(out.lst)
str(Tocs)
summary(is.na(Tocs$Value)) ## many missing values for country?
Tocs.sum = plyr::ddply(Tocs, .(Value), summarize, OCS_mT=sum(Total_OCS_mT, na.rm=TRUE), Total_Area_km2=sum(Area_km2, na.rm=TRUE))
Tocs.sum$ADM0_NAME = plyr::join(Tocs.sum, hb.leg)$ADM0_NAME
Tocs.sum$ADM1_NAME_filename = plyr::join(Tocs.sum, hb.leg, by="Value")$ADM1_NAME_filename
Tocs.sum$OCS_t_ha = round(Tocs.sum$OCS_mT*1e6/(Tocs.sum$Total_Area_km2*100))
Tocs.sum$OCS_mT = round(Tocs.sum$OCS_mT, 1)
write.csv(Tocs.sum, "Status_OCS_100cm_per_Country.csv")
## Total stock 0-100cm:
sum(Tocs.sum$OCS_mT[!is.na(Tocs.sum$Value)])

sfInit(parallel=TRUE, cpus=48)
sfExport("summary_OCS_tiles", "tileS.tbl", "tS.sel")
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(dplyr)
sfLibrary(utils)
out.lst <- sfClusterApplyLB(as.numeric(tS.sel), function(x){ summary_OCS_tiles(x, tileS.tbl=tileS.tbl, ocs="/data/GEOG/OCSTHA_M_30cm_300m_sin.tif") })
sfStop()
Tocs30 = rbind.fill(out.lst)
str(Tocs30)
Tocs30.sum = plyr::ddply(Tocs30, .(Value), summarize, OCS_mT=sum(Total_OCS_mT, na.rm=TRUE), Total_Area_km2=sum(Area_km2, na.rm=TRUE))
Tocs30.sum$ADM0_NAME = plyr::join(Tocs30.sum, hb.leg)$ADM0_NAME
Tocs30.sum$ADM1_NAME_filename = plyr::join(Tocs30.sum, hb.leg, by="Value")$ADM1_NAME_filename
Tocs30.sum$OCS_t_ha = round(Tocs30.sum$OCS_mT*1e6/(Tocs30.sum$Total_Area_km2*100))
Tocs30.sum$OCS_mT = round(Tocs30.sum$OCS_mT, 1)
write.csv(Tocs30.sum, "Status_OCS_30cm_per_Country.csv")
## Total stock 0-30cm:
sum(Tocs30.sum$OCS_mT[!is.na(Tocs30.sum$Value)])

## Compare with 1km data:
grid1km.sin = readGDAL("/data/aggregated/1km/OCSTHA_M_30cm_1km_sin.tif")
sum(grid1km.sin$band1*1000^2/10000, na.rm=TRUE)/1e6
## 1360 Pg C
rm(grid1km.sin); gc()

## List of property maps ----
tif.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
## 318
#tif.lst[grep("WRB",tif.lst)]

## Resample to 1 km (takes 20 mins) ----
sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "aggr_SG")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(tif.lst, aggr_SG)
sfStop()

sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "aggr_SG")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(tif.lst, aggr_SG, tr=1/20, tr.metric=5000, out.dir="/data/aggregated/5km/", ti="250m", tn="5km")
sfStop()
sfInit(parallel=TRUE, cpus=48)
sfExport("tif.lst", "aggr_SG")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(tif.lst, aggr_SG, tr=1/10, tr.metric=10000, out.dir="/data/aggregated/10km/", ti="250m", tn="10km")
sfStop()

## create PNGs (5km res)
out.lst <- list.files(path="/data/aggregated/5km", pattern=glob2rx("*_5km_ll.tif$"), full.names=TRUE, recursive=TRUE)

plot_SG <- function(i, res=150, zlim=c(0,40), width=7200, height=2987, ylim=c(-60, 85), replace=TRUE, country, soil.legends){
  out.file = paste0(normalizeFilename(basename(i)), ".png")
  if(replace==TRUE){
    r <- raster(i)
    isF = length(grep(pattern="WRB", basename(i)))>0 | length(grep(pattern="USDA", basename(i)))>0
    if(isF){
      col_scale = R_pal[["bpy_colors"]]
      breaks = c(seq(zlim[1], zlim[2], length=20), 100)
    } else {
      if(length(grep(pattern="ORC", basename(i)))>0){
        col_scale = SAGA_pal[[1]]
        breaks = c(0, 1, 2, 4, 6, 8, 11, 14, 18, 21, 25, 32, 40, 54, 70, 100, 140, 180, 240, 300, 600)
      }
      for(j in c("PHIHOX","PHIKCL","BLDFIE","CECSOL","SNDPPT","SLTPPT","CLYPPT")){
        if(length(grep(pattern=j, basename(i)))>0){
          col_scale = soil.legends[[j]]$COLOR
          breaks = c(soil.legends[[j]]$MIN[1], soil.legends[[j]]$MAX)
        }
      }
      if(length(grep(pattern="CRFVOL", basename(i)))>0){
        col_scale = SAGA_pal[[2]]
        breaks = c(0, 1, 2, 4, 6, 8, 10, 14, 19, 25, 32, 40, 54, 70, 88, 120, 190, 270, 430, 600, 1000)/10
      }
      if(length(grep(pattern="BDT", basename(i)))>0){
        col_scale = SAGA_pal[[1]]
        breaks = c(0, 1, 2, 4, 6, 8, 14, 20, 34, 46, 64, 85, 105, 140, 180, 220, 310, 450, 650, 950, 1800)*100
      }
      if(length(grep(pattern="BDR", basename(i)))>0){
        col_scale = rev(R_pal[["heat_colors"]])
        breaks = seq(0,260,length=21)
      }
      if(length(grep(pattern="BDRLOG", basename(i)))>0){
        col_scale = R_pal[["heat_colors"]]
        breaks = seq(0,100,length=21)
      }
    }
    ri = cut(r, breaks=breaks, include.lowest=TRUE, right=FALSE)
    png(file = out.file, res = res, width = width, height = height, type="cairo")
    #dev.new(width = 20, height = 8.3)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    image(ri, col=col_scale, ylim=ylim)
    lines(country, col="black")
    legend("bottomleft", legend=signif(rev(breaks[-1]), 3), fill=rev(col_scale), horiz=FALSE, pt.cex=2)
    dev.off()
  }
}

sfInit(parallel=TRUE, cpus=24)
sfExport("out.lst", "plot_SG", "country", "soil.legends")
sfLibrary(rgdal)
sfLibrary(raster)
sfLibrary(plotKML)
sfLibrary(grDevices)
out <- sfClusterApplyLB(out.lst, function(i){ try( plot_SG(i, country=country, soil.legends=soil.legends) ) } )
sfStop()

for(k in c("TAXNWRB","TAXOUSDA")){
  ## cluster:
  cl <- makeCluster(getOption("cl.cores", 48))
  ## Most probable class:
  cl.lst <- list.files("/data/GEOG/aggr", pattern=k, full.names=TRUE, recursive=TRUE)
  cl.legend <- read.csv(paste0("/data/models/", k, "/", k, "_legend.csv"))
  cl.legend$COLOR <- rgb(red=cl.legend$R/255, green=cl.legend$G/255, blue=cl.legend$B/255)
  cl.levs <- sapply(cl.lst, function(x){strsplit(x, "_")[[1]][2]})
  cl.tbl <- join(data.frame(Filename=cl.lst, Shortened_name=cl.levs, int=1:length(cl.levs)), cl.legend, typ="left")
  ## load all predictions in memory:
  probs <- stack(cl.lst)
  probs <- as(probs, "SpatialGridDataFrame")
  ## select only pixels with values:
  sel.pix <- !is.na(probs@data[,1])
  ## rank classes (takes 2+ mins):
  ranks <- data.frame(parallel::parApply(cl, -probs@data[sel.pix,1:length(cl.lst)], MARGIN=1, rank) )
  gc()
  ## most probable classes:
  probs@data[sel.pix,"cl1"] <- unlist(parallel::mclapply(ranks, function(x){which(x==1)}, mc.cores=48))
  writeGDAL(probs["cl1"], "TAXNWRB_1st_5km_ll.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(cl.tbl$Group)), colorTable=list(col.tbl$COLOR))
  writeGDAL(probs["cl2"], "TAXNWRB_2nd_5km_ll.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(cl.tbl$Group)), colorTable=list(col.tbl$COLOR))
  writeGDAL(probs["cl3"], "TAXNWRB_3rd_5km_ll.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(paste(cl.tbl$Group)), colorTable=list(col.tbl$COLOR))
  stopCluster(cl)
}

## Rename some files:
tif.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*.tif$"), full.names = TRUE, recursive = FALSE)
file.rename(from=tif.lst[grep("BLD", tif.lst)], to=gsub("BLD","BLDFIE",tif.lst[grep("BLD", tif.lst)]))
file.rename(from=tif.lst[grep("CECSUM", tif.lst)], to=gsub("CECSUM","CECSOL",tif.lst[grep("CECSUM", tif.lst)]))
## final list:
tif250m.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*_250m_ll.tif$"), full.names = TRUE, recursive = FALSE)
write.csv(data.frame(FileName=basename(tif250m.lst)), "tif250m.lst.csv")
## 252 files

## Legends for soil types ----
TAXOUSDA.leg <- read.csv("/data/models/TAXOUSDA/TAXOUSDA_legend.csv")
TAXOUSDA.leg = TAXOUSDA.leg[!is.na(TAXOUSDA.leg$R),]
cat(paste0(sapply(1:nrow(TAXOUSDA.leg), function(i){paste0(TAXOUSDA.leg$Number[i], ': \"', TAXOUSDA.leg$Group[i], '\"')}), collapse=", "), file = "TAXOUSDA_legend.txt")
TAXNWRB.leg <- read.csv("/data/models/TAXNWRB/TAXNWRB_legend.csv")
TAXNWRB.leg = TAXNWRB.leg[!is.na(TAXNWRB.leg$R),]
cat(paste0(sapply(1:nrow(TAXNWRB.leg), function(i){paste0(TAXNWRB.leg$Number[i], ': \"', TAXNWRB.leg$Group[i], '\"')}), collapse=", "), file = "TAXNWRB_legend.txt")

## Generate Inspire metadata files for each zipped GEOTIFF ----
mdSG.s <- read.csv("META_GEOTIFF_Stacked.csv", stringsAsFactors = FALSE)
doc = xmlInternalTreeParse("Metadata_template_7d_250m_ll.xml")
xm7d = xmlRoot(doc)
xm.ls = unlist(xmlToList("Metadata_template_7d_250m_ll.xml"))
xm.df = data.frame(Fields=attr(xm.ls, "names"), Values=xm.ls)
write.csv(xm.df, "Metadata_template_7d_250m_ll.csv")

## Function to generate INSPIRE compatible XMLs for Geonetwork (https://goo.gl/V84S5h):
SoilGrids250m_XML = function(i, template.xml="Metadata_template_7d_250m_ll.xml", mdSG.s, out.dir="/data/GEOG/"){
  out.xml = paste0(out.dir, mdSG.s$GSIF_id[i], mdSG.s$XML_ext[i], ".xml")
  if(!file.exists(out.xml)){
    doc = xmlInternalTreeParse(template.xml)
    xm7d = xmlRoot(doc)
    xmlValue(xm7d[["fileIdentifier"]][[1]]) = mdSG.s$gmd.fileIdentifier[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["title"]][[1]]) = mdSG.s$gmd.title[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["identifier"]][["RS_Identifier"]][["code"]][[1]]) = mdSG.s$GSIF_id[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["edition"]][[1]]) = stringr::str_sub(mdSG.s$gco.Date[i], 1, 7)
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["date"]][["CI_Date"]][["date"]][[1]]) = mdSG.s$gco.Date[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["abstract"]][[1]]) = mdSG.s$gmd.abstract[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["graphicOverview"]][["MD_BrowseGraphic"]][["fileName"]][[1]]) = mdSG.s$graphicOverview[i]
    xmlValue(xm7d[["identificationInfo"]][["MD_DataIdentification"]][["descriptiveKeywords"]][["MD_Keywords"]][[2]][[1]]) = mdSG.s$keyword2[i]
    ## Files available:
    removeNodes(xm7d[["distributionInfo"]][["MD_Distribution"]][["transferOptions"]][["MD_DigitalTransferOptions"]])
    x = list.files(path=out.dir, pattern=glob2rx(mdSG.s$Filenames[i]))
    if(length(x)==1){ xt = paste0("Download GeoTIFF map for soil stratum") } 
    if(length(x)==3){ xt = paste0("Download GeoTIFF map for 0-", c(30,100,200), " cm depth interval") }
    if(length(x)==6){ xt = paste0("Download GeoTIFF map for ", c(0,5,15,30,60,100),"-", c(5,15,30,60,100,200), " cm depth interval") }
    if(length(x)==7){ xt = paste0("Download GeoTIFF map for ", c(0,5,15,30,60,100,200), " cm depth") }
    if(length(x)>7){ xt = paste0("Download GeoTIFF map for probability") } 
    CI_o <- sprintf('<gmd:onLine><gmd:CI_OnlineResource><gmd:linkage><gmd:URL>ftp://ftp.soilgrids.org/data/recent/%s</gmd:URL></gmd:linkage><gmd:protocol><gco:CharacterString>WWW:DOWNLOAD-1.0-ftp--download</gco:CharacterString></gmd:protocol><gmd:name><gco:CharacterString>%s</gco:CharacterString></gmd:name><gmd:function><gmd:CI_OnLineFunctionCode codeList="http://standards.iso.org/ittf/PubliclyAvailableStandards/ISO_19139_Schemas/resources/codelist/ML_gmxCodelists.xml#CI_OnLineFunctionCode" codeListValue="download"/></gmd:function></gmd:CI_OnlineResource></gmd:onLine>', x, xt)
    ## WMS available:
    CI_w <- sprintf('<gmd:onLine><gmd:CI_OnlineResource><gmd:linkage><gmd:URL>http://data.isric.org/geoserver/sg250m/wms?</gmd:URL></gmd:linkage><gmd:protocol><gco:CharacterString>OGC:WMS</gco:CharacterString></gmd:protocol><gmd:name><gco:CharacterString>%s</gco:CharacterString></gmd:name><gmd:description><gco:CharacterString>%s</gco:CharacterString></gmd:description><gmd:function><gmd:CI_OnLineFunctionCode codeList="http://standards.iso.org/ittf/PubliclyAvailableStandards/ISO_19139_Schemas/resources/codelist/ML_gmxCodelists.xml#CI_OnLineFunctionCode" codeListValue="information"/></gmd:function></gmd:CI_OnlineResource></gmd:onLine>', gsub(".tif", "", gsub("_ll.tif", "", x)), gsub("Download GeoTIFF map for ", "", xt))
    pl1 = newXMLNode("gmd:MD_DigitalTransferOptions", parent=xm7d[["distributionInfo"]][["MD_Distribution"]][["transferOptions"]])
    parseXMLAndAdd(c(CI_o, CI_w), parent=pl1, nsDefs  = c(gmd="http://www.isotc211.org/2005/gmd", gco="http://www.isotc211.org/2005/gco")) 
    ## Write to a file:
    saveXML(xm7d, out.xml) 
  }
}
x = sapply(1:nrow(mdSG.s), SoilGrids250m_XML, mdSG.s=mdSG.s)
file.copy("/data/models/TAXOUSDA/TAXOUSDA_legend.csv", "/data/GEOG/TAXOUSDA_250m_ll.tif.csv")
file.copy("/data/models/TAXNWRB/TAXNWRB_legend.csv", "/data/GEOG/TAXNWRB_250m_ll.tif.csv")

## Old code:
mdSG <- read.csv("META_GEOTIFF_1B.csv")
data(landmask)
gridded(landmask) <- ~x+y
proj4string(landmask) <- "+proj=longlat +datum=WGS84"
landmask <- as(landmask["mask"], "SpatialPixelsDataFrame")
landmask <- landmask[landmask@coords[,2]> -62 & landmask@coords[,2]< 87.4,]

for(i in 1:301){
 out.xml.file = paste0("X:/SoilGrids250m/GEOG/", mdSG$FileName[i], ".xml") # "/data/GEOG/"
  if(!file.exists(out.xml.file)){
   xmd <- spMetadata(landmask, out.xml.file=out.xml.file,
    md.type="INSPIRE",
    CI_Citation_title = paste(mdSG$SERIES_NAME[i], ":", mdSG$ATTRIBUTE_LABEL[i], ":", mdSG$ATTRIBUTE_TITLE[i]),
    ## SoilGrids250m : PHIHOX_M_sl1 : Soil pH x 10 in H2O  at depth 0.00 m
    #CI_Online_resource_URL = paste(mdSG$DOWNLOAD_FTP_URL[i], mdSG$FileName[i], ".gz", sep=""),
    CI_Electronic_mail_address = mdSG$CITATION_ADDRESS[i],
    CI_Organisation_name = mdSG$CITATION_ORIGINATOR[i],
    Date_stamp = as.Date(mdSG$PUBLICATION_DATE[i], format="%Y-%m-%d"),
    CI_Citation_date = as.Date(mdSG$PUBLICATION_DATE[i], format="%Y-%m-%d"),
    MD_Thesaurus_date = as.Date(mdSG$PUBLICATION_DATE[i], format="%Y-%m-%d"),
    CI_Unique_name = mdSG$ATTRIBUTE_LABEL[i],
    MD_Abstract = mdSG$ATTRIBUTE_TITLE[i],
    MD_Organisation_name = "ISRIC - World Soil Information",
    MD_Electronic_mail_address = "tom.hengl@isric.org",
    MD_Keyword = c("soil", paste(mdSG$KEYWORD1[i]), paste(mdSG$KEYWORD2[i])),
    MD_Use_limitations = "Open Data Commons Open Database License (ODbL). This means that: You must attribute any public use of the database, or works produced from the database, in the manner specified in the ODbL. For any use or redistribution of the database, or works produced from it, you must make clear to others the license of the database and keep intact any notices on the original database. Share-Alike: If you publicly use any adapted version of this database, or works produced from an adapted database, you must also offer that adapted database under the ODbL. Keep open: If you redistribute SoilGrids geotifs, or an adapted version of it, then you may use technological measures that restrict the work (such as DRM) as long as you also redistribute a version without such measures.",
    MD_Other_restrictions = "http://www.isric.org/content/disclaimer-soilgrids",
    MD_Equivalent_scale = "200000",
    MD_Resolution = 250,
    Time_period_begin = "1950-01-01", ## as.Date("1950-01-01", format="%Y-%m-%d"),
    Time_period_end = "2015-12-31",  ## as.Date("2005-12-31", format="%Y-%m-%d")
    Extent_West_Longitude = -180,
    Extent_East_Longitude = 180,
    Extent_South_Latitude = -62,
    Extent_North_Latitude = 87.4,
    DQ_Lineage_statement = paste("Values (", ifelse(mdSG$CONFIDENCE_INTERVAL[i]=="M", "mean value", "upper (U) or lower (L) confidence limit"), ") predicted using statistical modelling. Measurement unit: ", mdSG$ATTRIBUTE_UNITS_OF_MEASURE[i], ". The automated soil mapping system 'SoilGrids' is explained in detail at: ", mdSG$PROJECT_URL[i], sep="")
  )
 }
}

## categories separately because some classes can get 'kicked-out' during modelling:
m_TAXOUSDA <- readRDS("/data/models/TAXOUSDA/mnetX_TAXOUSDA.rds")
USDAf <- data.frame(Group=m_TAXOUSDA$finalModel$lev, levs=1:length(m_TAXOUSDA$finalModel$lev))
USDA.col <- read.csv("/data/models/TAXOUSDA/TAXOUSDA_legend.csv")
USDA.col <- USDA.col[!is.na(USDA.col$R),]
USDA.col$COLOR <- rgb(USDA.col$R/255, USDA.col$G/255, USDA.col$B/255) 
rm(m_TAXOUSDA)
USDA.col <- join(USDAf, USDA.col[,c("Group","COLOR")], type="left")

m_TAXNWRB <- readRDS("/data/models/TAXNWRB/mnetX_TAXNWRB.rds")
WRBf <- data.frame(Group=m_TAXNWRB$finalModel$lev, levs=1:length(m_TAXNWRB$finalModel$lev))
WRB.col <- read.csv("/data/models/TAXNWRB/TAXNWRB_legend.csv")
WRB.col <- WRB.col[!is.na(WRB.col$R),]
WRB.col$COLOR <- rgb(WRB.col$R/255, WRB.col$G/255, WRB.col$B/255) 
rm(m_TAXNWRB)
WRB.col <- join(WRBf, WRB.col[,c("Group","COLOR")], type="left")

USDA.col$MIN <- USDA.col$levs; USDA.col$MAX <- USDA.col$levs+1
WRB.col$MIN <- WRB.col$levs; WRB.col$MAX <- WRB.col$levs+1
legend.sdl <- list(WRB.col, USDA.col)
names(legend.sdl) = c("TAXNWRB","TAXOUSDA")

## SLD files for classes:
for(j in 1:length(legend.sdl)){
  sld.file <- file(set.file.extension(names(legend.sdl)[j], ".sld"), "w", blocking=FALSE)
  l1 = newXMLNode("StyledLayerDescriptor", attrs=c("xsi:schemaLocation" = "http://www.opengis.net/sld StyledLayerDescriptor.xsd", version="1.0.0"), namespaceDefinitions=c("http://www.opengis.net/sld", "xsi" = "http://www.w3.org/2001/XMLSchema-instance", "ogc" = "http://www.opengis.net/ogc", "gml" = "http://www.opengis.net/gml"))
  l2 <- newXMLNode("NamedLayer", parent = l1)
  l3 <- newXMLNode("Name", "SoilGrids250m", parent = l2)
  l3b <- newXMLNode("UserStyle", parent = l2)
  l4 <- newXMLNode("Title", paste("SoilGrids250m", names(legend.sdl)[j], sep="_"), parent = l3b)
  l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
  l5 <- newXMLNode("Rule", parent = l4b)
  l6 <- newXMLNode("RasterSymbolizer", parent = l5)
  Group = paste(legend.sdl[[j]]$Group)
  levs = legend.sdl[[j]]$MAX   
  l7 <- newXMLNode("ColorMap", attrs=c(type="intervals"), parent = l6)
  txt <- sprintf('<ColorMapEntry color="%s" quantity="%.0f" label="%s" opacity="%.1f"/>', c("#FFFFFF", strtrim(legend.sdl[[j]]$COLOR, 7)), c(1, levs), c("NODATA", Group), c(0.0, rep(.7, length(legend.sdl[[j]]$COLOR))))  
  parseXMLAndAdd(txt, l7)
  saveXML(l1, sld.file)
  close(sld.file)
}

var.lst <- list("BDTICM", "BDRICM", "BDRLOG")
title.lst <- list("BDTICM_M_250m","BDRICM_M_250m", "BDRLOG_M_250m")  
name.lst <- list(paste(mdSG$VARIABLE_NAME[which(mdSG$FileName=="BDTICM_M_250m_ll.tif")]), paste(mdSG$VARIABLE_NAME[which(mdSG$FileName=="BDRICM_M_250m_ll.tif")]), paste(mdSG$VARIABLE_NAME[which(mdSG$FileName=="BDRLOG_M_250m_ll.tif")]))  
nodata.lst <- list(-99999, 255, 255, -32768)
levs.lst <- list(c(0, 30, 90, 150, 240, 360, 580, 880, 1050, 1200, 1350, 1650, 1900, 2100, 24000, 3200, 5400, 9100, 15000, 54121), c(0, 30, 60, 90, 120, 160, 190, 200), c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
col.lst <- list(colorRampPalette(brewer.pal(name = "YlGnBu", n=9))(length(levs.lst[[1]])), colorRampPalette(rev(brewer.pal(name = "Oranges", n=9)))(length(levs.lst[[2]])), colorRampPalette(brewer.pal(name = "YlOrRd", n=9))(length(levs.lst[[3]])))
  
## ## SLD files for properties:
for(j in 1:length(title.lst)){
  sld.file <- file(paste0(var.lst[[j]], ".sld"), "w", blocking=FALSE)
  l1 = newXMLNode("StyledLayerDescriptor", attrs=c("xsi:schemaLocation" = "http://www.opengis.net/sld StyledLayerDescriptor.xsd", version="1.0.0"), namespaceDefinitions=c("http://www.opengis.net/sld", "xsi" = "http://www.w3.org/2001/XMLSchema-instance", "ogc" = "http://www.opengis.net/ogc", "gml" = "http://www.opengis.net/gml"))
  l2 <- newXMLNode("NamedLayer", parent = l1)
  l3 <- newXMLNode("Name", "SoilGrids250m", parent = l2)
  l3b <- newXMLNode("UserStyle", parent = l2)
  l4 <- newXMLNode("Title", title.lst[[j]], parent = l3b)
  l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
  l5 <- newXMLNode("Rule", parent = l4b)
  l6 <- newXMLNode("RasterSymbolizer", parent = l5)
  l7 <- newXMLNode("Geometry", parent = l6)
  l8 <- newXMLNode("PropertyName", name.lst[[j]], parent = l7)
  l9 <- newXMLNode("Opacity", 1, parent = l7)
  l10 <- newXMLNode("ColorMap", attrs=c(type="intervals"), parent = l6)
  txt <- sprintf('<ColorMapEntry color="%s" label="%s" opacity="%.1f" quantity="%.1f"/>', c("#FFFFFF", strtrim(col.lst[[j]], 7)), c("NODATA", levs.lst[[j]]), c(0, rep(0.7, length(levs.lst[[j]]))), c(nodata.lst[[j]], levs.lst[[j]]))  
  parseXMLAndAdd(txt, l10)
  saveXML(l1, sld.file)
  close(sld.file)
}

## Add soil legends to GSIF:
data("soil.legends")
names(soil.legends) = c("ORCDRC","PHIHOX","PHIKCL","BLDFIE","CECSOL","SNDPPT","SLTPPT", "CLYPPT","CRFVOL","TAXOUSDA","TAXGWRB")
legUSDA <- read.csv("../models/TAXOUSDA/TAXOUSDA_legend.csv")
legUSDA <- legUSDA[!is.na(legUSDA$B),]
legWRB <- read.csv("../models/TAXNWRB/TAXNWRB_legend.csv")
legWRB <- legWRB[!is.na(legWRB$B),]
soil.legends$TAXOUSDA <- legUSDA[,c("Number","Group","Generic")]
soil.legends$TAXOUSDA$COLOR <- rgb(red=legUSDA$R/255, green=legUSDA$G/255, blue=legUSDA$B/255)
soil.legends$TAXNWRB <- legWRB[,c("Number","Group","Shortened_name")]
soil.legends$TAXNWRB$Generic <- legWRB$WRB_group
soil.legends$TAXNWRB$COLOR <- rgb(red=legWRB$R/255, green=legWRB$G/255, blue=legWRB$B/255)
save(soil.legends, file="D:/Rdev/GSIF/pkg/data/soil.legends.rda", compress="xz")

## Soil correlation tables:
csv.lst <- c("../../profiles/CanSIS/CAN_classes.csv", "../../profiles/China/cleanup_SU_SYM90.csv", "../profs/TAXOUSDA/TAXOUSDA_GreatGroups.csv", "../profs/TAXNWRB/WRB_versions.csv", "../profs/TAXNWRB/cleanup_SU_SYM74.csv", "../profs/TAXOUSDA/USDA_Great_Group_2_FAO.csv", "../profs/TAXNWRB/Full_data_FAO74_US_CPSS.csv") ## "../../profiles/Alterra/cleanup_Alterra.csv",, "../../profiles/Radambrasil/cleanup_RadamBrasil.csv",
soil.classes <- lapply(csv.lst, read.csv)
names(soil.classes) = c("Canadian", "FAO1990.WRB", "USDA_GreatGroups", "WRB_versions", "FAO1974.WRB", "USDA.WRB", "Soils_World")
soil.classes$Soils_World$SoilUnitDCPCS <- iconv(soil.classes$Soils_World$SoilUnitDCPCS, to="ASCII", sub="byte")
soil.classes$Soils_World$SoilGroupDCPCS <- iconv(soil.classes$Soils_World$SoilGroupDCPCS, to="ASCII", sub="byte")
soil.classes$Soils_World$SoilGroupCPCS <- iconv(soil.classes$Soils_World$SoilGroupCPCS, to="ASCII", sub="byte")
soil.classes$Soils_World$CPSS.from.report..tentative. <- iconv(soil.classes$Soils_World$CPSS.from.report..tentative., to="ASCII", sub="byte")
save(soil.classes, file="D:/Rdev/GSIF/pkg/data/soil.classes.rda", compress="xz")
DGC <- which(soil.classes$Canadian$CSSC_Great_Groups=="Dark Gray Chernozem")
soil.classes$Canadian[DGC,]
soil.classes$Soils_World[2,]

