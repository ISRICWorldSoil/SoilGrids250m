rm(list = ls(all = TRUE))
library(sp)
a.dir <- "/data/shang009/big/soildepth2"# dir of the project
setwd(a.dir)

country <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country$names, ":"), function(x) x[1])
country = as(map2SpatialPolygons(country, IDs=IDs), "SpatialLines")
b_poly <- as(extent(c(-180,180,-65,75)), "SpatialPolygons")
country = gIntersection(country, b_poly, byid = T)
proj4string(country) = "+proj=longlat +datum=WGS84"

tvar <- c("BDTICM", "BDRICM", "BDRLOG")
for(i in 2:3)
{
  src_d <- paste0(a.dir, paste0("/", tvar[i], "_M_1km_ll.tif"))
  out_d <- paste0(a.dir, paste0("/", tvar[i], "_M_10km_ll.tif"))
  system(paste("gdalwarp  -tr 0.08333333 0.08333333 -r average -overwrite",
               src_d, out_d))   
}
#system("gzip -k BD*10km*")

plotList <- NULL
g <- readGDAL("BDTICM_M_10km_ll.tif")
mean(g$band1,na.rm=T)
p.at <- sapply(sapply(seq(2,4.4,0.2), function(x)10^x), log1p)
plotList[[1]] <- spplot(g, zcol = "band1", scales=list(draw = FALSE), at = p.at,main = "Absolute depth to bedrock (cm)", colorkey = list(space = "right", height = 0.3), formula = as.formula(log1p(band1)~x+y),col.regions = SAGA_pal[[1]],sp.layout=list(country, col = 'black', lwd = 0.5))
plotList[[1]]$legend$right$args$key$labels$at <- c(sapply(sapply(seq(2,4.4,0.8), function(x)10^x), log1p))
plotList[[1]]$legend$right$args$key$labels$labels <- signif(sapply(seq(2,4.4,0.8), function(x)10^x),2)


g <- readGDAL("BDRICM_M_10km_ll.tif")
sum(g$band1>200,na.rm=T)/sum(!is.na(g$band1))
p.at <- seq(0,250,25)
plotList[[2]] <- spplot(g, zcol = "band1", colorkey = list(space = "right", height = 0.3), scales=list(draw = FALSE), at = p.at,main = "Censored depth to bedrock (cm)",col.regions = SAGA_pal[[1]],sp.layout=list(country, col = 'black', lwd = 0.5))

g <- readGDAL("BDRLOG_M_10km_ll.tif")
p.at <- seq(0,100,10)
plotList[[3]] <- spplot(g, zcol = "band1", colorkey = list(space = "right", height = 0.3),scales=list(draw = FALSE), at = p.at, main = "Probability of the occurrence of R horizon within 250 cm (%)",col.regions = SAGA_pal[["SG_COLORS_YELLOW_RED"]],sp.layout=list(country, col = 'black', lwd = 0.5))

bitmap(paste0("./pics/Fig_BDTICM.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[1]]
dev.off()

bitmap(paste0("./pics/Fig_BDRICM.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[2]]

dev.off()
bitmap(paste0("./pics/Fig_BDRLOG.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[3]]
dev.off()
