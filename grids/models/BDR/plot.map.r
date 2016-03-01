rm(list = ls(all = TRUE))
library(sp)
a.dir <- "/data/shang009/big/soildepth2"# dir of the project
setwd(a.dir)

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
p.at <- sapply(sapply(seq(0,5.5,0.5), function(x)10^x), log1p)
plotList[[1]] <- spplot(g, zcol = "band1", scales=list(draw = TRUE), at = p.at,main = "Absolute depth to bedrock (cm)", colorkey = list(space = "right", height = 0.4), formula = as.formula(log1p(band1)~x+y),col.regions = SAGA_pal[[1]])
plotList[[1]]$legend$right$args$key$labels$at <- sapply(sapply(seq(0,5,1), function(x)10^x), log1p)
plotList[[1]]$legend$right$args$key$labels$labels <- c(0,10,100,1000,10000,"100000")


g <- readGDAL("BDRICM_M_10km_ll.tif")
sum(g$band1>200,na.rm=T)/sum(!is.na(g$band1))
p.at <- seq(0,250,25)
plotList[[2]] <- spplot(g, zcol = "band1", colorkey = list(space = "right", height = 0.4), scales=list(draw = TRUE), at = p.at,main = "Censored depth to bedrock (cm)",col.regions = SAGA_pal[[1]])

g <- readGDAL("BDRLOG_M_10km_ll.tif")
p.at <- seq(0,100,10)
plotList[[3]] <- spplot(g, zcol = "band1", colorkey = list(space = "right", height = 0.4),scales=list(draw = TRUE), at = p.at, main = "Probability of the occurrence of R horizon within 250 cm (%)",col.regions = SAGA_pal[["SG_COLORS_YELLOW_RED"]])

bitmap(paste0("./pics/Fig_BDTICM.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[1]]
#do.call(grid.arrange, c(plotList, nrow = 3))
dev.off()

bitmap(paste0("./pics/Fig_BDRICM.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[2]]
#do.call(grid.arrange, c(plotList, nrow = 3))
dev.off()
bitmap(paste0("./pics/Fig_BDRLOG.tif"), width = 7.48, height = 4, units = "in", res =300, type = "tiff32nc")
plotList[[3]]
#do.call(grid.arrange, c(plotList, nrow = 3))
dev.off()
