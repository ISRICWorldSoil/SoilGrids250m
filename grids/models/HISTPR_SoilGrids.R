## Distribution of organic soils (probability of histosols) based on SoilGrids250m
## Tom.Hengl@isric.org

library(rgdal)
fao.lst <- c("Sapric.Histosols", "Hemic.Histosols", "Fibric.Histosols", "Cryic.Histosols", "Histic.Albeluvisols")
usda.lst <- c("Saprists", "Hemists", "Folists", "Fibrists")

histosol.prob <- function(i, in.path, fao.lst, usda.lst){
  out.p <- paste0(in.path, "/", i, "/HISTPR_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- c(paste0(in.path, "/", i, "/TAXNWRB_", fao.lst, "_", i, ".tif"), paste0(in.path, "/", i, "/TAXOUSDA_", usda.lst, "_", i, ".tif"))
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    names(s) <- c(fao.lst, usda.lst)
    gc()
    s$HISTPR <- (s@data[,"Sapric.Histosols"] + s@data[,"Saprists"])/2 + (s@data[,"Hemic.Histosols"] + s@data[,"Hemists"])/2 + (s@data[,"Fibric.Histosols"] + s@data[,"Fibrists"])/2 + s@data[,"Cryic.Histosols"] + s@data[,"Histic.Albeluvisols"] + s@data[,"Folists"]
    writeGDAL(s["HISTPR"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

in.path="G:/SoilGrids250m/predicted"
i = "NA_060_036"