## predict numeric variables using a model "gm" and write GeoTifs out - SoilGrids250m
## by: Tom.Hengl@isric.org

wrapper.predict_n <- function(i, gm1, sd=c(2.5, 7.5, 22.5, 45, 80, 150), varn, in.path, out.path, z.min=0, z.max){ #gm2,
  m <- readRDS(paste0(in.path, "/", i, "/", i, ".rds"))
  for(j in 1:length(sd)){
    out.c <- paste0(out.path, "/", i, "/", varn, "_sd", j, "_", i, ".tif")
    if(!file.exists(out.c)){
      ## predict:
      m$DEPTH <- sd[j]
      m$v <- predict(gm1, m, na.action = na.pass)
      if(!missing(z.min)){
        m$v <- ifelse(m$v < z.min, z.min, m$v)
      }
      if(!missing(z.max)){
        m$v <- ifelse(m$v > z.max, z.max, m$v)
      }
      writeGDAL(m["v"], out, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      gc()
    }
  }
}