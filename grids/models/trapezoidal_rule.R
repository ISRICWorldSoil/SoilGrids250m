## Trapezoidal rule (vertical aggregation function)
## https://www.wikiwand.com/en/Trapezoidal_rule
## Tom.Hengl@isric.org

#d = c(0,5,15,30,60,100,200)
#y = c(4.5,5.0,5.3,5.0,5.5,5.8,6.2)
#y1 = sum(c((d[2]-d[1])*(y[2]+y[1])/2, (d[3]-d[2])*(y[3]+y[2])/2, (d[4]-d[3])*(y[4]+y[3])/2))/d[4]
## same thing:
#weighted.mean((y[1:3]+y[2:4])/2, diff(d[1:4]))

## Aggregate values using variable depths:
agg_layers <- function(i, varn, d = c(0,10,35), in.path="/data/predicted", ot="Int16", dstnodata=-32768){
  in.lst <- paste0(path=in.path, "/", i, "/", varn, "_M_sl", 1:length(d), "_", i, ".tif")
  out.tif <- paste0(path=in.path, "/", i, "/", varn, "_M_agg", d[length(d)]-d[1], "cm_", i, ".tif")
  if(!file.exists(out.tif)){
    s <- stack(in.lst)
    s <- as(s, "SpatialGridDataFrame")
    ## Trapezoidal rule
    for(j in 1:(length(d)-1)){
      s@data[,paste0("sd",j)] <- rowMeans(s@data[,j:(j+1)])
    }
    s$sum_sd <- rowSums(as.matrix(s@data[,paste0("sd",1:(length(d)-1))]) %*% diag(diff(d)))/sum(diff(d))
    writeGDAL(s["sum_sd"], out.tif, type=ot, mvFlag=dstnodata, options="COMPRESS=DEFLATE")
  }
}
