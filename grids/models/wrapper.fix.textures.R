## Function to correct values for textures fractions (force sum to 100%)
## Tom.Hengl@isric.org

wrapper.fix.textures <- function(i, n.lst=c("SNDPPT","SLTPPT","CLYPPT"), SND.lst){
  tex.lst <- sapply(n.lst, function(x){gsub(pattern="SNDPPT", replacement=x, SND.lst[i])})
  if(file.exists(tex.lst[1])){
    x = readGDAL(tex.lst[1])
    x@data[,2] <- readGDAL(tex.lst[2])$band1
    x@data[,3] <- readGDAL(tex.lst[3])$band1
    names(x) <- n.lst
    gc()
    sums <- rowSums(x@data)
    xr = range(sums, na.rm=TRUE)
    gc()
    if(xr[1]<99|xr[2]>101){
      x$SLTPPT <- round(x$SLTPPT / sums * 100, 0)
      x$SNDPPT <- round(x$SNDPPT / sums * 100, 0)
      x$CLYPPT <- round(x$CLYPPT / sums * 100, 0)
      unlink(tex.lst[1]); unlink(tex.lst[2]); unlink(tex.lst[1])
      writeGDAL(x[n.lst[1]], tex.lst[1], "GTiFF", mvFlag=255, type="Byte", options="COMPRESS=DEFLATE")
      writeGDAL(x[n.lst[2]], tex.lst[2], "GTiFF", mvFlag=255, type="Byte", options="COMPRESS=DEFLATE")
      writeGDAL(x[n.lst[3]], tex.lst[3], "GTiFF", mvFlag=255, type="Byte", options="COMPRESS=DEFLATE")
      gc()
    }
  }
}

SND.lst <- list.files(path="/data/predicted", pattern="SNDPPT", recursive = TRUE, full.names = TRUE)
str(SND.lst)

sfInit(parallel=TRUE, cpus=48)
sfLibrary(rgdal)
sfLibrary(sp)
sfExport("SND.lst", "wrapper.fix.textures")
x <- sfClusterApplyLB(1:length(SND.lst), wrapper.fix.textures, SND.lst=SND.lst)
sfStop()
gc()
rm(x)