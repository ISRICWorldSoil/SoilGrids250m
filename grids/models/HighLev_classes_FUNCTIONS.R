##  function to derive dominant class in a reference group level: 
aggregateRGclass <- function(p, list){
 
 ### most probable class: 
 dfp <- data.frame(p[,1])
 
 for(z in 1:length(list)){
  names(p[grepl(list[z], names(p), perl=TRUE)])
  agg <- apply(X=p[names(p[grepl(list[z], names(p), perl=TRUE)])], MARGIN=1, FUN=sum)
  summary(agg)
  dfp <- cbind(dfp, agg)
 }
 
 dfp<- dfp[, 2:ncol(dfp)]
 colnames(dfp) <- list
 
 #### identify more probable class
 dfp$agg <- apply(X=dfp, MARGIN=1, FUN=
                    function(x){
                     pred <- names(x[which.max(x)])
                     return(pred)
                    })
 
 dfp$agg <- as.character(dfp$agg)
 tail(dfp)
 return(data.frame(agg=dfp$agg))} 

##  function to derive dominant class in a reference group level in the maps per tile; it is using the previous function and 2 external LUT tables (LUT_WRB, and LUT_USDA): 
aggregateRG_Maps <- function(i, in.path, out.path, varn, extent){ 
 
#  out.c <- paste0(out.path, "/tiles/",  varn, "_", i, "_RefG.tif") 
 out.c <- paste0(out.path, "/tiles/", i, "/", varn, "_", i, "_RefG.tif") 
 cl.lst <- list.files(path=(paste("./tiles/", i, sep="")), pattern=glob2rx(paste0(varn, "_*.tif$")), full.names=TRUE)
 library(stringr)
 mask <- grep(paste0(varn, "_",  substr(i, start = 0 , stop = 2)) , cl.lst)
 m <- raster::stack(cl.lst[-mask]) 
 names(m) <- sapply(names(m), function(x){strsplit(x, "_")[[1]][2]}) 
 m <- as(m, "SpatialGridDataFrame") 
 
 if(varn=="TAXNWRB"){
  lut<- read.csv(file = "./LUT_WRB.csv")
  rg <- c("Acrisols", "Albeluvisols", "Alisols", "Andosols","Anthrosols","Arenosols","Calcisols","Cambisols","Chernozems","Cryosols","Durisols","Ferralsols","Fluvisols","Gleysols","Gypsisols","Histosols","Kastanozems","Leptosols","Lixisols","Luvisols","Nitisols","Phaeozems","Planosols","Plinthosols","Podzols","Regosols","Solonchaks","Solonetz","Stagnosols","Technosols","Umbrisols","Vertisols")
  # set up first number to 0, if not the index when catNames equal to list it is not respect
  legend <- data.frame(agg=rg, int=0:31)
  
  
 }
 
 if(varn=="TAXOUSDA"){
  lut <- read.csv(file = "./LUT_USDA.csv")
  rg <- c("Alfisols","Andisols","Aridisols","Entisols","Gelisols","Histosol","Inceptisols" ,"Mollisols","Oxisol","Spodosols" ,"Ultisols","Vertisols")
  legend <- data.frame(agg=rg, int=0:11)
  
 }
 dum <- m@data
 reclass <- join(data.frame(level1=names(dum)), lut, by="level1")
 reclass$level3 <- paste(reclass$level1, reclass$agg, sep=" ")
 names(m@data) <- reclass$level3
 predicted <- aggregateRGclass(p=m@data, list=unique(reclass$agg))
 predicted <- join(predicted, lut, by="agg", match="first")
 m@data <- predicted
 writeGDAL(m["int"], out.c, drivername = "Gtiff",  mvFlag="-9999", options="COMPRESS=DEFLATE", catNames=list(paste((legend$agg))))
}
