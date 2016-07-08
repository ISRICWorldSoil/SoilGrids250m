#### initialization
setwd("E:/SoilGrids/")
# setwd("Y:/Maria/SOILGRIDS")
library(plyr)
library(xlsx)
library(mda)
library(plyr)
source("HighLev_classes_FUNCTIONS.R")

#### set variables
extent <- c("AF", "EU", "NA", "SA","OC")
varn.lst <- c("TAXNWRB", "TAXOUSDA")
var.lst <- c("WRB", "USDA")

#### calculate per tile
for (y in 1:length(extent)) {
 
 tiles.lst <- list.files(path=(paste("./tiles", sep="")), pattern=glob2rx(paste0(extent[y], "*")), full.names=FALSE) ## glob2rx("*.tif$"))

  for (c in 1:length(tiles.lst)) {
   tile <- tiles.lst[c]
 
    for (u in 1:length(varn.lst)) {
    varn <- varn.lst[u]
    test <- aggregateRG_Maps(i=tile, in.path=tile, out.path=getwd(), varn=varn)
        }
           }
              }
# end of script;
