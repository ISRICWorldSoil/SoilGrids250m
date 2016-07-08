#### initialization
setwd("E:/SoilGrids")
# setwd("Y:/Maria/SOILGRIDS")
library(plyr)
library(xlsx)
library(mda)
source("HighLev_classes_FUNCTIONS.R")

#### set variables
var.lst <- c("WRB", "USDA")

#### aggregation per system

for (u in 1:length(var.lst)) {

 varn <- var.lst[u]
 
 if(varn=="WRB"){
 ## create a lut to agg  WRB classes to reference groups
 download.file(paste0("http://gsif.isric.org/zipped/test.", varn, ".rda"), destfile=paste0(getwd(), paste0("/test.", varn, ".rda",sep=""))) 
 level1 <- c("Acrisols", "Albeluvisols", "Alisols", "Andosols","Anthrosols","Arenosols","Calcisols","Cambisols","Chernozems","Cryosols","Durisols","Ferralsols","Fluvisols","Gleysols","Gypsisols","Histosols","Kastanozems","Leptosols","Lixisols","Luvisols","Nitisols","Phaeozems","Planosols"  ,"Plinthosols","Podzols","Regosols","Solonchaks","Solonetz","Stagnosols","Technosols","Umbrisols","Vertisols")
 lut <- data.frame(int=0:31, agg=level1)
 load(paste0("./test.", varn, ".rda"))
 dum <-   test.WRB
 
 }

 if(varn=="USDA"){
 download.file(paste0("http://gsif.isric.org/zipped/test.", varn, ".rda"), destfile=paste0(getwd(), paste0("/test.", varn, ".rda",sep=""))) 
 ## create a lut for agg USDA classes 
 level1 <- c("Aqualfs","Cryalfs","Udalfs","Ustalfs","Xeralfs","Aquands","Cryands","Gelands","Torrands","Udands","Ustands","Vitrands","Xerands","Argids","Calcids","Cambids","Cryids","Durids","Gypsids","Salids","Aquents","Arents","Fluvents","Orthents","Psamments","Histels","Orthels","Turbels","Fibrists","Folists","Hemists","Saprists","Anthrepts","Aquepts","Cryepts","Gelepts","Ochrepts","Udepts","Ustepts","Xerepts","Albolls","Aquolls","Borolls","Cryolls","Gelolls","Rendolls","Udolls","Ustolls","Xerolls","Aquox","Perox","Torrox","Udox","Ustox","Aquods","Cryods","Gelods","Humods","Orthods","Aquults","Humults","Udults","Ustults","Xerults","Aquerts","Cryerts","Torrerts","Uderts","Usterts","Xererts")
  agg <- 
  c("Alfisols","Alfisols","Alfisols","Alfisols","Alfisols","Andisols","Andisols","Andisols","Andisols","Andisols","Andisols","Andisols","Andisols"   ,"Aridisols","Aridisols","Aridisols","Aridisols","Aridisols","Aridisols","Aridisols","Entisols","Entisols","Entisols","Entisols","Entisols","Gelisols","Gelisols","Gelisols","Histosol","Histosol","Histosol","Histosol","Inceptisols","Inceptisols","Inceptisols","Inceptisols","Inceptisols","Inceptisols","Inceptisols","Inceptisols","Mollisols","Mollisols","Mollisols","Mollisols","Mollisols","Mollisols","Mollisols","Mollisols","Mollisols","Oxisol","Oxisol","Oxisol","Oxisol","Oxisol","Spodosols","Spodosols","Spodosols","Spodosols","Spodosols","Ultisols","Ultisols","Ultisols","Ultisols","Ultisols","Vertisols","Vertisols","Vertisols","Vertisols","Vertisols","Vertisols")
 lut <- data.frame(level1=level1, agg=agg)
 load(paste0("./test.", varn, ".rda"))
 dum <-   test.USDA
 ## reclassify the names for the soil classes for USDA, that it is necessary for the function
 reclass <- join(data.frame(level1=names(dum[[2]])), lut, by="level1")
 reclass$level3 <- paste(reclass$level1, reclass$agg, sep=" ")
 names(dum[[2]]) <- reclass$level3
 reclass <- join(data.frame(level1=names(dum[[1]])), lut, by="level1")
 reclass$level3 <- paste(reclass$level1, reclass$agg, sep=" ")
 names(dum[[1]]) <- reclass$level3
}
 
predicted <- aggregateRGclass(dum[[2]],unique(lut$agg))
observed <- aggregateRGclass(dum[[1]], unique(lut$agg))

#### overall purity and confussion matrix
m <- data.frame(
 pred = predicted$agg,
 obs = observed$agg
)

# set factor levels
m$pred <- factor(m$pred, levels=unique(lut$agg) )
m$obs <- factor(m$obs, levels=unique(lut$agg) )

#### overall purity
cf <- confusion(m$pred, m$obs)
overall.pur <- round(sum(diag(cf))/sum(cf)*100,2)
overall.pur

#### map unit purity and class representation
# map unit purity (user's accuracy)
m$indicator <- NA
m$indicator <- ifelse(m$obs==m$pred, 1, 0)
purity <- tapply(m$indicator, INDEX=m$pred, FUN=mean)
(purity <- round(purity*100,1))

# class representation (producer's accuracy)
clasrep <- tapply(m$indicator, INDEX=m$obs, FUN=mean)
(clasrep <- round(clasrep*100,1))

# compile RSG specific purities
pur <- data.frame(
 mup = purity,
 cr = clasrep
)

## compile results in table
write.xlsx(cf, file = paste0("./cm", varn, ".xlsx"), sheetName="ConfusionMatrix")
write.xlsx(pur, file = paste0("./cm", varn, ".xlsx"), sheetName="RSG accuracies", append=TRUE)
write.xlsx(overall.pur, file = paste0("./cm", varn, ".xlsx"), sheetName="Overall accuracy", append=TRUE)
}
# end of script;