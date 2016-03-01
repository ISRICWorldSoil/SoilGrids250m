library(plyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
library(snowfall)
library(utils)
library(plotKML)
library(GSIF)

a.dir <- "/home/shang009/big/"# dir of the project
m.dir <- paste0(a.dir, "/soildepth2")
setwd(m.dir)
source("./code/cv_functions.R")
source(paste0(a.dir, "soildepth/code/head/functions.r"))

t.vars <- c("BDRICM", "BDRLOG", "BDTICM")
des <- read.csv("SoilGrids250m_COVS250m.csv")
pr.lst <- des$WORLDGRIDS_CODE
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ LATWGS84 + ', paste(pr.lst, collapse="+")))})

## H2O package more suited for large data (http://www.r-bloggers.com/benchmarking-random-forest-implementations/)
library(h2o)
## reset to use all cores:
localH2O = h2o.init(nthreads = -1, assertion = FALSE)

###takes several minutes
ovA <- read.csv("ov.BDR_SoilGrids250m.csv")
ovA$BDRLOG <-as.factor(ovA$BDRLOG)



####to speed up, only use a subset of the data
#ov<-ovA[1:10000,]
sp <- ovA[,c("X.1","LONWGS84", "LATWGS84")]
sp <- sp[!is.na(sp$LONWGS84)& !is.na(sp$LATWGS84),]
coordinates(sp) <-  ~ LONWGS84+LATWGS84
sp1 <- subsp(sp,c(0.01,0.01),1)
ov<- ovA[sp1$X.1,]
#ov<- ovA
rm(ovA, sp, sp1)




cat("Results of Cross-validation:\n\n", file="resultsCV_BDR.txt")
cv_lst <- rep(list(NULL), length(t.vars))
#for(j in 1:length(t.vars)){
for(j in 2:2){
  if(!file.exists(paste0("CV_", t.vars[j], ".rda"))){
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file="resultsCV_BDR.txt", append=TRUE)
    cat("\n", file="resultsCV_BDR.txt", append=TRUE)
    if(t.vars[j] == "BDRLOG")
    {
      
      cv_lst[[j]] <- cv_factor(formulaString.lst[[j]], rmatrix=ov, nfold=10, idcol="SOURCEID", cpus=2, h2o=TRUE)
      sink(file="resultsCV_BDR.txt", append=TRUE, type="output")
      print(cv_lst[[j]]$Classes)
      cat("\nConfusion.matrix:\n", file="resultsCV_BDR.txt", append=TRUE)
      print(cv_lst[[j]]$Confusion.matrix)
      cat("\n\nCohen.Kappa:\n", file="resultsCV_BDR.txt", append=TRUE)
      print(cv_lst[[j]]$Cohen.Kappa)
      cat("\n\npurity:\n", file="resultsCV_BDR.txt", append=TRUE)
      print(cv_lst[[j]]$Purity)
    }else{
      cv_lst[[j]] <- cv_numeric(formulaString.lst[[j]], rmatrix=ov[!is.na(ov[, t.vars[j]]),], nfold=10, idcol="SOURCEID", h2o=TRUE, Log=TRUE)
      sink(file="resultsCV_BDR.txt", append=TRUE, type="output")
      print(cv_lst[[j]]$Summary)
    }
    
    cat("\n", file="resultsCV_BDR.txt", append=TRUE)
    sink()
    assign(paste0("CV_", t.vars[j]), cv_lst[[j]])
    save(list=paste0("CV_", t.vars[j]), file=paste0("CV_", t.vars[j], ".rda"))
  }
}

h2o.shutdown(prompt = F)
