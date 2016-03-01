# title         : cv.plot.R
# purpose       : plot cross validation, call by cv.n25.r;
# reference     :
# producer      : Prepared by W. Shangguan 
# address       : In Beijing.
# inputs        : 
# outputs       : png files ;
# remarks 1     : Takes ca 1.2 hrs to run with 5 cpus for randomforest in the defualt setting in use 
library(hexbin)
library(lattice)
library(gridExtra)

a.dir <- "/home/shang009/big/"# dir of the project
m.dir <- paste0(a.dir, "/soildepth2/cv")
setwd(m.dir)
t.vars <- c("BDRICM", "BDRLOG", "BDTICM")
source("../code/plot_hexbin.R")

plt.names <- c("Depth to bedrock (up to 250 cm)", "Occurrence of the R horizon", "Absolute depth to bedrock (in cm)")
names(plt.names) = t.vars
breaks.lst <- list(c(seq(0,250,length=50)), seq(0,1,length=50), seq(0,50000))
names(breaks.lst) = t.vars
plt.log <- c(FALSE, FALSE, TRUE)
names(plt.log) = t.vars

for(j in c(1,3)){
  plot_hexbin(j, breaks.lst[[t.vars[j]]], main=plt.names[t.vars[j]], in.file=paste0("CV_", t.vars[j], ".rda"), log.plot=plt.log[t.vars[j]])
}

load("./CV_BDRLOG.rda")
pred <- CV_BDRLOG$CV_residuals$Predicted
obs <- CV_BDRLOG$CV_residuals$Observed
#pred <- CV_BDRLOG$Predicted[,1]
#obs <- CV_BDRLOG$Observed[,1]
AUC <- performance( prediction(pred, obs), measure="auc")@y.values[[1]]
ROC <- performance( prediction(pred, obs), measure="tpr", x.measure = "fpr")

png(file = "plot_CV_BDRLOG.png", res = 150, width=850, height=850, type="cairo")
plot(ROC)
lines(c(0,1),c(0,1))
text(0.6,0.2,paste0("ROC curve AUC = ", signif(AUC,2)))
dev.off()

