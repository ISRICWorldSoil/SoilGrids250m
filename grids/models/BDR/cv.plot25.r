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
data(R_pal, package = "plotKML")
val.c <- NULL
plotList <- list(NULL)
log.flag <- 0
j<-1
for(j in 1:2)
{

  ## derive ME & RMSE:
    if(log.flag == 1)
    {
        cv.lst[[j]]$rf.predlog <- cv.lst[[j]]$rf.pred
        cv.lst[[j]]$measlog <- cv.lst[[j]]$meas
        cv.lst[[j]]$rf.pred <- expm1(cv.lst[[j]]$rf.predlog)
        cv.lst[[j]]$meas <- expm1(cv.lst[[j]]$measlog)
    }else if(log.flag == 0){
        cv.lst[[j]]$rf.predlog <- log1p(cv.lst[[j]]$rf.pred)
        cv.lst[[j]]$measlog <- log1p(cv.lst[[j]]$meas)
    }
    cv.lst[[j]]$rf.predcl <- toclass(cv.lst[[j]]$rf.pred, dclass)
    cv.lst[[j]]$meascl <- toclass(cv.lst[[j]]$meas, dclass)
    val.c[[j]] <- getcor(cv.lst[[j]]$meascl, cv.lst[[j]]$rf.predcl, length(dclass$class))
    tvar <- paste0(tbl$ATTRIBUTE_LABEL[j])
    tmp <- cbind(cv.lst[[j]]$rf.pred,cv.lst[[j]]$meas)
    tmp <- tmp[complete.cases(tmp), ]
    if(j==2) tmp <- tmp[tmp[,2]<250,]
    d.meas <- min(tmp[,2], na.rm=TRUE)
    out <- tmp[,1]+ifelse(d.meas==0, 1, d.meas)
    meas <- tmp[,2]+ifelse(d.meas==0, 1, d.meas)   
    if(j==1)
    {
        plotList[[j]] <- hexbinplot(out~meas,
                 colramp=colorRampPalette(R_pal[["bpy_colors"]]), main= "Absolute depth to bedrock (cm)",
                 xlab="measured", ylab="predicted",
                 type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8,
                 scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)),
                 xlim=range(meas, na.rm=TRUE), ylim=range(meas, na.rm=TRUE),
                 #xlim=c(0,10000), ylim=c(0,10000),
                 asp=1, xbins=25, density=40, panel=pfun, colorcut=c(0,0.002,0.01,0.03,0.07,0.15,0.25,0.5,1))
    }else if(j==2){
         
          plotList[[j]] <- hexbinplot(out~meas,
                 colramp=colorRampPalette(R_pal[["bpy_colors"]]), main = "Censored depth to bedrock (cm)",
                 xlab="measured", ylab="predicted",
                 type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8,
                 scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)),
                 xlim=range(meas, na.rm=TRUE), ylim=range(meas, na.rm=TRUE),
                 #xlim=c(0,10000), ylim=c(0,10000),
                 asp=1, xbins=25, density=40, panel=pfun
                 , colorcut=c(0,0.002,0.01,0.03,0.07,0.15,0.25,0.5,1)
                 )
    }
#    tmp <- tmp[tmp[,2]<1000,]     
#    plotList[[j]] <- hexbinplot(tmp[,1] ~ tmp[,2],
#             colramp=colorRampPalette(R_pal[["bpy_colors"]]), main = tvar,
#             xlab="measured", ylab="predicted",
#             type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8,
#             xlim=range(tmp[,2]), ylim=range(tmp[,2]),
#             asp=1, xbins=25, density=40, panel=pfun)
}
#plot(plotList[[1]])


do.call(grid.arrange, c(plotList, ncol=2))
dev.copy(png,paste0("./pics/cv_p", PC.flag, "_a", arti.flag, "_s", soil.flag, fit.name, "log.png"),  width = 800, height = 480, units = "px")
dev.off()

