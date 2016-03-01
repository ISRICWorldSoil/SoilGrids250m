## function for the ensemble predictions for Cross Validation (runs in parallel hence can be run on large data sets)
## by: Tom.Hengl@isric.org, Gerard.Heuvelink@wur.nl and Maria.RuiperezGonzales@wur.nl

list.of.packages <- c("nnet", "plyr", "ROCR", "randomForest", "plyr", "parallel", "psych", "mda", "h2o", "dismo", "grDevices", "snowfall", "hexbin", "lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## --------------------------------------------------------------
## Classes:
## --------------------------------------------------------------

## prediction error for predicting probs:
cv_factor <- function(formulaString, rmatrix, nfold, idcol, cpus=nfold, h2o){ 
  varn <- all.vars(formulaString)[1]
  sel <- dismo::kfold(rmatrix, k=nfold, by=rmatrix[,varn])
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting...", sep=""))
  ## run in parallel:
  if(h2o==TRUE){
    cpus = 1
  } else {
    if(missing(cpus)){ 
      cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
    }
  }
  if(cpus==1){
    out <- list()
    for(j in 1:nfold){ 
      out[[j]] <- predict_parallel(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, h2o=h2o)
    }
  } else {
    snowfall::sfInit(parallel=TRUE, cpus=cpus)
    snowfall::sfExport("ensemble.predict","idcol","formulaString","rmatrix","sel","varn","predict_parallel")
    snowfall::sfLibrary(package="ROCR", character.only=TRUE)
    snowfall::sfLibrary(package="nnet", character.only=TRUE)
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="randomForest", character.only=TRUE)
    out <- snowfall::sfLapply(1:nfold, function(j){predict_parallel(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, h2o=h2o)})
    snowfall::sfStop()
  }
  
  
  ## calculate totals per soil type
  N.tot <- plyr::rbind.fill(lapply(out, function(x){x[["n.l"]]}))
  N.tot <- colSums(N.tot, na.rm = T)
  
  ## mean error per soil type:
  mean.error <- plyr::rbind.fill(lapply(out, function(x){x[["error.l"]]}))
  mean.error <- colSums(mean.error)/N.tot
  error <- plyr::rbind.fill(lapply(out, function(x){x[["error"]]}))
  obs <- plyr::rbind.fill(lapply(out, function(x){ as.data.frame(x[["obs.pred"]][[1]])}))
  pred <- plyr::rbind.fill(lapply(out, function(x){ as.data.frame(x[["obs.pred"]][[2]])}))
  
  #get rid of NA in pred
  error <- error[which(!is.na(pred[,1])&!is.na(pred[,2])),]
  obs <- obs[which(!is.na(pred[,1])&!is.na(pred[,2])),]
  pred <- pred[which(!is.na(pred[,1])&!is.na(pred[,2])),]
  
  ## Get the most probable class:
  cl <- parallel::makeCluster(getOption("cl.cores", cpus))
  ranks.pred <- parallel::parApply(cl, pred, MARGIN=1, which.max)
  ranks.obs <- parallel::parApply(cl, obs, MARGIN=1, which.max)
  parallel::stopCluster(cl)
  
  ## derive confusion matrix:
  cf <- mda::confusion(names(pred)[ranks.pred], names(obs)[ranks.obs])
  c.kappa <- psych::cohen.kappa(cf)
  purity <- sum(diag(cf))/sum(cf)*100  
  ## Accuracy for Binomial var [http://www.r-bloggers.com/evaluating-logistic-regression-models/]: 
  library(ROCR)
  ROC <- performance( prediction(pred[,1], obs[,1]), measure="tpr", x.measure = "fpr")
  TPR <- sapply(1:ncol(obs), function(c){mean(performance( prediction(pred[,c], obs[,c]), measure="tpr")@y.values[[1]])})
  AUC <- sapply(1:ncol(obs), function(c){performance( prediction(pred[,c], obs[,c]), measure="auc")@y.values[[1]]})
  cv.r <- list(obs, pred, error, data.frame(ME=mean.error, TPR=TPR, AUC=AUC, N=N.tot), cf, c.kappa, purity, ROC)
  names(cv.r) <- c("Observed", "Predicted", "CV_residuals", "Classes", "Confusion.matrix", "Cohen.Kappa", "Purity", "ROC")
  return(cv.r)
}

## factor-type vars:
ensemble.predict <- function(formulaString, s.train, s.test, MaxNWts = 19000, h2o=h2o, ...){ 
  ## drop empty levels to avoid errors:
  s.train[,all.vars(formulaString)[1]] <- droplevels(s.train[,all.vars(formulaString)[1]])
  
  if(h2o==TRUE){
    ## select only complete point pairs
    train.hex <- as.h2o(s.train[complete.cases(s.train[,all.vars(formulaString)]),all.vars(formulaString)], destination_frame = "train.hex")    
    test.hex <- as.h2o(s.test[,all.vars(formulaString)], destination_frame = "test.hex")
    gm2 <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString)), training_frame=train.hex)
    probs2 <- as.data.frame(h2o.predict(gm2, test.hex, type="probs", na.action=na.pass))[, c("p0","p1")]
    colnames(probs2) <- substring(names(probs2),2)
    probs <- probs2
    h2o.removeAll()
  } else {
    gm1 <- nnet::multinom(formulaString, s.train, MaxNWts = MaxNWts)
    probs1 <- predict(gm1, s.test, type="probs", na.action = na.pass) ## nnet
    gm2 <- randomForest(formulaString, data=s.train, na.action= na.omit)
    probs2 <- predict(gm2, s.test, type="prob", na.action = na.pass) ## randomForest 
    if(length(levels(s.train$BDRLOG))==2)
    {
      probs1 <- cbind(1-probs1,probs1)
      colnames(probs1) <- gm1$lev
    }
    lt <- list(probs1[,gm1$lev], probs2[,gm1$lev])
    probs <- Reduce("+", lt) / length(lt)  
  } 
  return(probs)
}   

## ensemble prediction in parallel (for parallelization):
predict_parallel <- function(j, sel, varn, formulaString, rmatrix, idcol,h2o){
  s.train <- rmatrix[!sel==j,]
  s.test <- rmatrix[sel==j,]
  n.l <- plyr::count(s.test[,varn])
  n.l <- data.frame(matrix(n.l$freq, nrow=1, dimnames = list(1, paste(n.l$x))))
  probs <- ensemble.predict(formulaString=formulaString, s.train=s.train, s.test=s.test, h2o=h2o)
  names <- colnames(probs)
  obs <- data.frame(lapply(names, FUN=function(i){ifelse(s.test[, varn]==i, 1, 0)}))
  names(obs) = names
  obs.pred <- list(as.matrix(obs[,names]), probs[,names])
  error <- Reduce("-", obs.pred)
  error.l <- as.data.frame(t(signif(colSums(error, na.rm=T), 3)))
  ## copy ID of the point
  error <- as.data.frame(error)
  error[,idcol] <- paste(s.test[,idcol])
  out <- list(n.l, obs.pred, error, error.l)
  names(out) <- c("n.l", "obs.pred", "error", "error.l")
  return(out)
}

## --------------------------------------------------------------
## Properties:
## --------------------------------------------------------------

## predict soil properties in parallel:
predict_parallelP <- function(j, sel, varn, formulaString, rmatrix, idcol, h2o){ 
  s.train <- rmatrix[!sel==j,]
  s.test <- rmatrix[sel==j,]
  n.l <- dim(s.test)[1]
  if(h2o==TRUE){
    ## select only complete point pairs
    train.hex <- as.h2o(s.train[complete.cases(s.train[,all.vars(formulaString)]),all.vars(formulaString)], destination_frame = "train.hex")
    gm1 <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString)), training_frame=train.hex) 
    gm2 <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString)), training_frame=train.hex)
    test.hex <- as.h2o(s.test[,all.vars(formulaString)], destination_frame = "test.hex")
    v1 <- as.data.frame(h2o.predict(gm1, test.hex, na.action=na.pass))$predict
    gm1.w = gm1@model$training_metrics@metrics$r2
    v2 <- as.data.frame(h2o.predict(gm2, test.hex, na.action=na.pass))$predict
    gm2.w = gm2@model$training_metrics@metrics$r2
    ## mean prediction based on accuracy:
    pred <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
    gc()
    h2o.removeAll()
  } else {
    gm <- randomForest(formulaString, data=s.train, na.action=na.omit)
    pred <- predict(gm, s.test, na.action = na.pass) 
  }
  obs.pred <- as.data.frame(list(s.test[,varn], pred))
  names(obs.pred) = c("Observed", "Predicted")
  obs.pred[,idcol] <- s.test[,idcol]
  obs.pred$fold = j
  return(obs.pred)
}

cv_numeric <- function(formulaString, rmatrix, nfold, idcol, cpus=nfold, h2o=FALSE, Log=FALSE){ 
  varn = all.vars(formulaString)[1]
  sel <- dismo::kfold(rmatrix, k=nfold)  
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting...", sep=""))
  if(nfold > nrow(rmatrix)){ 
    stop("'nfold' argument must not exceed total number of points") 
  }
  if(h2o==TRUE){
    cpus = 1
  } else {
    if(missing(cpus)){ 
      cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
    }
  }
  if(cpus==1){
    out <- list()
    for(j in 1:nfold){ 
      out[[j]] <- predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, h2o=h2o)
    }
  } else {
    snowfall::sfInit(parallel=TRUE, cpus=cpus)
    snowfall::sfExport("predict_parallelP","idcol","formulaString","rmatrix","sel","varn","h2o")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="randomForest", character.only=TRUE)
    out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, h2o=h2o)})
    snowfall::sfStop()
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE) 
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  #R.squared = 1-sum((out$Observed - out$Predicted)^2, na.rm=TRUE)/(var(out$Observed, na.rm=TRUE)*sum(!is.na(out$Observed)))
  R.squared = 1-var(out$Observed - out$Predicted, na.rm=TRUE)/var(out$Observed, na.rm=TRUE)
  if(Log==TRUE){
    ## If the variable is log-normal then logR.squared is probably more correct
    logRMSE = sqrt(mean((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE))
    #logR.squared = 1-sum((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE)/(var(log1p(out$Observed), na.rm=TRUE)*sum(!is.na(out$Observed)))
    logR.squared = 1-var(log1p(out$Observed) - log1p(out$Predicted), na.rm=TRUE)/var(log1p(out$Observed), na.rm=TRUE)
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, R.squared=R.squared, logRMSE=logRMSE, logR.squared=logR.squared)) 
  } else {
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, R.squared=R.squared))
  }
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}

## correlation plot:
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.abline(0,1,lty=1,lw=2,col="black")
  panel.abline(0+RMSE,1,lty=2,lw=2,col="black")
  panel.abline(0-RMSE,1,lty=2,lw=2,col="black")
}
