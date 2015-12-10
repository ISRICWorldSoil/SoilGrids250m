## function for the ensemble predictions for Cross Validation
## by: Maria.RuiperezGonzales@wur.nl and Tom.Hengl@isric.org

## --------------------------------------------------------------
## Classes:
## --------------------------------------------------------------

## prediction error for predicting probs:
cv_factor <- function(formulaString, rmatrix, nfold, idcol, cpus=nfold){ 
  varn <- all.vars(formulaString)[1]
  sel <- dismo::kfold(rmatrix, k=nfold, by=rmatrix[,varn])
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting...", sep=""))
  ## run in parallel:
  if(missing(cpus)){ 
    cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
  }
  snowfall::sfInit(parallel=TRUE, cpus=cpus)
  snowfall::sfExport("ensemble.predict","idcol","formulaString","rmatrix","sel","varn","predict_parallel")
  snowfall::sfLibrary(package="ROCR", character.only=TRUE)
  snowfall::sfLibrary(package="nnet", character.only=TRUE)
  snowfall::sfLibrary(package="plyr", character.only=TRUE)
  snowfall::sfLibrary(package="randomForest", character.only=TRUE)
  out <- snowfall::sfLapply(1:nfold, function(j){predict_parallel(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol)})
  snowfall::sfStop()
  ## calculate totals per soil type
  N.tot <- Reduce("+", lapply(out, function(x){x[["n.l"]]}))
  mean.error <- Reduce("/", list(Reduce("+", lapply(out, function(x){x[["error.l"]]})), N.tot))
  error <- plyr::rbind.fill(lapply(out, function(x){x[["error"]]}))
  obs <- plyr::rbind.fill(lapply(out, function(x){ as.data.frame(x[["obs.pred"]][[1]])}))
  pred <- plyr::rbind.fill(lapply(out, function(x){ as.data.frame(x[["obs.pred"]][[2]])}))
  RMSE <- sapply(1:ncol(obs), function(c){mean(performance( prediction(pred[,c], obs[,c]), measure="tpr")@y.values[[1]])})
  AUC <- sapply(1:ncol(obs), function(c){performance( prediction(pred[,c], obs[,c]), measure="auc")@y.values[[1]]})
  cv.r <- list(obs, pred, error, data.frame(TPR=do.call(rbind, lapply(out, function(x){x[["rmse.l"]]})), AUC=do.call(rbind, lapply(out, function(x){x[["acc.l"]]}))), data.frame(ME=mean.error, TPR=RMSE, AUC=AUC, N=N.tot))
  names(cv.r) <- c("Observed", "Predicted", "CV_residuals", "Points", "Classes")
  cv.r[["Points"]][,idcol] <- error[,idcol]
  return(cv.r)
}

## factor-type vars:
ensemble.predict <- function(formulaString, s.train, s.test, MaxNWts = 19000, ...){ 
 gm1 <- nnet::multinom(formulaString, s.train, MaxNWts = MaxNWts)
 gm2 <- randomForest(formulaString, data=s.train, ...)
 probs1 <- predict(gm1, s.test, type="probs", na.action = na.pass) ## nnet
 probs2 <- predict(gm2, s.test, type="prob", na.action = na.pass) ## randomForest
 lt <- list(probs1[,gm1$lev], probs2[,gm1$lev])
 probs <- Reduce("+", lt) / length(lt)
 return(probs)
}   

## prediction in a loop:
predict_parallel <- function(j, sel, varn, formulaString, rmatrix, idcol){
  s.train <- rmatrix[!sel==j,]
  s.test <- rmatrix[sel==j,]
  N <- plyr::count(s.test[,varn])
  n.l <- N$freq
  probs <- ensemble.predict(formulaString=formulaString, s.train=s.train, s.test=s.test)
  names <- colnames(probs)
  obs <- data.frame(lapply(names, FUN=function(i){ifelse(s.test[, varn]==i, 1, 0)}))
  names(obs) = names
  obs.pred <- list(as.matrix(obs[,names]), probs[,names])
  error <- Reduce("-", obs.pred)
  error.l <- signif(colSums(error), 3)
  ## copy ID of the point
  error <- as.data.frame(error)
  error[,idcol] <- paste(s.test[,idcol])
  ## Accuracy for Binomial var [http://www.r-bloggers.com/evaluating-logistic-regression-models/]:
  pred.l <- lapply(1:nrow(obs.pred[[2]]), function(i){prediction(obs.pred[[2]][i,], obs.pred[[1]][i,])})
  ## prediction error per point
  rmse.l <- do.call(rbind, lapply(1:length(pred.l), function(i){mean(performance( pred.l[[i]], measure="tpr")@y.values[[1]])}))
  acc.l <- do.call(rbind, lapply(1:length(pred.l), function(i){performance( pred.l[[i]], measure="auc")@y.values[[1]]}))
  out <- list(n.l, obs.pred, error, error.l, rmse.l, acc.l)
  names(out) <- c("n.l", "obs.pred", "error", "error.l", "rmse.l", "acc.l")
  return(out)
}

## --------------------------------------------------------------
## Properties:
## --------------------------------------------------------------

predict_parallelP <- function(j, sel, varn, formulaString, rmatrix, idcol){ 
  s.train <- rmatrix[!sel==j,]
  s.test <- rmatrix[sel==j,]
  n.l <- dim(s.test)[1]
  gm <- randomForest(formulaString, data=s.train, na.action=na.omit)
  pred <- predict(gm, s.test, na.action = na.pass) 
  error <- list(s.test[,varn], as.data.frame(pred))
  error.l <- colSums(Reduce("-", error))      
  error.labs <- colSums(abs(Reduce("-", error)))      
  error.l2 <- colSums((Reduce("-", error))**2)
  obs.pred <- as.data.frame(error)
  names(obs.pred) = c("Observed", "Predicted")
  obs.pred[,idcol] <- s.test[,idcol]
  obs.pred$fold = j
  out <- list(obs.pred, error.l, error.labs, error.l2, n.l)
  names(out) = c("Observed", "Error", "AbsError", "SquaredError", "N")
  return(out)
}

cv_numeric <- function(formulaString, rmatrix, nfold, idcol, cpus=nfold){ 
  varn = all.vars(formulaString)[1]
  sel <- dismo::kfold(rmatrix, k=nfold)  
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting...", sep=""))
  if(nfold > nrow(rmatrix)){ 
    stop("'nfold' argument must not exceed total number of points") 
  }
  if(missing(cpus)){ 
    cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
  }
  snowfall::sfInit(parallel=TRUE, cpus=cpus)
  snowfall::sfExport("predict_parallelP","idcol","formulaString","rmatrix","sel","varn")
  snowfall::sfLibrary(package="plyr", character.only=TRUE)
  snowfall::sfLibrary(package="randomForest", character.only=TRUE)
  out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol)})
  snowfall::sfStop()
  ## calculate totals
  N <- Reduce("+", lapply(out, function(x){x[["N"]]}))
  mean.error <- Reduce("/", list(Reduce("+", lapply(out, function(x){x[["Error"]]})), N))
  mean.error2 <- Reduce("/", list(Reduce("+", lapply(out, function(x){x[["SquaredError"]]})), N))
  mae.error <- Reduce("/", list(Reduce("+", lapply(out, function(x){x[["AbsError"]]})), N))
  rmse <- sqrt(Reduce("/", list(Reduce("+", lapply(out, function(x){x[["SquaredError"]]})), N)))
  cv.r <- list(plyr::rbind.fill(lapply(out, function(x){x[["Observed"]]})), data.frame(ME=mean.error, MAE=mae.error, RMSE=rmse))
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}