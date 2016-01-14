## Function to fit models using nfold CV:
cv_nfold <- function(nf, rmat, nfold.x, fm.g, tvar, sub.N, cell.size, n){
   set.seed(10002)
       #nf <-1
       t.rmat <- rmat[!(nfold.x==nf),]
       v.rmat <- rmat[nfold.x==nf,]
       v.rmat <- v.rmat[!is.na(v.rmat[,1]),]
       coordinates(t.rmat) <- ~LONWGS84 +LATWGS84.1
       proj4string(t.rmat) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
       #get rid of spatial clustering
       t.rmat <-  as.data.frame(GSIF::sample.grid(t.rmat, cell.size = c(cell.size, cell.size), n = n)$sub)
       ## To speed things up, select only fixed random sub-sample:
       if(sub.N > nrow(t.rmat)){ sub.N <- nrow(t.rmat) }
       sel.s  <- sample(1:nrow(t.rmat), sub.N)
       t.rmat <- t.rmat[sel.s, ]
    vs.rfs<-NULL
#   options(rf.cores=11, mc.cores=11)
#   vs.rfs <- var.select(fm.g, t.rmat, ntree =100, nsplit = 10)
#   #save(vs.rfs, file = paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
#   #load(paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))  
#   fm.g <- as.formula(paste0(tvar," ~ ",  paste(vs.rfs$topvars, collapse="+"))) 
   ####Random forests modelling
   dfs <- t.rmat[,all.vars(fm.g)]
   dfs.hex <- as.h2o(dfs[complete.cases(dfs),], destination_frame = "dfs.hex")
   m.grid <- as.h2o(v.rmat,  destination_frame="m.grid")
  fname <- paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".txt") 
  if(file.exists(fname))
  {
   mrf_path <- read.table(fname, as.is = T)$x
   m.rf <- h2o.loadModel(mrf_path)  
   rf.pred <- as.data.frame(h2o.predict(m.rf, m.grid, na.action=na.pass))$predict
  }else
  {
       try( m.rf <- h2o.randomForest(y=1, x=2:length(all.vars(fm.g)), training_frame=dfs.hex) ) 
       if(!class(.Last.value)[1]=="try-error"&!is.null(m.rf)){
         mrf_path = h2o.saveModel(m.rf, path="./cv/", force=TRUE)
         write.table(mrf_path, file=fname)
         rf.pred <- as.data.frame(h2o.predict(m.rf, m.grid, na.action=na.pass))$predict
       } else {
         rf.pred <- rep(NA, nrow(v.rmat))
       } 
   } 
  
    rm(m.rf)
    gc()
 if(m.flag == 1)
 {  
    ###deep learning
    fname <- paste0("./cv/", tvar, "m.dl_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".txt") 
    if(file.exists(fname))
    {
       mdl_path <- read.table(fname, as.is = T)$x
       m.dl <- h2o.loadModel(mdl_path) 
       dl.pred <- as.data.frame(h2o.predict(m.dl,m.grid, na.action=na.pass))$predict 
    }else
    {
         try(m.dl <- h2o.deeplearning(y=1, x=2:length(all.vars(fm.g)), training_frame=dfs.hex))
         if(!class(.Last.value)[1]=="try-error"&!is.null(m.dl)){
          mdl_path = h2o.saveModel(m.dl, path="./cv/", force=TRUE)
          write.table(mdl_path, file=fname)
          dl.pred <- as.data.frame(h2o.predict(m.dl,m.grid, na.action=na.pass))$predict       
        } else {
          dl.pred <- rep(NA, nrow(v.rmat))
        }
    }
    rm(m.dl)   
    if(tvar == "BDRLOG")
    {
        rk.pred <- NA 
    }else {
        # Fit linear regression model using stepwise selection
          fname <- paste0("./cv/",tvar, "m.rk_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda") 
        if(file.exists(fname))
        {
            load(fname)  
            try(rk.pred <- predict(m.rk, v.rmat)) 
        }else
        {
            m.rk <- lm(fm.g, t.rmat)
            m.rk <- step(m.rk)
            save(m.rk, file = fname)
            try(rk.pred <- predict(m.rk, v.rmat)) 
        } 
        rm(m.rk)      
        gc()
    }
 } 
    out.df <- data.frame(  rf.pred = rf.pred,
    meas=eval(fm.g[[2]], v.rmat), longitude=v.rmat$LONWGS84, latitude=v.rmat$LATWGS84)
    if(m.flag == 1)
    {
        out.df$rk.pred  = rk.pred
        out.df$dl.pred  = dl.pred
    } 
   return(out.df)
}
