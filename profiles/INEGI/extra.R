

inegi <- read.csv("INEGIv1_1km.csv")
## takes 10 secs to read
str(inegi)
summary(inegi$X_COORD)
summary(inegi$Y_COORD)
## some coordiantes have been messed up:
inegi[!(inegi$X_COORD>0&inegi$Y_COORD>0),]
sp.xy <- inegi[,c("IDENTIFI", "X_COORD", "Y_COORD")]
coordinates(sp.xy) <- ~ X_COORD + Y_COORD
proj4string(sp.xy) <- CRS(mx.csy)
sp.ll <- as.data.frame(spTransform(sp.xy, CRS("+proj=longlat +datum=WGS84")))
str(sp.ll)
sp.ll[which(!(inegi$X_COORD>0&inegi$Y_COORD>0)), c("x", "y")] <- NA
inegi$LONWGS84 <- sp.ll$x
inegi$LATWGS84 <- sp.ll$y

## ID column:
plyr:::nunique(inegi$IDENTIFI)
inegi$SOURCEID <- as.factor(paste("MX", inegi$IDENTIFI, sep="_"))
## check if there are duplicates:
length(levels(inegi$SOURCEID))

## get horizon colours:
inegi.cols <- inegi[,c("SOURCEID","A_DENOM","A_COLORH","A_COLORS","E_COLOR","B_COLORH")]
inegi.cols <- rename(inegi.cols, c("A_DENOM"="SYMHOR_A","A_COLORH"="MCOMNS_A","A_COLORS"="MCDMNS_A","E_COLOR"="MCOMNS_B","B_COLORH"="MCOMNS_C"))
horizons <- getHorizons(inegi.cols, idcol="SOURCEID", sel=c("SYMHOR","MCOMNS","MCDMNS"))

## rename columns:
inegi.f <- rename(inegi, c("HLIMSUPE"="UHDICM", "HLIMINFE"="LHDICM", "FAO68"="TAXGWRB", "ARCILLA"="CLYPPT", "LIMO"="SLTPPT", "ARENA"="SNDPPT", "PH"="PHIHO5", "MO"="ORCDRC", "CIC"="CEC"))

## subset to complete data:
sel.c <- !is.na(inegi.f$LHDICM)&!is.na(inegi.f$UHDICM)&!is.na(inegi.f$LONWGS84)&!is.na(inegi.f$LATWGS84)& !(rowSums(inegi.f[,c("UHDICM","LHDICM")])==0)
inegi.f <- inegi.f[sel.c,]
## filter missing values:
inegi.f$PHIHO5 <- ifelse(inegi.f$PHIHO5<2, NA, inegi.f$PHIHO5)
inegi.f$SNDPPT <- ifelse(inegi.f$SNDPPT==0, NA, inegi.f$SNDPPT)
inegi.f$SLTPPT <- ifelse(inegi.f$SLTPPT==0, NA, inegi.f$SLTPPT)
inegi.f$CLYPPT <- ifelse(inegi.f$CLYPPT==0, NA, inegi.f$CLYPPT)
inegi.f$ORCDRC <- inegi.f$ORCDRC * 10/1.724
summary(inegi.f$ORCDRC)
## View(inegi.f)

## convert to SPC class:
depths(inegi.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(inegi.f) <- ~ LONWGS84 + LATWGS84 + TAXGWRB
## check if there are still some duplicates
sp <- inegi.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sp::zerodist(sp)
## no spatial duplicates  :)
plyr:::nunique(inegi.f@site$SOURCEID)
plyr:::nunique(paste(sp@coords[,1], sp@coords[,2], sep=""))
inegi.f@site$SOURCEDB = "MexSPDB"
## Depth to bedrock i.e. 'R' horizon:
sel.r1 <- grep(pattern="x", inegi.f@horizons$LIM_ROCA, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="x", inegi.f@horizons$LIM_REG, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="R", inegi.f@horizons$HSIMBOLO, ignore.case=FALSE, fixed=FALSE)
sel.r <- unique(sel.r1, sel.r2, sel.r3)
inegi.f@horizons$BDRICM <- NA
inegi.f@horizons$BDRICM[sel.r] <- inegi.f@horizons$PROFUNDI[sel.r]
bdr.d <- aggregate(inegi.f@horizons$BDRICM, list(inegi.f@horizons$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
m <- merge(inegi.f@site, bdr.d, all.y=FALSE)
m$BDRICM <- ifelse(m$BDRICM<2, NA, m$BDRICM)
summary(m$BDRICM)


# prepare the dataset:
wsp14 <- list(sites=m, horizons=inegi.f@horizons[,c("SOURCEID","UHDICM","LHDICM","PHIHO5","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC")])
str(wsp14)
wsp14$sites$TAXGWRB <- as.character(wsp14$sites$TAXGWRB)
lapply(as.list(wsp14$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp14$horizons), function(x){sum(!is.na(x))})
save(wsp14, file="../wsp14.rda", compress="xz")


## plotting horizon data:
inegi.xy <- inegi[,c("IDENTIFI", "X_COORD", "Y_COORD", "NHORIZON", "HLIMSUPE", "HLIMINFE", "PH", "ARCILLA", "MO", "FAO68")]
inegi.xy <- inegi.xy[!inegi.xy$NHORIZON==0,]
for(j in c("PH", "ARCILLA")){ inegi.xy[inegi.xy[,j]==0&!is.na(inegi.xy[,j]),] = NA }
str(inegi.xy)
inegi.xy <- inegi.xy[!is.na(inegi.xy$IDENTIFI),]
inegi.xy$IDENTIFI <- as.factor(paste("ID", inegi.xy$IDENTIFI, sep=""))
summary(inegi.xy$FAO68)
inegi.spc <- inegi.xy
## convert to geosamples:
depths(inegi.spc) <- IDENTIFI ~ HLIMSUPE + HLIMINFE
site(inegi.spc) <- ~ X_COORD + Y_COORD + FAO68
coordinates(inegi.spc) <- ~X_COORD + Y_COORD
proj4string(inegi.spc) <- CRS(mx.csy)
## convert to geosamples:
inegi.geo <- as.geosamples(inegi.spc)
save(inegi.geo, file="inegi.geo.rda")

inegi.xy1 <- inegi.xy[inegi.xy$NHORIZON==1,]
inegi.xy1 <- inegi.xy1[!(inegi.xy1$X_COORD==0|inegi.xy1$Y_COORD==0),]
coordinates(inegi.xy1) <- ~ X_COORD + Y_COORD
proj4string(inegi.xy1) <- CRS(mx.csy)
save(inegi.xy1, file="inegi.xy1.rda")
#diff(t(inegi.xy1@bbox))
#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#kml(inegi.xy, shape=shape, colour=PH, colour_scale=R_pal[["pH_pal"]], points_names="", balloon=TRUE)

# end of script;