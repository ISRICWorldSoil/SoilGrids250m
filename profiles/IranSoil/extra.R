depths(Iransoil) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(Iransoil) <- ~ LONWGS84 + LATWGS84 + TAXNUSDA
## check if there are still some duplicates
Iransoil@site$SOURCEID[duplicated(Iransoil@site$SOURCEID)]
## spatial duplicates:
sp2 <- Iransoil@site
coordinates(sp2) <- ~ LONWGS84 + LATWGS84
sp::zerodist(sp2)
## 7 duplicates!

## prepare a SITE and hor table table:
wsp24 <- list(sites=Iransoil@site, horizons=Iransoil@horizons[,c("SNDPPT","SLTPPT","CLYPPT","PHIHO5","BLD","ORCDRC")])
str(wsp24)
wsp24$sites$SOURCEDB = "IranSoil"
lapply(as.list(wsp24$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp24$horizons), function(x){sum(!is.na(x))})
wsp24$sites$TAXNUSDA <- as.character(wsp24$sites$TAXNUSDA)
save(wsp24, file="../wsp24.rda", compress="xz")