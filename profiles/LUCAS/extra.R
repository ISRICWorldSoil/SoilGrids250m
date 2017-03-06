## convert to SPC class:
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TIMESTRR
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
str(sel)

## add missing columns:
profs.f@site$SOURCEDB <- "LUCAS"

## export:
wsp11 <- list(sites=profs.f@site, horizons=profs.f@horizons[,c("SOURCEID","UHDICM","LHDICM","CRFVOL","PHIHO5","PHIKCL","ORCDRC","SNDPPT","SLTPPT","CLYPPT","CEC")])
str(wsp11)
lapply(as.list(wsp11$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp11$horizons), function(x){sum(!is.na(x))})
save(wsp11, file="../wsp11.rda", compress="xz")