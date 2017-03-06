## create soil profiles:
profs.f <- join(SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXGWRB")], HOR.s[,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT","SLTPPT","PHIHO5","CEC","BLD","ORCDRC")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TAXGWRB + SOURCEDB
## check if there are still some duplicates
profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
str(sel)
## no duplicates

## prepare a SITE and hor table table:
wsp22 <- list(sites=profs.f@site, horizons=profs.f@horizons)
str(wsp22)
lapply(as.list(wsp22$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp22$horizons), function(x){sum(!is.na(x))})
save(wsp22, file="../wsp22.rda", compress="xz")