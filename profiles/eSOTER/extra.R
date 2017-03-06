
## create soil profiles:
profs.f <- join(SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","TAXGWRB","TAXNUSDA")], HOR.s[,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHO5","CEC","BLD","ORCDRC")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TAXGWRB + TAXNUSDA + SOURCEDB + TIMESTRR
## check if there are still some duplicates
profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
str(sel)
## 55 duplicates!

## prepare a SITE and hor table table:
wsp21 <- list(sites=sites.m, horizons=profs.f@horizons)
str(wsp21)
lapply(as.list(wsp21$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp21$horizons), function(x){sum(!is.na(x))})
save(wsp21, file="../wsp21.rda", compress="xz")
