profs.f <- join(site[sel.s,c("SOURCEID", "LONWGS84", "LATWGS84", "TAXGWRB")], horizon[sel.h,c("SOURCEID", "UHDICM", "LHDICM", "PHIHO5", "SNDPPT", "SLTPPT", "CLYPPT", "ORCDRC")], type='left')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TAXGWRB
## check if there are still some duplicates
profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sp::zerodist(sp)

## prepare a SITE and hor table table:
wsp5 <- list(sites=profs.f@site, horizons=profs.f@horizons)
str(wsp5)
lapply(as.list(wsp5$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp5$horizons), function(x){sum(!is.na(x))})
wsp5$sites$TAXGWRB <- as.character(wsp5$sites$TAXGWRB)
save(wsp5, file="../wsp5.rda", compress="xz")