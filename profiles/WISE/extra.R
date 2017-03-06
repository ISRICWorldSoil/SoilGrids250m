
profs.f <- join(SITE.sm[,c("SOURCEID", "LONWGS84", "LATWGS84", "TIMESTRR", "TAXGWRB", "TAXNUSDA")], HORIZON[sel.h,c("SOURCEID","UHDICM","LHDICM","CLYPPT","SNDPPT", "SLTPPT", "CRFVOL","PHIHO5","PHIKCL","BLD","ORCDRC")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TIMESTRR + TAXGWRB + TAXNUSDA + BDRICM
## check if there are still some duplicates
#profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]

## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sel <- sp::zerodist(sp)
str(sel)
## TH: 2504 duplicate points?!

# prepare a SITE and hor table table:
wsp8 <- list(sites=profs.f@site, horizons=profs.f@horizons)
str(wsp8)
wsp8$sites$TAXGWRB <- as.character(wsp8$sites$TAXGWRB)
wsp8$sites$TAXNUSDA <- as.character(wsp8$sites$TAXNUSDA)
lapply(as.list(wsp8$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp8$horizons), function(x){sum(!is.na(x))})
save(wsp8, file="../wsp8.rda", compress="xz")