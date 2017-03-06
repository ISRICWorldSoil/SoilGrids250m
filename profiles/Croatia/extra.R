hor1 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR1_OZN, UHDICM=as.numeric(site$HOR1_GOR), LHDICM=as.numeric(site$HOR1_DON), CRFVOL=as.numeric(site$HOR1_MKP),  PHIHO5=as.numeric(site$HOR1_PH1), PHIKCL=as.numeric(site$HOR1_PH2), SNDPPT=as.numeric(site$HOR1_MSP)+as.numeric(site$HOR1_MKP), SLTPPT=as.numeric(site$HOR1_MP), CLYPPT=as.numeric(site$HOR1_MG), ORCDRC=as.numeric(site$HOR1_HUM)*10/1.724)
hor2 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR2_OZN, UHDICM=as.numeric(site$HOR2_GOR), LHDICM=as.numeric(site$HOR2_DON), CRFVOL=as.numeric(site$HOR2_MKP),  PHIHO5=as.numeric(site$HOR2_PH1), PHIKCL=as.numeric(site$HOR2_PH2), SNDPPT=as.numeric(site$HOR2_MSP)+as.numeric(site$HOR2_MKP), SLTPPT=as.numeric(site$HOR2_MP), CLYPPT=as.numeric(site$HOR2_MG), ORCDRC=as.numeric(site$HOR2_HUM)*10/1.724)
hor3 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR3_OZN, UHDICM=as.numeric(site$HOR3_GOR), LHDICM=as.numeric(site$HOR3_DON), CRFVOL=as.numeric(site$HOR3_MKP),  PHIHO5=as.numeric(site$HOR3_PH1), PHIKCL=as.numeric(site$HOR3_PH2), SNDPPT=as.numeric(site$HOR3_MSP)+as.numeric(site$HOR3_MKP), SLTPPT=as.numeric(site$HOR3_MP), CLYPPT=as.numeric(site$HOR3_MG), ORCDRC=as.numeric(site$HOR3_HUM)*10/1.724)
hor4 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR4_OZN, UHDICM=as.numeric(site$HOR4_GOR), LHDICM=as.numeric(site$HOR4_DON), CRFVOL=as.numeric(site$HOR4_MKP), PHIHO5=as.numeric(site$HOR4_PH1), PHIKCL=as.numeric(site$HOR4_PH2), SNDPPT=as.numeric(site$HOR4_MSP)+as.numeric(site$HOR4_MKP), SLTPPT=as.numeric(site$HOR4_MP), CLYPPT=as.numeric(site$HOR4_MG), ORCDRC=as.numeric(site$HOR4_HUM)*10/1.724)
hor5 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR5_OZN, UHDICM=as.numeric(site$HOR5_GOR), LHDICM=as.numeric(site$HOR5_DON), CRFVOL=as.numeric(site$HOR5_MKP), PHIHO5=as.numeric(site$HOR5_PH1), PHIKCL=as.numeric(site$HOR5_PH2), SNDPPT=as.numeric(site$HOR5_MSP)+as.numeric(site$HOR5_MKP), SLTPPT=as.numeric(site$HOR5_MP), CLYPPT=as.numeric(site$HOR5_MG), ORCDRC=as.numeric(site$HOR5_HUM)*10/1.724)
hor6 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR6_OZN, UHDICM=as.numeric(site$HOR6_GOR), LHDICM=as.numeric(site$HOR6_DON), CRFVOL=as.numeric(site$HOR6_MKP), PHIHO5=as.numeric(site$HOR6_PH1), PHIKCL=as.numeric(site$HOR6_PH2), SNDPPT=as.numeric(site$HOR6_MSP)+as.numeric(site$HOR6_MKP), SLTPPT=as.numeric(site$HOR6_MP), CLYPPT=as.numeric(site$HOR6_MG), ORCDRC=as.numeric(site$HOR6_HUM)*10/1.724)
hor7 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR7_OZN, UHDICM=as.numeric(site$HOR7_GOR), LHDICM=as.numeric(site$HOR7_DON), CRFVOL=as.numeric(site$HOR7_MKP), PHIHO5=as.numeric(site$HOR7_PH1), PHIKCL=as.numeric(site$HOR7_PH2), SNDPPT=as.numeric(site$HOR7_MSP)+as.numeric(site$HOR7_MKP), SLTPPT=as.numeric(site$HOR7_MP), CLYPPT=as.numeric(site$HOR7_MG), ORCDRC=as.numeric(site$HOR7_HUM)*10/1.724)
hor8 <- data.frame(SOURCEID=site$SOURCEID, HORIZON=site$HOR8_OZN, UHDICM=as.numeric(site$HOR8_GOR), LHDICM=as.numeric(site$HOR8_DON), CRFVOL=as.numeric(site$HOR8_MKP), PHIHO5=as.numeric(site$HOR8_PH1), PHIKCL=as.numeric(site$HOR8_PH2), SNDPPT=as.numeric(site$HOR8_MSP)+as.numeric(site$HOR8_MKP), SLTPPT=as.numeric(site$HOR8_MP), CLYPPT=as.numeric(site$HOR8_MG), ORCDRC=as.numeric(site$HOR8_HUM)*10/1.724)
## bind together
horizon <- do.call(rbind, list(hor1, hor2, hor3, hor4, hor5, hor6, hor7, hor8))

profs.f <- join(site[sel.s,c("SOURCEID", "LONWGS84", "LATWGS84", "TIMESTRR", "TAXGWRB", "SOURCEDB", "BDRICM")], horizon[sel.h,c("SOURCEID","UHDICM","LHDICM","PHIHOX","PHIKCL","SNDPPT","SLTPPT","CLYPPT","ORCDRC","CRFVOL")], type='inner')
depths(profs.f) <- SOURCEID ~ UHDICM + LHDICM
# extract site data
site(profs.f) <- ~ LONWGS84 + LATWGS84 + TIMESTRR + TAXGWRB + SOURCEDB + BDRICM
## check if there are still some duplicates
profs.f@site$SOURCEID[duplicated(profs.f@site$SOURCEID)]
## spatial duplicates:
sp <- profs.f@site
coordinates(sp) <- ~ LONWGS84 + LATWGS84
sp::zerodist(sp)
## 18 duplicates!

# prepare a SITE and hor table table:
wsp4 <- list(sites=profs.f@site, horizons=profs.f@horizons)
wsp4$sites$TAXGWRB <- as.character(wsp4$sites$TAXGWRB)
str(wsp4)
lapply(as.list(wsp4$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp4$horizons), function(x){sum(!is.na(x))})
save(wsp4, file="../wsp4.rda", compress="xz")

## Visualize all points:
library(plotKML)
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(HRSGDB.ll, shape = shape, labels = PROF_ID, kmz = TRUE) 
