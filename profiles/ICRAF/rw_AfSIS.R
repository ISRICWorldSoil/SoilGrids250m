## AfSIS Sentinel site data from 2014



# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

load("afsp.spc.rda")
str(afsp.spc)
afsp.spc$SAMPLEID <- make.unique(as.character(afsp.spc$SOURCEID))
afsp.spc$DEPTH <- afsp.spc$UHDICM + (afsp.spc$LHDICM - afsp.spc$UHDICM)/2
SPROPS.AfSIS <- afsp.spc[afsp.spc$SOURCEDB=="Af_soilspec"|afsp.spc$SOURCEDB=="Af_LDSF",c("SOURCEID","UHDICM","LHDICM","DEPTH","BLD","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","ORCDRC","LONWGS84","LATWGS84", "SAMPLEID")]
SPROPS.AfSIS$SOURCEDB <- "AfSIS_Sentinel_Sites"
SPROPS.AfSIS <- SPROPS.AfSIS[!is.na(SPROPS.AfSIS$LONWGS84) & !is.na(SPROPS.AfSIS$LATWGS84),]
View(SPROPS.AfSIS)
## 18,055
save(SPROPS.AfSIS, file="SPROPS.AfSIS.rda")
