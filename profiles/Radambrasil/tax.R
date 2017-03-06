## BR classification clean tables

## Download the database:
download.file("http://www.esalq.usp.br/gerd/BrazilSoilDB_08VI05.xls", destfile="BrazilSoilDB_08VI05.xls")

library(aqp)
library(plyr)
library(GSIF)
library(maptools)
library(gdata)
#perl <- gdata:::findPerl("perl")
perl = "C:/Perl64/bin/perl.exe"
tmp <- read.xls("BrazilSoilDB_08VI05.xls", perl=perl, sheet=2)
## takes 2-3 mins to read!
profs <- tmp

## WRB classes added by Alessandro Rosa / Eliana de Souza:
#taxa <- read.csv("BR_taxa_correlation.csv") 
taxa <- read.xls("BR_taxa_correlation.xls", perl=perl, sheet=7) 

taxa.BR1 <- taxa[taxa$year_source==1990,c("year_source","taxon_source","taxon_destination")]
taxa.BR1$taxon_source_BR <- taxa.BR1$taxon_source
taxa.BR2 <- taxa[taxa$year_destination==2007,c("year_source","taxon_source","taxon_destination","year_destination")]
taxa.BR2$taxon_destination_WRB <- taxa.BR2$taxon_destination
taxa.BR <- merge(y=taxa.BR2[,c("taxon_source","taxon_destination_WRB")], x=taxa.BR1[,c("taxon_destination","taxon_source_BR")], by.y="taxon_source", by.x="taxon_destination", all.x=TRUE)
write.csv(taxa.BR, file="taxa_BR.csv")

## Match taxa:
levels(profs$SoilClass)
profs$taxon_source_BR1990 <- as.factor(gsub("_", " ", profs$SoilClass))
levs <- levels(taxa$taxon_source_BR1990)
levs.s <- nchar(levs)
profs$taxon_source_BR1990_ <- NA
for(j in 1:length(levs)){
  tax.t <- levs[order(levs.s, decreasing=FALSE)][j]
  if(tax.t %in% c("Latossolo","Litossolo","Podzol","Regossolo","Cambissolo","Planossolo","Vertissolo","Solonchak","Rendzina")){
    sel <- grep(tax.t, paste(profs$taxon_source_BR1990), ignore.case=TRUE, fixed=FALSE)
  } else {
    sel <- agrep(tax.t, paste(profs$taxon_source_BR1990), ignore.case=TRUE, fixed=FALSE, max.distance=.2)
  } 
  if(length(sel)>0){
    profs[sel,"taxon_source_BR1990_"] = tax.t
  }
}
taxa.levs <- levels(as.factor(paste(profs$taxon_source_BR1990, profs$taxon_source_BR1990_, sep=" : ")))
 write.csv(taxa.levs, file="tax_levs.csv")

