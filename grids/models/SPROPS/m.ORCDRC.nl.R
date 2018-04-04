# This script aims to run the prediction for Organic Carbon restrictied to the Netherlands
# The computation should be light enough to run in a desktop computer.
# It allows a first contact with the code and the necessary data sources.
# This script was built from the main m.SPROPS.R script.
#
# Author: Luís de Sousa
###############################################################################

# These lines should be run in a root session
list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "ranger", "xgboost", "mda", "psych", "stringr", "caret", "plotKML", "maptools", "maps", "stringr", "R.utils", "grDevices", "GSIF")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Location of the data folder - note that this script was originally intended to to find code and data all in the same path
# Usually, data files (e.g. model outputs) are not in the code repository because they are too large
data_path <- "~/SoilGridsData/"

# Location of the local code repository
repo_path <- "~/git/SoilGrids250m/"

# On production, this must be adapted to either the data or the code folder
setwd("~/git/SoilGrids250m/grids/models/SPROPS")

# Load required libraries
library(plyr)
library(stringr)
library(sp)
library(rgdal)
library(devtools)
library(xgboost) ## xgboost_0.6-4
library(ranger) ## ranger_0.6.7
library(caret)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(snowfall)
library(utils)
library(plotKML)
library(R.utils)
library(GSIF)
library(parallel)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "gdalwarp"
gdalbuildvrt = "gdalbuildvrt"
system("gdal-config --version")
source(paste(repo_path, "grids/models/wrapper.predict_cs.R", sep=""))
source(paste(repo_path, "grids/models/saveRDS_functions.R", sep=""))
source(paste(repo_path, "grids/models/mosaick_functions_ll.R", sep=""))
source(paste(repo_path, "grids/models/extract_tiled.R", sep=""))
## metadata:
metasd <- read.csv(paste(repo_path, 'grids/GEOG/META_GEOTIFF_1B.csv', sep=""), stringsAsFactors = FALSE)
sel.metasd = names(metasd)[-sapply(c("FileName","VARIABLE_NAME"), function(x){grep(x, names(metasd))})]
## covariates:
des <- read.csv(paste(data_path, "models/SoilGrids250m_COVS250m.csv", sep=""))
mask_value <- as.list(des$MASK_VALUE)
names(mask_value) = des$WORLDGRIDS_CODE

## points:
load(paste(data_path, "profs/SPROPS/SPROPS.pnts.rda", sep="")) ## spatial locations only
## 173,806 points
load(paste(data_path, "profs/SPROPS/all.pnts.rda", sep=""))

## spatia overlay (20 mins):
## This must be made for the entire points dataset, otherwise some properties may disappear.
if(!file.exists("ovA.rds"))
{
	tile.pol = rgdal::readOGR(paste(data_path, "models/tiles_ll_100km.shp", sep=""), "tiles_ll_100km")
	ov <- extract.tiled(x=SPROPS.pnts, tile.pol=tile.pol, path=paste(data_path, "tt/SoilGrids250m/predicted250m", sep=""), ID="ID", cpus=10)
	ovA <- join(all.pnts, ov, type="left", by="LOC_ID")
	## Save to avoid computing a second time
	saveRDS(ovA, "ovA.rds")
} 
else # If the file is already present then ovA was already computed
{ 
	loadRDS("ovA.rds")
}

