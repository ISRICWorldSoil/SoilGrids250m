## Prepare SoilGrids for import into GEE
## by tom.hengl@isric.org and milan.kili11@gmail.com

library(rgdal)
library(snowfall)
var.tbl = read.csv("list_gge_names.csv")
str(var.tbl)

bind_rename <- function(var, gee.name, bands){
  in.files = list.files(path="/data/GEOG", pattern=glob2rx(paste0(var, "*.tif$")), full.names = TRUE)
  if(length(in.files)>1 & bands>1){
    if(!file.exists(paste0(gee.name, '.tif'))){
      out.tmp <- tempfile(fileext = ".txt")
      cat(in.files, sep="\n", file=out.tmp)
      vrt = paste0(var, ".vrt")
      system(paste0('gdalbuildvrt -separate ', vrt, ' -input_file_list ', out.tmp))
      system(paste0('gdal_translate --config GDAL_CACHEMAX 99000 ', vrt,' ', gee.name, '.tif -co \"COMPRESS=DEFLATE\" -co BIGTIFF=YES'))
    }
  } else {
    system(paste0('cp ', in.files, ' /data/GEOG/gee/', gee.name, '.tif'))
  }
}

sfInit(parallel=TRUE, cpus=5)
sfExport("var.tbl", "bind_rename")
sfLibrary(rgdal)
out <- sfClusterApplyLB(1:nrow(var.tbl), function(i){bind_rename(var=paste(var.tbl[i,"SoilGrids_code"]), gee.name=paste(var.tbl[i,"Asset_name"]), bands=var.tbl[i,"Bands"])})
sfStop()

## band names:
x = paste(sapply(basename(list.files(path="/data/GEOG", pattern=glob2rx("TAXNWRB_*_*_*.tif$"), full.names = TRUE)), function(i){ gsub("\\.", "_", strsplit(strsplit(i, "TAXNWRB_")[[1]][2], "_250m")[[1]][1])}), collapse=",")
cat(x, sep="\n", file="WRBN_list_of_bands.txt")
x = paste(sapply(basename(list.files(path="/data/GEOG", pattern=glob2rx("TAXOUSDA_*_*_*.tif$"), full.names = TRUE)), function(i){ gsub("\\.", "_", strsplit(strsplit(i, "TAXOUSDA_")[[1]][2], "_250m")[[1]][1])}), collapse=",")
cat(x, sep="\n", file="TAXO_list_of_bands.txt")


## Install earthengine command line (https://developers.google.com/earth-engine/command_line)
## Generate code for each asset (name of bands can be long strings)

# 1.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Depth_to_bedrock_up_to_200cm --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Depth_to_bedrock_up_to_200cm.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 2.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Occurrence_R_horizon --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Occurrence_R_horizon.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 3.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Absolute_depth_to_bedrock --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Absolute_depth_to_bedrock.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 4.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_bulk_density --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Soil_bulk_density.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 5.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_CEC --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Soil_CEC.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 6.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Clay_content --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Clay_content.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 7.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Coarse_fragments --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Coarse_fragments.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 8.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_organic_carbon_stock --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Soil_organic_carbon_stock.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sd_000_005,sd_005_015,sd_015_030,sd_030_060,sd_060_100,sd_100_200")

# 9.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_organic_carbon_density --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Soil_organic_carbon_density.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 10.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_organic_carbon_content --nodata_value=-32768 --pyramiding_policy=mean gs://0soilgrids250m/Soil_organic_carbon_content.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 11.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_pH_H2O --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Soil_pH_H2O.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 12.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_pH_KCl --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Soil_pH_KCl.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 13.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Sand_content --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Sand_content.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 14.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Silt_content --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Silt_content.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 15.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_texture_class --nodata_value=255 --pyramiding_policy=sample gs://0soilgrids250m/Soil_texture_class.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands sl_000,sl_005,sl_015,sl_030,sl_060,sl_100,sl_200")

# 16.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/WRB_soil_type_most_probable --nodata_value=255 --pyramiding_policy=sample gs://0soilgrids250m/WRB_soil_type_most_probable.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 17.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/WRB_soil_type_probabilities --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/WRB_soil_type_probabilities.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands Acric_Ferralsols,Acric_Plinthosols,Albic_Arenosols,Albic_Luvisols,Alic_Nitisols,Aluandic_Andosols,Aric_Regosols,Calcaric_Regosols,Calcic_Chernozems,Calcic_Gleysols,Calcic_Gypsisols,Calcic_Histosols,Calcic_Kastanozems,Calcic_Luvisols,Calcic_Solonetz,Calcic_Vertisols,Cryic_Histosols,Cutanic_Alisols,Endogleyic_Cambisols,Endogleyic_Planosols,Ferralic_Arenosols,Ferralic_Cambisols,Fibric_Histosols,Gleyic_Luvisols,Gleyic_Podzols,Gleyic_Solonetz,Gypsic_Solonchaks,Haplic_Acrisols,Haplic_Acrisols__Alumic_,Haplic_Acrisols__Ferric_,Haplic_Acrisols__Humic_,Haplic_Albeluvisols,Haplic_Alisols,Haplic_Andosols,Haplic_Arenosols,Haplic_Arenosols__Calcaric_,Haplic_Calcisols,Haplic_Calcisols__Sodic_,Haplic_Cambisols,Haplic_Cambisols__Calcaric_,Haplic_Cambisols__Chromic_,Haplic_Cambisols__Dystric_,Haplic_Cambisols__Eutric_,Haplic_Cambisols__Humic_,Haplic_Cambisols__Sodic_,Haplic_Chernozems,Haplic_Cryosols,Haplic_Ferralsols,Haplic_Ferralsols__Rhodic_,Haplic_Ferralsols__Xanthic_,Haplic_Fluvisols,Haplic_Fluvisols__Arenic_,Haplic_Fluvisols__Calcaric_,Haplic_Fluvisols__Dystric_,Haplic_Fluvisols__Eutric_,Haplic_Gleysols,Haplic_Gleysols__Dystric_,Haplic_Gleysols__Eutric_,Haplic_Gypsisols,Haplic_Kastanozems,Haplic_Leptosols,Haplic_Leptosols__Eutric_,Haplic_Lixisols,Haplic_Lixisols__Chromic_,Haplic_Lixisols__Ferric_,Haplic_Luvisols,Haplic_Luvisols__Chromic_,Haplic_Luvisols__Ferric_,Haplic_Nitisols__Rhodic_,Haplic_Phaeozems,Haplic_Planosols__Dystric_,Haplic_Planosols__Eutric_,Haplic_Podzols,Haplic_Regosols__Dystric_,Haplic_Regosols__Eutric_,Haplic_Regosols__Sodic_,Haplic_Solonchaks,Haplic_Solonchaks__Sodic_,Haplic_Solonetz,Haplic_Umbrisols,Haplic_Vertisols,Haplic_Vertisols__Eutric_,Hemic_Histosols,Histic_Albeluvisols,Hypoluvic_Arenosols,Leptic_Cambisols,Leptic_Luvisols,Leptic_Phaeozems,Leptic_Regosols,Leptic_Umbrisols,Lithic_Leptosols,Lixic_Plinthosols,Luvic_Calcisols,Luvic_Chernozems,Luvic_Phaeozems,Luvic_Planosols,Luvic_Stagnosols,Mollic_Gleysols,Mollic_Leptosols,Mollic_Solonetz,Mollic_Vertisols,Petric_Calcisols,Petric_Durisols,Plinthic_Acrisols,Protic_Arenosols,Rendzic_Leptosols,Sapric_Histosols,Solodic_Planosols,Stagnic_Luvisols,Turbic_Cryosols,Umbric_Albeluvisols,Umbric_Ferralsols,Umbric_Gleysols,Vertic_Cambisols,Vertic_Luvisols,Vetic_Acrisols,Vitric_Andosols,Vitric_Cryosols")

# 18.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/USDA_soil_type_most_probable --nodata_value=255 --pyramiding_policy=sample gs://0soilgrids250m/USDA_soil_type_most_probable.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 19.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/USDA_soil_type_probabilities --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/USDA_soil_type_probabilities.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands Albolls,Aqualfs,Aquands,Aquents,Aquepts,Aquerts,Aquods,Aquolls,Aquox,Aquults,Arents,Argids,Borolls,Calcids,Cambids,Cryalfs,Cryands,Cryepts,Cryids,Cryods,Cryolls,Durids,Fibrists,Fluvents,Folists,Gelands,Gelods,Gypsids,Hemists,Histels,Humods,Humults,Ochrepts,Orthels,Orthents,Orthods,Perox,Psamments,Rendolls,Salids,Saprists,Torrands,Torrerts,Torrox,Turbels,Udalfs,Udands,Udepts,Uderts,Udolls,Udox,Udults,Ustalfs,Ustands,Ustepts,Usterts,Ustolls,Ustox,Ustults,Vitrands,Xeralfs,Xerands,Xerepts,Xererts,Xerolls,Xerults")

# 20.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Available_soil_water_capacity --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Available_soil_water_capacity.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00 --bands h1_000,h1_005,h1_015,h1_030,h1_060,h1_100,h1_200,h2_000,h2_005,h2_015,h2_030,h2_060,h2_100,h2_200,h3_000,h3_005,h3_015,h3_030,h3_060,h3_100,h3_200,tS_000,tS_005,tS_015,tS_030,tS_060,tS_100,tS_200,wp_000,wp_005,wp_015,wp_030,wp_060,wp_100,wp_200")

# 21.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Organic_soil_type_probability --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Organic_soil_type_probability.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 22.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_sodic_grade --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Soil_sodic_grade.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")

# 23.
system("earthengine upload image --asset_id=users/tomhengl/SoilGrids250m/Soil_subsoil_acidity_grade --nodata_value=255 --pyramiding_policy=mean gs://0soilgrids250m/Soil_subsoil_acidity_grade.tif --property '(string)name=SoilGrids250m' --time_start 1970-01-01T12:00:00 --time_end 2014-01-01T12:00:00")


#   usage: earthengine upload image [-h] [--wait [WAIT]] [--asset_id ASSET_ID]
# [--last_band_alpha]
# [--nodata_value NODATA_VALUE]
# [--pyramiding_policy PYRAMIDING_POLICY]
# [--bands BANDS] [--crs CRS]
# [--property PROPERTY]
# [--time_start TIME_START]
# [--time_end TIME_END]
# src_files [src_files ...]
# 
# Uploads an image from Cloud Storage to Earth Engine. See docs for "asset set"
# for additional details on how to specify asset metadata properties.
