
plotCovsSoilGrids250m <- function(s, path, ATTRIBUTE_TITLE, DESCRIPTION, EVI_range, ES_range, TD_range, TN_range, NBR4_range, NBR7_range, SNW_range){
  sfInit(parallel=TRUE, cpus=40)
  sfLibrary(plotKML)
  sfClusterEval( plotKML.env(convert="convert", show.env=FALSE) )
  sfLibrary(XML)
  sfLibrary(rgdal)
  sfLibrary(raster)
  sfExport("s", "path", "ATTRIBUTE_TITLE", "DESCRIPTION", "EVI_range", "ES_range", "TD_range", "TN_range", "NBR4_range", "NBR7_range", "SNW_range", ".writePNG")
  x <- sfLapply(1:length(names(s)), .writePNG, s=s, path=path, ATTRIBUTE_TITLE=ATTRIBUTE_TITLE, DESCRIPTION=DESCRIPTION, EVI_range=EVI_range, ES_range=ES_range, TD_range=ES_range, TN_range=TN_range, NBR4_range=NBR4_range, NBR7_range=NBR7_range, SNW_range=SNW_range) 
  sfStop()
}

.writePNG <- function(i, s, path, ATTRIBUTE_TITLE, DESCRIPTION, EVI_range, ES_range, TD_range, TN_range, NBR4_range, NBR7_range, SNW_range){
  nm <- strsplit(names(s)[i], "_")[[1]][1]
  r <- raster(s, i)
  names(r)="colour"
  file.name <- paste0(path, names(s)[i], ".kml")
  if(length(grep(nm, pattern="DEMMRG5"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern="SLPMRG5"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=R_pal[["heat_colors"]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern="CRVMRG5"))>0|length(grep(nm, pattern="VBFMRG5"))>0|length(grep(nm, pattern="VDPMRG5"))>0|length(grep(nm, pattern="TWIMRG5"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern="DVMMRG5"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[2]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern="DVMMRG5"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=R_pal[["tex_pal"]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern="EX"))>0){
    if(length(grep(nm, pattern="EX1"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 03, 01)} 
    if(length(grep(nm, pattern="EX2"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="EX3"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="EX4"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 09, 28)}
    if(length(grep(nm, pattern="EX5"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="EX6"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_YELLOW_GREEN"]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=EVI_range, plot.legend=FALSE)
  }
  if(length(grep(nm, pattern="ES"))>0){
    if(length(grep(nm, pattern="ES1"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 03, 01)} 
    if(length(grep(nm, pattern="ES2"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="ES3"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="ES4"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 09, 28)}
    if(length(grep(nm, pattern="ES5"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="ES6"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=ES_range, plot.legend=FALSE)
  }
  if(length(grep(nm, pattern=glob2rx("I??MOD4")))>0){
    if(length(grep(nm, pattern="I01"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 02, 01)} 
    if(length(grep(nm, pattern="I02"))>0){ tbegin=ISOdate(2015, 02, 01); tend=ISOdate(2015, 03, 01)}
    if(length(grep(nm, pattern="I03"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 04, 01)}
    if(length(grep(nm, pattern="I04"))>0){ tbegin=ISOdate(2015, 04, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="I05"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 06, 01)}
    if(length(grep(nm, pattern="I06"))>0){ tbegin=ISOdate(2015, 06, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="I07"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 08, 01)}
    if(length(grep(nm, pattern="I08"))>0){ tbegin=ISOdate(2015, 08, 01); tend=ISOdate(2015, 09, 01)}
    if(length(grep(nm, pattern="I09"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 10, 01)}
    if(length(grep(nm, pattern="I10"))>0){ tbegin=ISOdate(2015, 10, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="I11"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2015, 12, 01)}
    if(length(grep(nm, pattern="I12"))>0){ tbegin=ISOdate(2015, 12, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=rev(SAGA_pal[["SG_COLORS_YELLOW_RED"]]), TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=NBR4_range, plot.legend=FALSE)
  }
  if(length(grep(nm, pattern=glob2rx("M??MOD4")))>0){
    if(length(grep(nm, pattern="M01"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 02, 01)} 
    if(length(grep(nm, pattern="M02"))>0){ tbegin=ISOdate(2015, 02, 01); tend=ISOdate(2015, 03, 01)}
    if(length(grep(nm, pattern="M03"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 04, 01)}
    if(length(grep(nm, pattern="M04"))>0){ tbegin=ISOdate(2015, 04, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="M05"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 06, 01)}
    if(length(grep(nm, pattern="M06"))>0){ tbegin=ISOdate(2015, 06, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="M07"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 08, 01)}
    if(length(grep(nm, pattern="M08"))>0){ tbegin=ISOdate(2015, 08, 01); tend=ISOdate(2015, 09, 01)}
    if(length(grep(nm, pattern="M09"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 10, 01)}
    if(length(grep(nm, pattern="M10"))>0){ tbegin=ISOdate(2015, 10, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="M11"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2015, 12, 01)}
    if(length(grep(nm, pattern="M12"))>0){ tbegin=ISOdate(2015, 12, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=rev(SAGA_pal[["SG_COLORS_YELLOW_BLUE"]]), TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=NBR7_range, plot.legend=FALSE)
  }
  if(length(grep(nm, pattern=glob2rx("T??MOD4")))>0){
    if(length(grep(nm, pattern="T01"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 02, 01)} 
    if(length(grep(nm, pattern="T02"))>0){ tbegin=ISOdate(2015, 02, 01); tend=ISOdate(2015, 03, 01)}
    if(length(grep(nm, pattern="T03"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 04, 01)}
    if(length(grep(nm, pattern="T04"))>0){ tbegin=ISOdate(2015, 04, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="T05"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 06, 01)}
    if(length(grep(nm, pattern="T06"))>0){ tbegin=ISOdate(2015, 06, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="T07"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 08, 01)}
    if(length(grep(nm, pattern="T08"))>0){ tbegin=ISOdate(2015, 08, 01); tend=ISOdate(2015, 09, 01)}
    if(length(grep(nm, pattern="T09"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 10, 01)}
    if(length(grep(nm, pattern="T10"))>0){ tbegin=ISOdate(2015, 10, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="T11"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2015, 12, 01)}
    if(length(grep(nm, pattern="T12"))>0){ tbegin=ISOdate(2015, 12, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=TD_range)
  }
  if(length(grep(nm, pattern=glob2rx("N??MOD4")))>0){
    if(length(grep(nm, pattern="N01"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 02, 01)} 
    if(length(grep(nm, pattern="N02"))>0){ tbegin=ISOdate(2015, 02, 01); tend=ISOdate(2015, 03, 01)}
    if(length(grep(nm, pattern="N03"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 04, 01)}
    if(length(grep(nm, pattern="N04"))>0){ tbegin=ISOdate(2015, 04, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="N05"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 06, 01)}
    if(length(grep(nm, pattern="N06"))>0){ tbegin=ISOdate(2015, 06, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="N07"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 08, 01)}
    if(length(grep(nm, pattern="N08"))>0){ tbegin=ISOdate(2015, 08, 01); tend=ISOdate(2015, 09, 01)}
    if(length(grep(nm, pattern="N09"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 10, 01)}
    if(length(grep(nm, pattern="N10"))>0){ tbegin=ISOdate(2015, 10, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="N11"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2015, 12, 01)}
    if(length(grep(nm, pattern="N12"))>0){ tbegin=ISOdate(2015, 12, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=TN_range)
  }
  if(length(grep(nm, pattern="TMDMOD3"))>0|length(grep(nm, pattern="TSDMOD3"))>0|length(grep(nm, pattern="TMNMOD3"))>0|length(grep(nm, pattern="TSNMOD3"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], raster_name=paste0("./", names(s)[i], ".png"))
  }
  if(length(grep(nm, pattern=glob2rx("SN??MOD4")))>0){
    if(length(grep(nm, pattern="SN1"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 03, 01)} 
    if(length(grep(nm, pattern="SN2"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="SN3"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="SN4"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 09, 28)}
    if(length(grep(nm, pattern="SN5"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="SN6"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_YELLOW_BLUE"]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=SNW_range, plot.legend=FALSE)
  }
  if(length(grep(nm, pattern="GLC"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_YELLOW_GREEN"]], raster_name=paste0("./", names(s)[i], ".png"), z.lim=c(0,100), plot.legend=FALSE)
  }
  if(length(grep(nm, pattern=glob2rx("P??MRG3")))>0){
    if(length(grep(nm, pattern="P01"))>0){ tbegin=ISOdate(2015, 01, 01); tend=ISOdate(2015, 02, 01)} 
    if(length(grep(nm, pattern="P02"))>0){ tbegin=ISOdate(2015, 02, 01); tend=ISOdate(2015, 03, 01)}
    if(length(grep(nm, pattern="P03"))>0){ tbegin=ISOdate(2015, 03, 01); tend=ISOdate(2015, 04, 01)}
    if(length(grep(nm, pattern="P04"))>0){ tbegin=ISOdate(2015, 04, 01); tend=ISOdate(2015, 05, 01)}
    if(length(grep(nm, pattern="P05"))>0){ tbegin=ISOdate(2015, 05, 01); tend=ISOdate(2015, 06, 01)}
    if(length(grep(nm, pattern="P06"))>0){ tbegin=ISOdate(2015, 06, 01); tend=ISOdate(2015, 07, 01)}
    if(length(grep(nm, pattern="P07"))>0){ tbegin=ISOdate(2015, 07, 01); tend=ISOdate(2015, 08, 01)}
    if(length(grep(nm, pattern="P08"))>0){ tbegin=ISOdate(2015, 08, 01); tend=ISOdate(2015, 09, 01)}
    if(length(grep(nm, pattern="P09"))>0){ tbegin=ISOdate(2015, 09, 01); tend=ISOdate(2015, 10, 01)}
    if(length(grep(nm, pattern="P10"))>0){ tbegin=ISOdate(2015, 10, 01); tend=ISOdate(2015, 11, 01)}
    if(length(grep(nm, pattern="P11"))>0){ tbegin=ISOdate(2015, 11, 01); tend=ISOdate(2015, 12, 01)}
    if(length(grep(nm, pattern="P12"))>0){ tbegin=ISOdate(2015, 12, 01); tend=ISOdate(2016, 01, 01)}
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=log1p(colour), layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_WHITE_BLUE"]], TimeSpan.begin=tbegin, TimeSpan.end=tend, raster_name=paste0("./", names(s)[i], ".png"), z.lim=log1p(c(0,450)*0.05), plot.legend=FALSE)
  }
  if(length(grep(nm, pattern="USG"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], raster_name=paste0("./", names(s)[i], ".png"), z.lim=c(0,100), plot.legend=FALSE)
  }
  if(length(grep(nm, pattern="GTD"))>0){
    kml(r, folder.name=names(s)[i], file.name=file.name, colour=colour, layer.name=nm, subfolder.name=ATTRIBUTE_TITLE[i], colour_scale=SAGA_pal[[1]], raster_name=paste0("./", names(s)[i], ".png"), plot.legend=FALSE)
  }
}
