
############## LIBRARIES ####################### #######################
library(devtools)
library(Rcpp)
library(ncdf4)
library(terra)
library(sf)
library(lubridate)
library(mesoclim)
library(mesoclimAddTrees)

##############  DIRECTORIES ####################### #######################
# Root directory relative to these data inputs
dir_od<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data"
dir_root<-"/Users/jonathanmosedale/Data"
dir_in<-file.path(dir_root,'mesoclim_inputs')
list.files(dir_in)

dir_ukcp<-dir_in
dir.exists(dir_ukcp)

dir_sst<-dir_in
dir.exists(dir_sst)

##############  FILEPATHS ####################### #######################
coast_file<-file.path(dir_in,'CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
file.exists(coast_file)
ukdtm_file<-file.path(dir_in,"uk_dtm.tif") # 50m dtm  raster
file.exists(ukdtm_file)
ukcpdtm_file<-file.path(dir_root,"mesoclim_inputs","orog_land-rcm_uk_12km_osgb.nc")
file.exists(ukcpdtm_file)



