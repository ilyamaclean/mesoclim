#####################################################################
# To create bias correction models for the whole of the UK for correction of future UKCP modeled data
# 1.A Prepare observed data for chosen years from ERA5 and convert to daily values.
# 1.B. Merge with HadUK data sources for same years.
# 2. Prepare modelled UKCP18rcm data for the same years
# 3. Calculate and save bias correction models for each variable
#####################################################################

library(devtools)
load_all()

library(terra)
library(sf)
library(mesoclim)
library(lubridate)
library(mesoclimAddTrees)
library(ncdf4)

#####################################################################
# Setup parameters etc
#####################################################################
dir_root<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data"

# Dir source of era5 hourly data and output for daily data
dir_era5<-file.path(dir_root,'era5')# raw 2018, 20,21,22, 24
dir_era5dly<-file.path(dir_root,'era5')# pre-prep 2010-20

# Dirs for haduk data
dir_haduk1km<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/haduk_1km_monthly"
dir_haduk12km<-"/Users/jonathanmosedale/Data/haduk_12km"
#dir_hadukcp<-file.path(dir_root,'haduk','haduk_vars')

# UKCP data
dir_ukcp<- '/Users/jonathanmosedale/Data/mesoclim_inputs' # 2 decades of data 2010-30
#dir_ukcpsst<- '/Users/jonathanmosedale/Data/mesoclim_inputs' # 2019-31

# dir for writinng models and data for corretion
dir_bmodels<-file.path(dir_root,'bias_correct_models')
dir_bcdata<-file.path(dir_root,'bc_data')

# DTM and land/sea mask files
era5lsm_file<-file.path(dir_era5,"lsm.tif")
era5geop_file<-file.path(dir_era5,"geopotential.tif")
ukcpdtm_file<-file.path(dir_ukcp,'orog_land-rcm_uk_12km_osgb.nc')

# Proportion of land/sea in ukcp 12km grid
lsm12km.r<-rast("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_12km.tif")
plot(lsm12km.r)


# ERA5 Landsea mask - convert 0 to NA - proportion of cell that is land
lsm_era5<-rast(era5lsm_file,subds='lsm')

# ERA5 Elevation (geometric height) from geopotential
#geop<-mask(rast(era5geop_file),ifel(lsm_era5==0,NA,lsm))
geop<-rast(era5geop_file)
elev_era5<-mesoclimAddTrees::elev_from_geop(geop)
plot(elev_era5)

# UKCP18rcm elevation
ukcpdtm<-rast(ukcpdtm_file)

#### Years for which historic observed and modelled data will be compared
years<-c(2020,2021,2022) # respectively among the coldest 2018, wettest, driest2022, hottest 2012 of recent years


#####################################################################
# 1.B. OBSERVED DATA - Get Haduk data - whole years resampled to UKCOP18rcm
# ONLY covers UK - not rest of Ireland
# Remove time from time variable (just date)
# USES FUNCTION: addtrees_hadukdata
#####################################################################
write_yr_tif<-TRUE # TRUE if yearly tifs required

#download_haduk(dir_haduk,  as.POSIXlt(paste0(years[1],"/01/01")),as.POSIXlt(paste0(years[length(years)],"/12/31")),
#               c("tasmax","tasmin"),freq = "day", cedausr, cedapwd, access_token)

for(yr in years){
  print(yr)
  startdate<-as.POSIXlt(paste0(yr,"/01/01"))
  enddate<-as.POSIXlt(paste0(yr,"/12/31"))

  # Get daily precipitation, tmax and tmin

  r<-import_hadukdata(startdate, enddate, NA, filepath=dir_haduk1km, haduk_res="1km", var = "rainfall")
  if(write_yr_tif) writeRaster(r,file.path(dir_haduk1km,paste0('rain_',yr,'.tif')),overwrite=TRUE )

  r<-import_hadukdata(startdate, enddate, NA, filepath=dir_haduk1km, haduk_res="1km", var = "tasmax")
  if(write_yr_tif) writeRaster(r,file.path(dir_haduk1km,paste0('tasmax_',yr,'.tif')) )

  r<-import_hadukdata(startdate, enddate, NA, filepath=dir_haduk1km, haduk_res="1km", var = "tasmin")
  if(write_yr_tif) writeRaster(r,file.path(dir_haduk1km,paste0('tasmin_',yr,'.tif')) )
}

######## Crop and aggregate to UKCP 12km grid then interpolate any coastal cells without values ########

rain.r<-rast(file.path(dir_haduk1km,paste0('rain_',years,'.tif')))
rain.r<-crop(extend(rain.r,ext(lsm12km.r)),ext(lsm12km.r))
rain12km<-terra::aggregate(rain.r,12,fun="mean",na.rm=TRUE)
rain12km<-mask( .spatinterp(rain12km),lsm12km.r, maskvalue=0 )
writeRaster(rain12km,file.path(dir_haduk12km,'rain_12km_3yrs.tif'),overwrite=TRUE)

tasmax.r<-rast(file.path(dir_haduk1km,paste0('tasmax_',years,'.tif')))
tasmax.r<-crop(extend(tasmax.r,ext(lsm12km.r)),ext(lsm12km.r))
tmax12km<-terra::aggregate(tasmax.r,12,fun="mean",na.rm=TRUE)
tmax12km<-mask( .spatinterp(tmax12km),lsm12km.r, maskvalue=0 )
writeRaster(tmax12km,file.path(dir_haduk12km,'tmax_12km_3yrs.tif'),overwrite=TRUE)

tasmin.r<-rast(file.path(dir_haduk1km,paste0('tasmin_',years,'.tif')))
tasmin.r<-crop(extend(tasmin.r,ext(lsm12km.r)),ext(lsm12km.r))
tmin12km<-terra::aggregate(tasmin.r,12,fun="mean",na.rm=TRUE)
tmin12km<-mask( .spatinterp(tmin12km),lsm12km.r, maskvalue=0 )
writeRaster(tmin12km,file.path(dir_haduk12km,'tmin_12km_3yrs.tif'),overwrite=TRUE)

# Observed list of all climate variables
obs_list<-list()
obs_list$tmax<-rast()
obs_list$tmin<-rast()
obs_list$prec<-rast()

obs_list$prec<-rain12km
obs_list$tmax<-tmax12km
obs_list$tmin<-tmin12km

#####################################################################
# Get Haduk data to observation list if yearly files already available
#####################################################################
#obs_list<-list()
#years<-seq(startyr,endyr,1)
#obs_list$tmax<-rast(file.path(dir_hadukcp,paste0('tasmax_',years,'.tif')))
#obs_list$tmin<-rast(file.path(dir_hadukcp,paste0('tasmin_',years,'.tif')))
#obs_list$prec<-rast(file.path(dir_hadukcp,paste0('rain_',years,'.tif')))


#####################################################################
# 1.A. OBSERVED DATA - Get ERA5 hourly variables - convert to daily and resample to UKCP
# If several years chosen then will create list of spatrasters for each variable for the entire time period
#####################################################################
# Variables to get from ERA5 data for bias correction
vars<-c('psl','huss','hurs','u10','v10','swrad','difrad','lwrad')

###### Prepare and calculate daily values for each file using .era5todaily and writes daily nc file
era5files<-paste0("era5_surface_ukeire_",years,".nc")
if(all(file.exists(file.path(dir_era5,era5files)))!=TRUE) stop(paste("ERA5 hourly file not found in directory",dir_era5,"!!!!"))

era5files<-era5files[order(era5files)]
era5day_list<-list()

t0<-now()
for (f in era5files){
  print(f)
  daily<-era5todaily(filein=file.path(dir_era5,f),pathout=dir_era5dly,landsea=lsm_era5,elev=elev_era5,output='singlenc',vars)
  keys<-unique(c(names(daily),names(era5day_list)))
  if(f==era5files[[1]]) era5day_list<-daily else era5day_list<-mapply(c, era5day_list[keys], daily[keys], SIMPLIFY = FALSE)
  print(now()-t0)
}

### Sample plots from last saved file
for(v in names(era5day_list)){
  plot(era5day_list[[v]][[12]],main=v)
}


###### Resample daily ERA5 to ukcp resolution ##########################
for (v in vars) obs_list[[v]]<-.resample(era5day_list[[v]],lsm12km.r, msk=TRUE,method="bilinear")

###### Calculate derived variables ##########################
# Re-calculate rel humidity? - using d2m, surface pressure and temperature
#tc<-as.array((obs_list$tasmax+obs_list$tasmin)/2)
#ea<-.satvap(as.array(obs_list$d2m)-273.15)
#es<-.satvap(tc)
#rh<-(ea/es)*100

# OR Recalculate rel humidity from spec humidity and new pk/tc
#n<-dim(obs_list$huss)[3]
#pk<- as.array(obs_list$psl / ((293-0.0065*ukcpdtm)/293)^5.26)
#rh<-suppressWarnings(converthumidity(.is(obs_list$huss), intype = "specific", outtype = "relative",tc,pk ))

# Very high rh values produced so impose limits
#rh[rh>100]<-100
#rh[rh<20]<-20
#obs_list$relhum<-.rast(rh,ukcpdtm)

# Calculate wind speed - convert to 2m above ground - KEEP AT 10m
wsp<-sqrt(obs_list$u10^2+obs_list$v10^2)
wsp[wsp<0.001]<-0.001
##wsp<-.windhgt(wsp, zi=10, zo=2)
obs_list[['windspeed']]<-wsp

# Calculate winddir
obs_list[['winddir']]<-(atan2(obs_list$u10,obs_list$v10)*(180/pi)+180)%%360
obs_list<-obs_list[-which(names(obs_list) %in% c('v10','u10'))]

# Change name psl to pres???
write_climdata(obs_list,filepath=file.path(dir_bcdata,'obs_list_2020_2022.Rds'),overwrite=TRUE)



#####################################################################
# 2. Get UKCP (model) data
#####################################################################
dtmc<-crop(ukcpdtm,lsm12km.r)
plot(dtmc)
collection<-'land-rcm'
domain<-'uk'
member<-'01'
rcp<-'rcp85'
startdate<-as.POSIXlt(paste0(years[1],'/01/01'))
enddate<-as.POSIXlt(paste0(years[length(years)],'/12/31'))

climdata<-ukcp18toclimarray(dir_ukcp,dtmc,startdate,enddate,collection,domain,member,temp_hgt=1.5, wind_hgt=10,toArrays=FALSE)

#. ukcp18rcm<-checkinputs(read_climdata(file.path(dir_out,'ukcp18rcm.Rds')), tstep = "day")

write_climdata(climdata,filepath=file.path(dir_bcdata,'ukcpdata_2020_2022.Rds'),overwrite=TRUE)

#####################################################################
# 3. Calculate biascorrect models
#####################################################################
climdata<-readRDS(file.path(dir_bcdata,'ukcpdata_2020_2022.Rds'))
climdata[ which(lapply(lapply(climdata,class),"[",1)=="PackedSpatRaster") ]<-lapply(climdata[which(lapply(lapply(climdata,class),"[",1)=="PackedSpatRaster") ],unwrap)
# Define area that can be corrected - all UKCP cells that possess land cells from Terrain50 (lsmask12km.r)
msk<-ifel(lsm12km.r==0,NA,lsm12km.r)

# Get variables and correct names in obs list if required
vars<-c("relhum","pres","lwrad","swrad","tmax","tmin","windspeed" )
names(obs_list)[which(names(obs_list)=='psl')]<-"pres"
names(obs_list)[which(names(obs_list)=='hurs')]<-"relhum"
vars %in% names(obs_list)
vars %in% names(climdata)

days<-nlyr(obs_list$tmax)
days==nlyr(climdata$tmax)

model_list<-list()
# t0<-now();test1000<-biascorrect(obs_list[[v]],  mask(modelled_list[[v]],msk), mod_out = TRUE, samplenum=1461);now()-t0
t0<-now()
for (v in vars) {
  print(v)
  model_list[[v]]<-biascorrect(mask(obs_list[[v]],msk), hist_mod=mask(climdata[[v]],msk), fut_mod = NA,mod_out = TRUE, rangelims = NA, samplenum=days)
}
print(now()-t0)

# Precipitation - probably need to check this more re thresholds??
v<-"prec"
model_list[[v]]<-precipcorrect(mask(obs_list[[v]],msk), hist_mod=mask(climdata[[v]],msk), mod_out = TRUE)

saveRDS(model_list,file.path(dir_bmodels,paste0('bias_correct_models_2020_2022.Rds')))

#####################################################################
# 4. Apply biascorrect models
#####################################################################
model_list[["prec"]]<-lapply(model_list[["prec"]],unwrap)

# Apply bc models
bc_list<-list()
for (v in vars) {
  bc_list[[v]]<-biascorrect_apply(fut_mod=mask(climdata[[v]],msk), biasmods=model_list[[v]], rangelims = NA)

  dif<-mask(obs_list[[v]],msk)-bc_list[[v]]
  stats<-c(app(dif,mean),app(dif,sd),(app(mask(obs_list[[v]],msk),max)-app(mask(obs_list[[v]],msk),min))-(app(bc_list[[v]],max)-app(bc_list[[v]],min)),mean(bc_list[[v]]))
  names(stats)<-c("Mean obs-bc_model","sd obs-bc_model","obs range - bc_model range","Mean bc")
  plot(stats)
}

# Apply bc model to precipitation - bit more exploration - decide on prec_threhold =0.01
v<-"prec"
bc_list[[v]]<-precipcorrect_apply(fut_mod=mask(climdata[[v]],msk), biasmods=model_list[[v]])
dif<-mask(obs_list[[v]],msk)-bc_list[[v]]
stats<-c(app(dif,mean),app(dif,sd),(app(mask(obs_list[[v]],msk),max)-app(mask(obs_list[[v]],msk),min))-(app(bc_list[[v]],max)-app(bc_list[[v]],min)),mean(bc_list[[v]]))
names(stats)<-c("Mean obs-bc_model","sd obs-bc_model","obs range - bc_model range","Mean bc")
plot(stats)








# Explore results. - USE RASTERVIS https://oscarperpinan.github.io/rastervis/!!!
plot(unwrap(model_list$prec$mu_tot),main='Total rainfall')
plot(unwrap(model_list$prec$mu_frac ),main='Proportion of rainy days')

cf_maxs<- c(app(obs_list$tmax,max),mask(app(ukcp18rcm$tmax,max),msk))
names(cf_maxs)<-c('Yr tmx HadUK Observed','Yr tmx UKCP18 downsscaled')
plot(cf_maxs)

# SW rad
cf<- c(app(obs_list$swrad,mean),mask(app(ukcp18rcm$swrad,mean),msk))
names(cf)<-c('Observed mean','UKCP18 mean')
plot(mask(cf,msk))




