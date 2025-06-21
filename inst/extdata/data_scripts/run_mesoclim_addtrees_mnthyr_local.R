# Set up libraries files and directories for mac/jasmin
#  setup_file<-"inst/extdata/data_scripts/setup_mesoclim_mac.R"
#  source(setup_file)
#  devtools::load_all()
############## LIBRARIES ####################### #######################
library(devtools)
library(Rcpp)
library(ncdf4)
library(terra)
library(sf)
library(lubridate)
library(mesoclim)
library(mesoclimAddTrees)
tstart<-now()

##############  DIRECTORIES - UPDATE THESE ####################### #######################
# Root directory relative to these data inputs
# dir_od<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data"
dir_root<-"/Users/jonathanmosedale/Data"
dir_in<-file.path(dir_root,'mesoclim_inputs')

### UKCP input data directory
dir_ukcp<-dir_in
dir.exists(dir_ukcp)

### Sea temperature input data directory
dir_sst<-file.path(dir_in,"sst")
dir.exists(dir_sst)

# Directory holding bias correction models
dir_bcmodels<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/bias_correct_models"
dir.exists(dir_bcmodels)

#### Output directory
dir_out<-"/Users/jonathanmosedale/Data/mesoclim_outputs/exmoor_conifers"

##############  FILEPATHS - UPDATE THESE ####################### #######################
# Coastal mask
coast_file<-file.path(dir_in,'CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
file.exists(coast_file)

# Fine resolution DTM for UK = Terrain50
ukdtm_file<-file.path(dir_in,"uk_dtm.tif") # 50m dtm  raster
file.exists(ukdtm_file)

# UKCP 12km dtm
ukcpdtm_file<-file.path(dir_root,"mesoclim_inputs","orog_land-rcm_uk_12km_osgb.nc")
file.exists(ukcpdtm_file)

###### Parcel file and label for outputs - can be a tif in which case parcel made of whole extent #######
dir_parcels<-file.path(dir_root,'mesoclim_inputs','land_parcels')

#parcels_file<-file.path(dir_parcels,'dartmoor_conifers.shp') # usual test area in Cornwall - coast effect
#parcels_file<-file.path(dir_parcels,'killerton_parcels.shp') # low elev variation
#parcels_file<-file.path(dir_parcels,'cairngorm_parcels.shp') # inland & high elev variation
#parcels_file<-file.path(dir_parcels,'exmoor_parcels.shp')  # Large area - high coast and elev effects
parcels_file<-file.path(dir_parcels,'exmoor_conifers.shp')  # Large area - high coast and elev effects
#parcels_file<-file.path(dir_parcels,'skye_parcels.shp')  # coast and high elev effect
#parcels_file<-file.path(dir_parcels,'southdevon_parcels.shp') # LARGE AREA
# parcels_file<-file.path(system.file("extdata/dtms/dtmf_inland.tif",package="mesoclim")) # inland tif example

print(paste("Parcels input file:",parcels_file))
file.exists(parcels_file)
print(paste("Output directory:",dir_out))
dir.exists(dir_out)

############## RUN PARAMETERS - UPDATE THESE ####################### #######################

#### Label for outputs
arealabel<-"exmoor_conifers"

#### Parcel identifier field
#parcel_id<-"gid" # CEH parcels
parcel_id<-"OBJECTID_1" # Scortland Nat Forest parcels
#parcel_id<-"id"

#### UKCP options and time period
member<-"01"
startdate<-as.POSIXlt('2011/01/01',tz="UTC")
enddate<-as.POSIXlt('2020/12/31',tz="UTC")

#### Bias correction?
bias_correct<-TRUE

#### Which outputs - graphical, parcel csv, mesoclim grids
outputs<-FALSE
parcel_output<-TRUE
mesoclim_output<-FALSE

print(paste("Model run:",member))
print(paste("Start date:",startdate))
print(paste("End date:",enddate))

############## FIXED PARAMETERS ####################### #######################
# These are fixed for ADDTREES analyses - shouldn't need to change
collection<-'land-rcm'
domain<-'uk'
res="12km"

############## 2 PREPARE INPUTS ####################### #######################

### Area of interest and elevation data - AOI for downscaling defined by parcel data

# Load UK fine resolution dtm & boundary(coast)
dtmuk<-terra::rast(ukdtm_file) # 50m resolution dtm of all UK
coast_v<-terra::project(terra::vect(coast_file),dtmuk)

# Load parcels file - whole area used if raster otherwise checks geometries - converts to terra vect
if(tools::file_ext(parcels_file)=="tif"){
  parcels_v<-vect(ext(rast(parcels_file)))
  crs(parcels_v)<-crs(rast(parcels_file))
} else{
  parcels_sf<-st_read(parcels_file)
  num_invalid<-length(which(!st_is_valid(parcels_sf)))
  if(num_invalid>0) print(paste("Correcting ",num_invalid,"geometries in parcel file provided:",parcels_file))
  parcels_sf<-st_make_valid(parcels_sf)
  num_invalid<-length(which(!st_is_valid(parcels_sf)))
  if(num_invalid>0) warning(paste("There remain ",num_invalid,"invalid geometries after correction!!!"))
  parcels_v<-terra::vect(parcels_sf)
  if(parcel_output & !parcel_id %in% names(parcels_sf)) stop(paste("Cannot find parcel ID variable,",parcel_id,"among variables of",parcels_file))
}
parcels_v<-terra::project(parcels_v,dtmuk)

## Generate local area vector and required dtms
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi - assumes no coastal etc effects beyond 12km of parcel locations
dtmc<-mesoclimAddTrees::get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
dtmm <- terra::mask(terra::crop(terra::crop(dtmuk, dtmc), dtmuk), coast_v)

# Plot dtmf and overlay parcels
if(outputs){
  plot(dtmuk);plot(vect(ext(dtmc)),col="red",add=T)
  terra::plot(dtmc,main='DTMs')
  #terra::plot(dtmm,add=TRUE)
  terra::plot(dtmf,add=TRUE)
  terra::plot(parcels_v,add=TRUE)
}

# Size of area and number of parcels that it will downscale for
print(paste("Downscaling for an area of",round(expanse(aoi,"km"),0),"km^2, generating output for",nrow(parcels_v),"parcels"))


###### Load bias correction models if required
bcmodels_file<-file.path(dir_bcmodels,paste0('bias_correct_models_',member,'_2020_2022.Rds')) # will need to be after definition of model member
if(bias_correct) model_list<-readRDS(bcmodels_file)

###### Calculate topographical properties - worth it if looping over several downscaling calls (eg multiple years)
# Windshelter coef
t0<-now()
wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)

# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)
#terrain shading
results<-calculate_terrain_shading(dtmf,steps=24,toArrays=FALSE)
skyview<-results$skyview
horizon<-results$horizon
print(paste0("Topographical processing = ",now()-t0))

####### Prepare climate and seas surface data - prepare for whole time period if up to ~ 10 years
t0<-now()

# Process climate data from UKCP18 regional files on ceda archive - PROVIDE dtmc AND aoi as ONE
climdata<-ukcp18toclimarray(dir_ukcp,dtmc,startdate,enddate,collection,domain,member,temp_hgt=2, wind_hgt=2)

# Process sea surface temp data
sstdata<-create_ukcpsst_data(dir_sst,startdate,enddate,dtmc,member)

dataprep_time<-now()-t0
print(paste("Time for preparing data =", format(dataprep_time)))

if(outputs){
  plot(project(sstdata[[1]],crs(dtmc)))
  plot(dtmm,add=T)
  plot(aoi,add=TRUE)
}

# Check data - plot summary figs if requested
if(outputs){
  checkinputs(climdata, tstep = "day",plots=outputs)
}
print(paste0("Climate data processing = ",now()-t0))

if(inherits(climdata$relhum,"array")) climdata$relhum<-ifelse(climdata$relhum>100,100,climdata$relhum)
if(inherits(climdata$relhum,"SpatRaster")) climdata$relhum<-ifel(climdata$relhum>100,100,climdata$relhum)

############## BIAS CORRECTION AND SPATIAL DOWNSCALE BY YEAR AND MONTH ####################### #######################
years<-unique(c(year(startdate):year(enddate)))
# By YEAR
for (yr in years){
  t0<-now()
  yrstart<-as.POSIXlt(paste0(yr,'/01/01'),tz=tz(climdata$tme))
  yrend<-as.POSIXlt(paste0(yr,'/12/31'),tz=tz(climdata$tme))
  climdata_y<-subset_climdata(climdata,yrstart,yrend)

  ### BIAS CORRECTION of year's data if requested
  if(bias_correct){
    climdata_y<-biascorrect_climdata(climdata_y,model_list,prec_thold =0.01,fillna=TRUE)
    if(outputs) for(v in c('tmin','tmax','relhum','pres','swrad','lwrad','windspeed','winddir','prec')) plot(climdata_y[[v]][[1]],main=v)
    print(paste0("Bias correction processing = ",now()-t0))
  }

  # SPATIAL DOWNSCALE BY month
  out<-list()
  mnths<-unique(month(climdata_y$tme))
  for(m in mnths){
    print(m)
    sdatetime<-as.POSIXlt(paste0(yr,'/',m,'/01'),tz=tz(climdata$tme))
    dys<-lubridate::days_in_month(sdatetime)
    edatetime<-as.POSIXlt(paste0(yr,'/',sprintf("%02d", m),'/',sprintf("%02d", dys)),tz=tz(climdata$tme))
    climdata_m<-subset_climdata(climdata_y,sdatetime,edatetime)
    out[[length(out)+1]]<-spatialdownscale(climdata=climdata_m,
                                sst=sstdata,
                                dtmf=dtmf, dtmm=dtmm, basins = basins, wca=wca, skyview=skyview, horizon=horizon, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2, include_tmean=TRUE,
                                rhmin = 20, pksealevel = TRUE, patchsim = FALSE,
                                terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)
  }

  # Create whole year rasters for each var from list of months
  mesoclimate<-list()
  mesoclimate$dtm<-out[[1]]$dtm
  mesoclimate$windheight_m<-out[[1]]$windheight_m
  mesoclimate$tempheight_m<-out[[1]]$tempheight_m
  mesoclimate$tme<-do.call(c,(lapply(out,"[[","tme")))
  append_vars<-names(out[[1]])[c(5:length(out[[1]]))]
  for(v in append_vars) mesoclimate[[v]]<-rast(unlist(lapply(out,"[",v)))
  for(v in append_vars) names(mesoclimate[[v]])<-terra::time(mesoclimate[[v]])

  downscale_time<-now()-t0
  print(paste("Time for downscaling single year =", format(downscale_time)))

  # Plot mean,max and min days for each variable if spatraster
  if(outputs){
    climvars<-names(mesoclimate)[which(unlist(lapply(mesoclimate,inherits,what="SpatRaster")))]
    for(var in climvars){
      print(var)
      r<-mesoclimate[[var]]
      names(r)<-rep(var,nlyr(r))
      plot_q_layers(r,vtext=var)
    }
  }

  # Write yearly grids or parcel csvs
  if(mesoclim_output){
    mesoclim_file<-file.path(dir_out,paste0('mesoclimate_',res(dtmf)[1],'m_',arealabel,'_',yr,'.Rds'))
    print(paste("Writing mesoclimate file:",mesoclim_file))
    write_climdata(mesoclimate,mesoclim_file,overwrite=TRUE)
  }
  if(parcel_output){
    tp<-now()
    tmean<-(mesoclimate$tmax+mesoclimate$tmin)/2
    # Convert from relative to specific humidity
    spechum<-converthumidity(h=mesoclim:::.is(mesoclimate$relhum),
                             intype='relative',outtype='specific',
                             tc=mesoclim:::.is(tmean),
                             pk=mesoclim:::.is(mesoclimate$pres) )

    # Calculate and write parcel values -
    parcel_list<- create_parcel_list(mesoclimate,parcels_v,id='gid',output_spechum=TRUE)
    if(yr==years[1]) ov<-"replace" else ov<-"append"
    write_parcels(parcel_list, dir_out, overwrite=ov)
    print(paste("Time for parcel calculation and writing =", format(now()-tp)))
  }
  print(paste("Total yearly processing time =", format(now()-t0)))

}# end year
print(paste("Run time =", format(now()-tstart)))



############## Get a parcel variable and plot a map of it   ####################### #######################
#  var_sf<-get_parcel_var(mesoclimate,'tmean', parcels_v,id='gid', stat='mean' )
#  map_parcel_var(var_sf[which(st_is_valid(var_sf)),], plotvar='tmean', idvar='gid')



