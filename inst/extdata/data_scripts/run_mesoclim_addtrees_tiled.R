############## Code executable on LOCAL machine with downloaded inputs #######################
# Set up libraries files and directories for mac/jasmin
setup_file<-"inst/extdata/data_scripts/setup_mesoclim_mac.R"
source(setup_file)

tstart<-now()
dir_bcmodels<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/bias_correct_models"

############## RUN PARAMETERS ####################### #######################
# Parcel file and label for outputs
parcels_file<-file.path(dir_root,'mesoclim_inputs','focal_parcels.shp') # s devon
arealabel<-"sdevon"

# UKCP options and time period
member<-"01"
startyear<-"2022"
endyear<-"2022"
startdate<-as.POSIXlt(paste0(startyear,'/01/01'),tz="UTC")
enddate<-as.POSIXlt(paste0(endyear,'/12/31'),tz="UTC")

# Bias correction?
bias_correct<-TRUE
bcmodels_file<-file.path(dir_bcmodels,'bias_correct_models_01_2020_2022.Rds') # will need to be after definition of model member

# Which outputs - graphical, parcel csv, mesoclim grids
dir_out<-file.path(dir_root,'mesoclim_outputs')
outputs<-TRUE
parcel_output<-TRUE
mesoclim_output<-TRUE

print(paste("Parcels input file:",parcels_file))
file.exists(parcels_file)
print(paste("Output directory:",dir_out))
dir.exists(dir_out)
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

# Load parcels file - check geometries and project to crs of output dtm (OS coords)
parcels_sf<-st_read(parcels_file)
num_invalid<-length(which(!st_is_valid(parcels_sf)))
if(num_invalid>0) print(paste("Correcting ",num_invalid,"geometries in parcel file provided:",parcels_file))
parcels_sf<-st_make_valid(parcels_sf)
num_invalid<-length(which(!st_is_valid(parcels_sf)))
if(num_invalid>0) warning(paste("There remain ",num_invalid,"invalid geometries after correction!!!"))

parcels_v<-terra::project(terra::vect(parcels_sf),dtmuk)

# Generate local area and dtm and wider extents
#. aoi<-vect(ext(c(255000,265000,43000,65000)))
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi
dtmc<-get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
#dtmm_res <- round(exp((log(terra::res(dtmc)[1]) + log(terra::res(dtmf)[1]))/2))
dtmm <- terra::mask(terra::crop(terra::crop(dtmuk, dtmc), dtmuk), coast_v)
#dtmm<-get_dtmm(dtmf,dtmc,dtmuk)

# Plot dtmf and overlay parcels
if(outputs){
  terra::plot(dtmc,main='DTMs')
  #terra::plot(dtmm,add=TRUE)
  terra::plot(dtmf,add=TRUE)
  terra::plot(parcels_v,add=TRUE)
}

# Size of area and number of parcels that it will downscale for
print(paste("Downscaling for an area of",round(expanse(aoi,"km"),0),"km^2, generating output for",nrow(parcels_v),"parcels"))


###### Load bias correction models if required
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

# Process sea surface temo data from ceda archive ??CHANGE OUTPUT TO PROJECTION OF DTMC/AOI? COMBINE WITH ABOVE?
sstdata<-create_ukcpsst_data(dir_sst,startdate,enddate,dtmc,member)

dataprep_time<-now()-t0
print(paste("Time for preparing data =", format(dataprep_time)))

if(outputs){
  plot(project(sstdata[[1]],crs(dtmc)))
  plot(dtmm,add=T)
  plot(aoi,add=TRUE)
}

# Check data - plot summary figs
if(outputs){
  #climdata<-checkinputs(climdata, tstep = "day")
}
print(paste0("Climate data processing = ",now()-t0))



############## 3 SPATIAL DOWNSCALE WITH TILING ####################### #######################
# Option to write yearly .tifs
# Option to write parcel .csv outputs  - appends new data to parcel file outputs with each loop
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
  # BY MONTHLY AND BY TILE
  t0<-now()
  mesoclimate<-spatialdownscale_tiles(climdata_y, sstdata, dtmf, dtmm,
                                basins = basins, wca=wca, skyview=skyview, horizon=horizon,
                                cad = TRUE, coastal = TRUE, thgto =2, whgto=2, include_tmean=TRUE,
                                rhmin = 20, pksealevel = TRUE, patchsim = FALSE,terrainshade = TRUE,
                                precipmethod = "Elev", fast = TRUE, noraincut = 0.01,
                                toArrays=FALSE, overlap=1000, sz=10000)
  downscale_time<-now()-t0
  print(paste("Time for downscaling single year =", format(downscale_time)))

  if(outputs){
    climvars<-c('tmean','tmin','tmax','relhum','pres','swrad','lwrad','windspeed','winddir','prec')
    for(var in climvars){
      print(var)
      r<-mesomonth[[var]]
      names(r)<-rep(var,nlyr(r))
      plot_q_layers(r,vtext=var)
    }
  }
  #write_climdata(mesoclimate,file.path(dir_out,'mesoclimate_1yr_focal.Rds'))

  ###### Calculate and write parcel outputs with each yearly iteration

  # Calculate parcel values
  t0<-now()
  parcel_list<- create_parcel_list(mesoclimate,parcels_v,id='gid',output_spechum=TRUE)
  # Change to append except for first year!!!!!
  if(yr==year(startdate)) write_parcels(parcel_list, dir_out, overwrite='replace') else  write_parcels(parcel_list, dir_out, overwrite='append')

  parcel_time<-now()-t0
  print(paste("Time for parcel calculation and writing =", format(parcel_time)))
  total_time<-now()-t0
  print(paste("Total parcel writing time =", format(parcel_time)))

} # end year loop

##############  Write a list of parcel ID and OSGB XY coordinates of centroid  #######################
parcels_sf<-st_as_sf(parcels_v)
parcel_centroids <- round(st_coordinates(st_centroid(st_make_valid(parcels_sf))),1)
parcels_txt<-data.frame("id"=parcels_sf$gid,"x"=parcel_centroids[,"X"],"y"=parcel_centroids[,"Y"])
write.table(parcels_txt, file.path(dir_out,"parcel_ids.csv"), sep = ",", row.names = FALSE,  col.names = TRUE, quote=FALSE)

############## Get a parcel variable and plot a map of it   #######################
var_sf<-get_parcel_var(mesoclimate,'swrad', parcels_v,id='gid', stat='mean' )
map_parcel_var(var_sf, plotvar='swrad', idvar='gid')

