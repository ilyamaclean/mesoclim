############## Code executable on LOCAL machine with downloaded inputs #######################
# Set up libraries files and directories for mac/jasmin
setup_file<-"inst/extdata/data_scripts/setup_mesoclim_mac.R"
source(setup_file)
tstart<-now()

############## RUN PARAMETERS ####################### #######################
#### Parcel file and label for outputs - can be a tif in which case parcel made of whole extent
# parcels_file<-file.path(dir_root,'CEH_Exmoor.shp') # exmoor
parcels_file<-file.path(dir_root,'mesoclim_inputs','killerton','land_parcels.shp') # killerton
# parcels_file<-file.path(dir_root,'mesoclim_inputs','focal_parcels.shp') # s devon
# parcels_file<-file.path(dir_root,'mesoclim_inputs','land_parcels.shp') # porthleven
# parcels_file<-file.path(dir_root,'mesoclim_inputs','wscotlandext.shp') # west scotland high elev dif & coast
# parcels_file<-file.path(system.file("extdata/dtms/dtmf_inland.tif",package="mesoclim")) # inland tif example

#### Label for outputs
arealabel<-"killerton"

#### UKCP options and time period
member<-"01"
startdate<-as.POSIXlt('2018/05/01',tz="UTC")
enddate<-as.POSIXlt('2018/05/31',tz="UTC")

#### Bias correction?
bias_correct<-TRUE

#### Output directory
dir_out<-"/Users/jonathanmosedale/Data/mesoclim_outputs"

#### Which outputs - graphical, parcel csv, mesoclim grids
outputs<-TRUE
parcel_output<-FALSE
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

  # SPATIAL DOWNSCALE BY year
  mesoclimate<-spatialdownscale(climdata=climdata_y,
                                sst=sstdata,
                                dtmf=dtmf, dtmm=dtmm, basins = basins, wca=wca, skyview=skyview, horizon=horizon, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2, include_tmean=TRUE,
                                rhmin = 20, pksealevel = TRUE, patchsim = FALSE,
                                terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)


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



