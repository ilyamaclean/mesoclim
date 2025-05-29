############## Code executable on LOCAL machine with downloaded inputs #######################
# Set up libraries files and directories for mac
setup_file<-"inst/extdata/data_scripts/setup_mesoclim_mac.R"
source(setup_file)

############## RUN PARAMETERS ####################### #######################
parcels_file<-file.path('/Users/jonathanmosedale/Downloads/land_parcels.shp 9/land_parcels.shp')
dir_out<-file.path(dir_root,'mesoclim_outputs')
member<-"01"
startyear<-"2022"
endyear<-"2022"
startdate<-as.POSIXlt(paste0(startyear,'/01/01'),tz="UTC")
enddate<-as.POSIXlt(paste0(endyear,'/12/31'),tz="UTC")

# Any plot or print outputs? set to FALSE for Jasmin runs
outputs<-TRUE

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

###### Area of interest and elevation data - AOI for downscaling defined by parcel data

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
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi
dtmc<-get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
dtmm<-get_dtmm(dtmf,dtmc,dtmuk)
# dtmw<-crop(dtmuk,dtmc)

# Plot dtmf and overlay parcels
if(outputs){
  terra::plot(dtmc,main='DTMs')
  terra::plot(dtmm,add=TRUE)
  terra::plot(dtmf,add=TRUE)
  terra::plot(parcels_v,add=TRUE)
}

# Size of area and number of parcels that it will downscale for
print(paste("Downscaling for an area of",round(expanse(aoi,"km"),0),"km^2, generating output for",nrow(parcels_v),"parcels"))


###### Calculate topographical properties - worth it if looping over several downscaling calls (eg multiple years)

# Windshelter coef
t0<-now()
wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)

# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)
print(now()-t0)


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
if(outputs) climdata<-checkinputs(climdata, tstep = "day")


############## 3 SPATIAL DOWNSCALE ####################### #######################
# If modelling many years will require looping for each year - appends new data to parcel file outputs with each loop
years<-unique(c(year(startdate):year(enddate)))
#yr<-years
for (yr in years){
  sdatetime<-as.POSIXlt(paste0(yr,'/01/01'))
  edatetime<-as.POSIXlt(paste0(yr,'/01/31'))

  # If cad = TRUE can take long time to calculate cold air drainage across large areas !!!
  t0<-now()
  mesoclimate<-spatialdownscale(climdata=subset_climdata(climdata,sdatetime,edatetime), sst=subset_climdata(sstdata,sdatetime,edatetime),
                                dtmf=dtmf, dtmm=dtmw, basins = basins, wca=wca, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2,
                                rhmin = 20, pksealevel = TRUE, patchsim = FALSE,
                                terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)

  downscale_time<-now()-t0
  print(paste("Time for downscaling single year =", format(downscale_time)))

  if(outputs){
    climvars<-c('tmin','tmax','relhum','pres','swrad','lwrad','windspeed','winddir','prec')
    for(var in climvars){
      print(var)
      r<-mesoclimate[[var]]
      names(r)<-rep(var,nlyr(r))
      plot_q_layers(r,vtext=var)
    }
  }
  #write_climdata(mesoclimate,file.path(dir_out,'mesoclimate_1yr_focal.Rds'))

  ###### Calculate and write parcel outputs with each yearly iteration

  # Convert from relative to specific humidity
  tmean<-(mesoclimate$tmax+mesoclimate$tmin)/2
  spechum<-converthumidity(h=mesoclim:::.is(mesoclimate$relhum),
                           intype='relative',outtype='specific',
                           tc=mesoclim:::.is(tmean),
                           pk=mesoclim:::.is(mesoclimate$pres) )

  # Calculate parcel values
  t0<-now()
  parcel_list<- create_parcel_list(mesoclimate,parcels_v,id='gid',output_spechum=TRUE)


  #write_parcels(parcel_list, dir_out, overwrite="replace")
  parcel_time<-now()-t0
  print(paste("Time for parcel calculation and writing =", format(parcel_time)))
  total_time<-now()-t0
  print(paste("Total parcel writing time =", format(parcel_time)))

} # end year loop


############## Get a parcel variable and plot a map of it   ####################### #######################
var_sf<-get_parcel_var(mesoclimate,'swrad', parcels_v,id='gid', stat='mean' )
map_parcel_var(var_sf[which(st_is_valid(var_sf)),], plotvar='swrad', idvar='gid')



############## CUT OUT: Download input files from sharepoint to local directory##############
library(Microsoft365R)
shrpt<-get_sharepoint_site('ADD-TREES Project')
data_shrpt<-shrpt$get_drive('Data Files')

# Filepath to vector file of parcels (produced by ellicitor app - defines AOI and Outputs)
data_shrpt$get_item("mesoclim_inputs/parcels")$download(dir_in)

# Filepath to UK coastline boundary polygon (not necessary if dtm already masked)
data_shrpt$get_item("mesoclim_inputs/boundaries")$download(dir_in)

# Filepath to fine resolution DTM of UK (OS Terrain50)
data_shrpt$get_item("mesoclim_inputs/dtm")$download(dir_in)

# Filepath to UKCP18 orography (coarse resolution DTM) file as downloaded from ceda
data_shrpt$get_item("mesoclim_inputs/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc")$download(file.path(dir_in,"orog_land-rcm_uk_12km_osgb.nc"))

# Directory for UKCP18 RCM input files
vars=c('clt','hurs','pr','prsn','psl','rls','rss','tasmax','tasmin','uas','vas')
vars=c('pr','prsn','psl','rls','rss','tasmax','tasmin','uas','vas')
decades<-.find_ukcp_decade(collection,startdate,enddate)

for(v in vars){
  f<-paste0(v,'_rcp85_',collection,'_',domain,'_',res,'_',member,'_day_',decades,'.nc')
  print(f)
  data_shrpt$get_item(paste0("mesoclim_inputs/ukcp18rcm/",f))$download(file.path(dir_in,f))
}

# Directory for UKCP18 Sea Surface temp input files
memberid<-ukcp18lookup$PP_ID[which(ukcp18lookup$Member_ID %in% member)]
start<-startdate %m-% months(1)
end<-enddate %m+% months(1)
yrs<-unique(c(year(start):year(end)))

sstfiles<-paste0("NWSClim_NWSPPE_",memberid,"_",yrs,"_gridT.nc")
for(f in sstfiles){
  print(f)
  data_shrpt$get_item(paste0("mesoclim_inputs/ukcp18sst/",f))$download(file.path(dir_in,f))
}
