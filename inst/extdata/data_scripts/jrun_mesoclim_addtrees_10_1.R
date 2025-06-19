#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############## Check User PARAMETERS passed via Rscript line command  ####################### #######################
# 1 parcel_file location
# 2 output directory
# 3 model memeber number eg "01"
# 4 start year
# 5 end year
if (length(args)==0) {
  stop("At least one argument must be supplied (parcel file)", call.=FALSE)
} else if (length(args)>5) {
  stop("Up to four arguments are required (parcel file, output dir, starttime, endtime)", call.=FALSE) }

if(!file.exists(args[1])) stop("Input file provided does NOT exist!!!")

maxyear<-"2080"
minyear<-"1981"
maxdate<-as.POSIXlt(paste0(maxyear,'/11/30'),tz="UTC")

out_default<-"/gws/nopw/j04/uknetzero/mesoclim/mesoclim_outputs"

if (length(args)<5) {
  warning("Using default end times!!")
  args[5] = 2080
}
if (length(args)<4) {
  warning("Using default start and end times!!")
  args[4] = 2024
}
if (length(args)<3){
  warning("Using default UKCP modelrun (01)!!")
  args[3] = "01"

}
if (length(args)<2) {
  warning("Using default output directory!!")
  args[2] = out_default
}

if(!dir.exists(args[2])){
  warning(paste("Output directory provided does NOT exist, using default dir",out_default))
  args[2]<-out_default
}
if(args[5]>maxyear){
  warning(paste("End year too high - changing to",maxyear))
  args[5]<-maxyr
}
if(args[4]<minyear){
  warning(paste("Start year too low - changing to",minyear))
  args[4]<-minyr
}

############## Assign PARAMETERS  ####################### #######################
parcels_file<-args[1]
dir_out<-args[2]
member<-args[3]
startyear<-args[4]
endyear<-args[5]
ftr_sdate<-as.POSIXlt(paste0(startyear,'/01/01'),tz="UTC")
ftr_edate<-as.POSIXlt(paste0(endyear,'/12/31'),tz="UTC")

print(paste("Parcels input file:",parcels_file))
print(paste("Output directory:",dir_out))
print(paste("Model run:",member))
print(paste("Start date:",ftr_sdate))
print(paste("End date:",ftr_edate))

############## Define Other Options  ####################### #######################
# Print/Plot outputs
outputs<-FALSE
# Parcel climate data output (one file per parcel)
parcel_output<-TRUE
# Yearly .tif files of climate rasters
mesoclim_output<-FALSE

### Bias correction - TRUE for most Jasmin runs
bias_correct<-TRUE

############## LIBRARIES ####################### #######################
dir_lib<-"/gws/nopw/j04/uknetzero/mesoclim/mesoclim_lib"
library(terra)
library(sf)
library(lubridate)
library(mesoclim, lib.loc=dir_lib)
library(mesoclimAddTrees, lib.loc=dir_lib)

terraOptions(tempdir = "/gws/nopw/j04/uknetzero/mesoclim/terra_storage")

############## 1A INPUT FILES & DIRECTORIES ####################### #######################
tstart<-now()

# basepath to badc/... ceda data can be set for testing but use "" for runs on Jasmin
ceda_basepath <-""

# Root directory relative to these data inputs
dir_root<-"/gws/nopw/j04/uknetzero/mesoclim"

# Filepath to vector file of parcels output by ellicitor app.
if(!file.exists(parcels_file)) stop("Cannot find parcels input file!!")

# Filepath to coarse res DTM matching UKCP data
ukcpdtm_file<-file.path(ceda_basepath,"badc","ukcp18","data","land-rcm","ancil","orog","orog_land-rcm_uk_12km_osgb.nc")
if(!file.exists(ukcpdtm_file)) stop("Cannot find UKCP18 DTM file for the UK!!")

# Filepath to fine resolution DTM of UK (OS Terrain50 -  which has advantages in having a .vrt raster that can be queried to enhance cropping and extraction of relevant area to match aoi etc. Subdirectories within hold actual tiled data.)
ukdtm_file<-file.path(dir_root,'mesoclim_inputs','dtm',"uk_dtm.tif") # 50m dtm virtual raster
if(!file.exists(ukdtm_file)) stop("Cannot find DTM file for the UK!!")

# Filepath to coastline boundary polygon (not necessary if dtm already masked)
coast_file<-file.path(dir_root,'mesoclim_inputs','boundaries','CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
if(!file.exists(coast_file)) stop("Cannot find coastal mask file for the UK!!")

# Directory for outputs - to which individual parcel .csv timeseries files are written.
if(!dir.exists(dir_out)) stop("Cannot find output directory!!")

# Directory holding bias correction models
dir_bcmodels<-file.path(dir_root,"mesoclim_inputs","bc_models")
if(!dir.exists(dir_bcmodels)) stop("Cannot find bias correction model directory!!")

# Remove any existing parcel files in dir_out
#pfiles<-list.files(dir_out,full.names=TRUE, pattern="parcel")
#file.remove(pfiles)

############## 2 PREPARE SPATIAL INPUTS ####################### #######################

# Load UK fine resolution dtm & boundary(coast) mask at 50m resolution
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

# If parcel_output write a list of parcel ID and OSGB XY coordinates of centroid
if(parcel_output){
  parcels_sf<-st_as_sf(parcels_v)
  parcel_centroids <- round(st_coordinates(st_centroid(st_make_valid(parcels_sf))),1)
  parcels_txt<-data.frame("id"=parcels_sf$gid,"x"=parcel_centroids[,"X"],"y"=parcel_centroids[,"Y"])
  write.table(parcels_txt, file.path(dir_out,"parcel_ids.csv"), sep = ",", row.names = FALSE,  col.names = TRUE, quote=FALSE)
}

## Generate local area vector and required dtms
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi - assumes no coastal etc effects beyond 12km of parcel locations
dtmc<-mesoclimAddTrees::get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate dtm for coatal/wind effects - same extent as dtmc
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

###### Calculate constant topographical properties   ####################### #######################
# Windshelter coef
wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)
# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)
# Terrain shading
results<-calculate_terrain_shading(dtmf,steps=24,toArrays=FALSE)
skyview<-results$skyview
horizon<-results$horizon
print(paste0("Topographical processing = ",now()-t0))

###### Load bias correction models if required  ####################### #######################
if(bias_correct){
  bcmodels_file<-file.path(dir_bcmodels,paste0('bias_correct_models_',member,'_2020_2022.Rds')) # will need to be after definition of model member
  message(paste("Loading bias correction models file",bcmodels_file))
  model_list<-readRDS(bcmodels_file)
}

dataprep_time<-now()-tstart
print(paste("Time for preparing data =", format(dataprep_time)))

############## 3a Prepare climate  data in 10 yr loop ####################### #######################

yr10seq<-seq(floor(year(ftr_sdate)),floor(year(ftr_edate)),10)

for(start_year in yr10seq){
  t0<-now()
  end_year<-min(c(start_year+9,year(ftr_edate)))
  startdate<-as.POSIXlt(paste0(start_year,'/01/01'))
  enddate<-as.POSIXlt(paste0(end_year,'/12/31'))

  # Process climate data from UKCP18 regional files on ceda archive
  climdata<-addtrees_climdata(dtmc=dtmc,startdate=startdate,enddate=enddate,collection='land-rcm',domain='uk',member=member,basepath=ceda_basepath)

  # Process sea surface temo data from ceda archive ??CHANGE OUTPUT TO PROJECTION OF DTMC/AOI? COMBINE WITH ABOVE?
  sstdata<-addtrees_sstdata(startdate=startdate,enddate=enddate,dtmc=climdata$dtm,member='01',basepath=ceda_basepath)

  if(outputs){
    plot(project(sstdata[[1]],crs(dtmc)))
    plot(dtmm,add=T)
    plot(aoi,add=TRUE)
  }

  # Check data - plot summary figs
  if(outputs) climdata<-checkinputs(climdata, tstep = "day")

  ############## 3b SPATIAL DOWNSCALE by year ####################### #######################

  ## Yearly downscaling
  years<-seq(start_year,end_year,1)

  for (yr in years){
    t0<-now()
    sdatetime<-as.POSIXlt(paste0(yr,'/01/01'))
    edatetime<-as.POSIXlt(paste0(yr,'/12/31'))
    climdata_y<-subset_climdata(climdata,sdatetime,edatetime)

    mesoclimate<-spatialdownscale(climdata=climdata_y,
                                  sst=sstdata,
                                  dtmf=dtmf, dtmm=dtmm, basins = basins, wca=wca, skyview=skyview, horizon=horizon, cad = TRUE,
                                  coastal = TRUE, thgto =2, whgto=2, include_tmean=TRUE,
                                  rhmin = 20, pksealevel = TRUE, patchsim = FALSE,
                                  terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)

    if(outputs){
      climvars<-c('tmin','tmax','tmean','relhum','pres','swrad','lwrad','windspeed','winddir','prec')
      for(var in climvars){
        print(var)
        r<-mesoclimate[[var]]
        names(r)<-rep(var,nlyr(r))
        plot_q_layers(r,vtext=var)
      }
    }

    ############## Outputs by Year ####################### #######################

    # Write yearly grids or parcel csvs
    if(mesoclim_output){
      mesoclim_file<-file.path(dir_out,paste0('mesoclimate_',res(dtmf)[1],'m_',arealabel,'_',yr,'.Rds'))
      print(paste("Writing mesoclimate file:",mesoclim_file))
      write_climdata(mesoclimate,mesoclim_file,overwrite=TRUE)
    }
    ##  Calculate and write parcel outputs
    if(parcel_output){
      # Calculate and write parcel values -
      parcel_list<- create_parcel_list(mesoclimate,parcels_v,id='gid',output_tmean=TRUE,output_spechum=TRUE)
      if(yr==years[1]) ov<-"replace" else ov<-"append"
      write_parcels(parcel_list, dir_out, overwrite=ov)
      print(paste("Time for parcel calculation and writing =", format(now()-tp)))
    }

    year_time<-now()-t0
    print(paste("Time for downscaling year",yr,"=", format(year_time)))
  } # yr

}# ten years

total_time<-now()-tstart
print(paste("Total time for completion =", format(total_time)))

#  var_sf<-get_parcel_var(mesoclimate,'tmax', parcels_v,id='gid', stat='mean' )
#  map_parcel_var(var_sf, plotvar='tmax', idvar='gid')
