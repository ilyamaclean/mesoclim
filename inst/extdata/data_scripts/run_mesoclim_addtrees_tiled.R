############## Code executable on LOCAL machine with downloaded inputs #######################
# Downloads input data from sharepoint to local drive dir - can all be deleted after modeling

############## LIBRARIES ####################### #######################
# install_gitgub("ilyamaclean/mesoclim")
# install_gitgub("jrmosedale/mesoclimAddTrees")
library(devtools)
library(terra)
library(sf)
library(mesoclimAddTrees)
library(lubridate)
library(mesoclim)

############## Downscale by TILES Function - see RCODE.R file ####################### #######################
spatialdownscale_tiles<-function(climdata, sst, dtmf, dtmm = NA, basins = NA, wca=NA, cad = TRUE,
                                 coastal = TRUE, thgto =2, whgto=2, rhmin = 20, pksealevel = TRUE,
                                 patchsim = TRUE, terrainshade = TRUE,precipmethod = "Elev",
                                 fast = TRUE, noraincut = 0, toArrays=FALSE, tilesize=c(100,100)){
  # Calculate tile set - imposes minimum ncol/rows per tile but size can vary
  xseq<-seq(1,ncol(dtmf),tilesize[1])
  yseq<-seq(1,nrow(dtmf),tilesize[2])
  minnum<-10 # Imposes min number of rows/cols & add tile or extends existing tile
  xdif<-ncol(dtmf)-xseq[length(xseq)]
  if(xdif>minnum) xseq<-c(xseq,ncol(dtmf)) else xseq[length(xseq)]<-xseq[length(xseq)]+xdif
  ydif<-nrow(dtmf)-xseq[length(yseq)]
  if(ydif>minnum) yseq<-c(yseq,nrow(dtmf)) else yseq[length(yseq)]<-yseq[length(yseq)]+ydif
  print(paste("Processing up to",(length(xseq)*length(yseq))/2,"tiles..."))

  climvars<-c('tmin','tmax','relhum','pres','swrad','lwrad','windspeed','winddir','prec','difrad')
  mesoclimate_tiles<-list()
  t<-1
  for(x in 1:(length(xseq)-1)){
    for(y in 1:(length(yseq)-1)){
      xmx<- xseq[x+1]
      xmn<- xseq[x]
      ymx<- yseq[y+1]
      ymn<- yseq[y]

      dtmf_tile<-dtmf[ymn:ymx,xmn:xmx,drop=FALSE]
      if(!all(is.na(values(dtmf_tile)))){ # Check it contains land cells (assumes already masked)
        print(paste("Downscaling tile",t))
        if(!is.logical(basins)) basins_tile<-basins[ymn:ymx,xmn:xmx,drop=FALSE]
        if(!is.logical(wca)) wca_tile<-wca[ ymn:ymx, xmn:xmx, ]
        # Downscale
        t0<-now()
        mesoclimate<-spatialdownscale(subset_climdata(climdata,sdatetime,edatetime), subset_climdata(sstdata,sdatetime,edatetime),
                                      dtmf_tile, dtmm, basins = basins_tile, wca=wca_tile, cad = TRUE,
                                      coastal = TRUE, thgto =2, whgto=2,
                                      rhmin = 20, pksealevel = TRUE, patchsim = TRUE,
                                      terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)

        # assign to list of tiles
        mesoclimate_tiles[[t]]<-mesoclimate
        downscale_time<-now()-t0
        print(paste("Time for downscaling single year for tile",t,"=", format(downscale_time)))
        t<-t+1
        rm(mesoclimate)
      }
    }
  }

  # Recreate single mesoclim list opf merged rasters
  mesoclimate<-list()
  mesoclimate$dtm<-do.call(merge, lapply(mesoclimate_tiles,`[[`, "dtm"))
  mesoclimate$tme<-mesoclimate_tiles[[1]]$tme
  mesoclimate$windheight_m<-mesoclimate_tiles[[1]]$windheight_m
  mesoclimate$tempheight_m<-mesoclimate_tiles[[1]]$tempheight_m
  for(v in climvars){
    mesoclimate[[v]]<-do.call(merge, lapply(mesoclimate_tiles,`[[`, v))
  }
  return(mesoclimate)
}

############## CLIMATE PARAMETERS ####################### #######################

# Start time for future climate timeseries.
startdate<-as.POSIXlt('2021/01/01')

# End time for future climate timeseries.
enddate<-as.POSIXlt('2029/12/31') # If using shared data folder use max value of as.POSIXlt('2039/12/31')

# Model run of UKCP18rcm to be downscaled.
member<-c('01')
# These are fixed for ADDTREES analyses - shouldn't need to change
collection<-'land-rcm'
domain<-'uk'
res="12km"
############## 1A INPUT FILES & DIRECTORIES ####################### #######################

# Any plot or print outputs? set to FALSE for Jasmin runs/slightly faster runs
#outputs<-FALSE
outputs<-TRUE

# Root directory - other dirs relative to this one or just use absolute paths below
dir_root<-"/Users/jonathanmosedale/Data"
dir_in<-file.path(dir_root,'mesoclim_inputs')
list.files(dir_in)


#### Input files and directories
#parcels_file<-file.path(dir_root,'lyme','land_parcels.shp') # elicitor app output file
parcels_file<-file.path(dir_in,'focal_parcels.shp') # elicitor app output file

file.exists(parcels_file)
coast_file<-file.path(dir_in,'CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
file.exists(coast_file)
ukdtm_file<-file.path(dir_in,"uk_dtm.tif") # 50m dtm  raster
file.exists(ukdtm_file)
ukcpdtm_file<-file.path(dir_root,"mesoclim_inputs","orog_land-rcm_uk_12km_osgb.nc")
file.exists(ukcpdtm_file)
dir_ukcp<-dir_in
dir.exists(dir_ukcp)
dir_sst<-dir_in
dir.exists(dir_sst)
# Directory for OUTPUTS - to which individual parcel .csv timeseries files are written.
dir_out<-file.path(dir_root,'mesoclim_outputs')  # output dir
#dir_out<-file.path(dir_root,'lyme')
dir.exists(dir_out)



############## 2 PREPARE INPUTS ####################### #######################

###### Area of interest and elevation data - AOI for downscaling defined by parcel data

# Load UK fine resolution dtm & boundary(coast)
dtmuk<-terra::rast(ukdtm_file) # 50m resolution dtm of all UK
coast_v<-terra::project(terra::vect(coast_file),dtmuk)

# Load parcels file and project to crs of output dtm (OS coords)
parcels_v<-terra::project(terra::vect(parcels_file),dtmuk)
#  parcels_v<-crop(parcels_v,c(260000,276719.799963001,35501.250040666,45000)) # To crop to about 1/3 of Peak district area

# Generate local area and dtm and wider extents
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi
dtmc<-get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
dtmm<-get_dtmm(dtmf,dtmc,dtmuk)

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


############## 3 SPATIAL DOWNSCALE - MAY REQUIRE TILING ####################### #######################
# If modelling many years will require looping for each year - appends new data to parcel file outputs with each loop
years<-unique(c(year(startdate):year(enddate)))
#yr<-years
for (yr in years){
  sdatetime<-as.POSIXlt(paste0(yr,'/01/01'))
  edatetime<-as.POSIXlt(paste0(yr,'/12/31'))

  # INTRODUCE TILING HERE !!!!
  t0<-now()
  mesoclimate<-spatialdownscale_tiles(subset_climdata(climdata,sdatetime,edatetime), subset_climdata(sstdata,sdatetime,edatetime),
                                dtmf, dtmm, basins = basins, wca=wca, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2,
                                rhmin = 20, pksealevel = TRUE, patchsim = TRUE,
                                terrainshade = TRUE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01,
                                tilesize=c(100,100))

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

