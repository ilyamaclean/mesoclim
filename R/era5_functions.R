#' @title convert era4 ncdf4 file to format required for model
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by [mcera5::request_era5()] to the correct formal required for subsequent modelling.
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtm a SpatRaster object of elevations covering the extent of the study area (see details)
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasters.
#' @return a list of the following:
#' \describe{
#'    \item{dtm}{Digital elevation of downscaled area in metres (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2)}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @details the model requires that input climate data are projected using a coordinate reference
#' system in which x and y are in metres. Since values returned by [mcera5::request_era5()]
#' are in lat long, the output data are reprojected using the coordinate reference system and
#' extent of dtm (but retain the approximate original grid resolution of the input climate data).
#' Returned climate data match the resolution, coordinate reference system and extent of `dtmc`.
#' @examples
#' ncfile<-system.file('extdata/era5/era5_surface_ukeire_2018_05.nc',package='mesoclim')
#' aoi<-terra::vect(terra::ext(-7.125,-2.875,49.375,51.625),crs='EPSG:4326')
#' dtm<-terra::rast(system.file('extdata/dtms/era5dtm.tif',package='mesoclim'))
#' era5input<-era5toclimarray(ncfile, dtm=NA, aoi=aoi)
#' plot_q_layers(.rast(era5input$temp,era5input$dtm))
#' checkinputs(era5input,'hour')
era5toclimarray <- function(ncfile, dtm=NA, aoi=NA, dtr_cor_fac = 1.285, toArrays=TRUE)  {
  # Get extent of aoi and project to same as ukcp data (lat lon) - check if any provided dtm covers
  if(class(aoi)[1]!='logical'){
    if(!class(aoi)[1] %in% c("SpatRaster","SpatVector","sf")) stop("Parameter aoi NOT of suitable spatial class ")
    if(class(aoi)[1]=="sf") aoi<-vect(aoi)
    target_crs<-crs(rast(ncfile[1]))
    aoi_e<-ext(terra::project(aoi,target_crs) )
    if(class(dtm)[1]!='logical'){
      dtm_e<-terra::ext(project(dtm,target_crs))
      if(aoi_e$xmin<dtm_e$xmin|aoi_e$ymin<dtm_e$ymin|aoi_e$xmax>dtm_e$xmax|aoi_e$ymax>dtm_e$ymax) stop('DTM provided small than aoi provided!!')
    }
  }

  # Get variables
  t2m<-rast(ncfile,subds = "t2m") # Air temperature (K)
  d2m<-rast(ncfile,subds = "d2m") # Dewpoint temperature (K)
  sp<-rast(ncfile,subds = "sp") # Surface pressure (Pa)
  u10<-rast(ncfile,subds = "u10") # U-wind at 10m (m/s)
  v10<-rast(ncfile,subds = "v10") # V-wind at 10m (m/s)
  tp<-rast(ncfile,subds = "tp") # Total precipitation (m)
  #msnlwrf<-rast(ncfile,subds = "msnlwrf")   # Mean surface net long-wave radiation flux (W/m^2)
  msdwlwrf<-rast(ncfile,subds = "msdwlwrf") # Mean surface downward long-wave radiation flux (W/m^2)
  fdir<-rast(ncfile,subds = "fdir") #  Total sky direct solar radiation at surface (W/m^2)
  ssrd<-rast(ncfile,subds = "ssrd") # Surface short-wave (solar) radiation downwards (W/m^2)
  lsm<-rast(ncfile,subds = "lsm") # Land sea mask
  geop<-rast(ncfile,subds = "z") # Geopotential

  # Create coarse-resolution dtm if dtm provided otherwise derive from geopotential
  if(class(dtm)[1]!='logical'){
    te<-terra::project(t2m[[1]],terra::crs(dtm))
    agf<-terra::res(te)[1]/terra::res(dtm)[1]
    dtmc<-terra::aggregate(dtm,fact=agf,fun=mean,na.rm=T)
  } else{
    dtmc<-geop[[1]]/9.80665
    dtmc<-ifel(dtmc<0,0,dtmc)
  }
  # IF aoi provided project dtm to same projection and extent
  if(class(aoi)[1]!='logical') dtmc<-terra::crop(terra::project(dtmc,crs(aoi)),aoi)

  # Apply coastal correction to temperature data
  tmn<-.ehr(.hourtoday(as.array(t2m)-273.15,min))
  mu<-(1-as.array(lsm))*dtr_cor_fac+1
  tc<-.rast(((as.array(t2m)-273.15)-tmn)*mu+tmn,t2m)

  # Calculate vapour pressure
  ea<-.rast(.satvap(as.array(d2m)-273.15),t2m)

  # Project and crop all variables to match dtmc
  tc<-terra::project(tc,dtmc)
  ea<-terra::project(ea,dtmc)
  sp<-terra::project(sp,dtmc)
  u10<-terra::project(u10,dtmc)
  v10<-terra::project(v10,dtmc)
  tp<-terra::project(tp,dtmc)
  #msnlwrf<-terra::project(msnlwrf,dtmc)
  msdwlwrf<-terra::project(msdwlwrf,dtmc)
  fdir<-terra::project(fdir,dtmc)
  ssrd<-terra::project(ssrd,dtmc)

  # Derive varies
  temp<-as.array(tc) # Temperature (deg c)
  relhum<-(as.array(ea)/.satvap(temp))*100 # Relative humidity (%)
  pres<-as.array(sp)/1000  # SEa-level surface pressure (kPa)
  swrad<-as.array(ssrd)/3600 # Downward shortwave radiation (W/m^2)
  difrad<-swrad-as.array(fdir)/3600 # Downward diffuse radiation (W/m^2)
  lwrad<-as.array(msdwlwrf) # Downward longwave radiation
  windspeed<-sqrt(as.array(u10)^2+as.array(v10)^2)*log(67.8*2-5.42)/log(67.8*10-5.42) # Wind speed (m/s)
  winddir<-as.array((terra::atan2(u10,v10)*180/pi+180)%%360) # Wind direction (deg from N)
  prec<-as.array(tp)*1000
  # Generate POSIXlt object of times
  tme<-as.POSIXlt(time(t2m),tz="UTC")
  # Create output list
  out<-list(dtm=dtmc, tme=tme,windheight_m=10,tempheight_m=2)
  climout<-list(temp=temp,relhum=relhum,pres=pres,swrad=swrad,difrad=difrad,
                lwrad=lwrad,windspeed=windspeed,winddir=winddir,prec=prec)
  # Convert climate data to SpatRasters if requested
  if(toArrays==FALSE) climout<-lapply(climout,.rast,tem=dtmc)
  # Output for returning
  out<-c(out,climout)
  return(out)
}
