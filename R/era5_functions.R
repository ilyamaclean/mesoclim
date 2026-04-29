#' @title convert era4 ncdf4 file to format required for model
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by mcera5 pkg function request_era5 to the correct formal required for subsequent modelling.
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtmc a SpatRaster object of ERA5 elevations of same or larger extent as aoi and same resolution as ncfile (see details)
#' @param lsm  a SpatRaster object of ERA5 landsea mask of same or larger extent as aoi and same resolution as ncfile
#' @param aoi a  SpatRaster, sf or vect of the area to which outputs are cropped & reprojected - if NA dtmc used to define extent
#' @param startdate POSIXlt class defining start date of required timeseries
#' @param enddate POSIXlt class defining end date of required timeseries
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasters.
#' @param zo= height of wind speed above ground in metres to be output - defaults to ERA5 10m
#' @return a list of the following:
#' \describe{
#'    \item{dtmc}{Digital elevation of area (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2)}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed}{at 2m above ground (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @keywords preprocess era5
#' @details the model requires that input climate data are projected using a coordinate reference
#' system in which x and y are in metres. The output data are reprojected using the coordinate reference system and
#' extent of aoi (but SHOULD retain the approximate original grid resolution of the input climate data).
#' Suitable lsm can be derived from ERA5 landsea mask and dtmc can be derived from ERA5 geopotential ancillary data.
#' @examples
#'  \dontrun{
#' ncfile<-'path_to_downloaded_era5_file'
#' aoi<-terra::vect(terra::ext(-7.125,-2.875,49.375,51.625),crs='EPSG:4326')
#' dtm<-terra::rast(system.file('extdata/dtms/era5dtm.tif',package='mesoclim'))
#' era5input<-era5toclimarray(ncfile, dtm=NA, aoi=aoi)
#' plot_q_layers(.rast(era5input$temp,era5input$dtm))
#' checkinputs(era5input,'hour')
#' }
era5toclimarray <- function(ncfile, dtmc, lsm=NA, aoi=NA, startdate, enddate, dtr_cor_fac = 1.285, toArrays=TRUE, zo=10)  {
  # Check ERA5 variables and get time, check within start and end dates and ensure dtmc = res and crs of data
  nc <- nc_open(ncfile)
  era5vars <- names(nc$var)
  t2m<-rast(ncfile,"t2m")
  if("expver" %in% era5vars) tme<-as.POSIXlt(nc$var$expver$dim[[1]]$vals,tz='GMT',origin="1970-01-01")
  if(!"expver" %in% era5vars) tme<-as.POSIXlt(time(t2m), tz = "UTC")
  if(length(which(tme>=startdate & tme<=enddate))==0) stop("ERA5 data provided not within start and end dates!!")
  dtmc<-project(dtmc,crs(t2m))
  if(res(dtmc)[1]!=res(t2m)[1]){
    agf <- terra::res(t2m)[1]/terra::res(dtmc)[1]
    dtmc <- terra::aggregate(dtmc, fact = agf, fun = mean, na.rm = T)
  }
  nc_close(nc)

  # Crop dtmc to aoi if latter provided - check lsm presence if coastal correction
  if(is.na(dtr_cor_fac)) coastalcorrect<-FALSE else coastalcorrect<-TRUE
  if (class(aoi)[1] != "logical"){
    if (!class(aoi)[1] %in% c("SpatRaster", "SpatVector",
                              "sf"))
      stop("Parameter aoi NOT of suitable spatial class ")
    if (class(aoi)[1] == "sf") aoi <- vect(aoi)
    
    # Check dtmc and lsm extent equals or exceeds aoi
    aoiera5<-project(aoi,crs(dtmc))
    if(class(lsm)[1] == "logical" & (dtr_cor_fac!=0 | class(dtr_cor_fac)[1] == "logical")){
      warning("No lsm provided so ignoring coastal correction!!!")
      coastalcorrect<-FALSE
    }
    if(ext(dtmc)<ext(aoiera5)) stop("dtmc smaller than aoi!!!")
    if(class(lsm)[1] != "logical"  & ext(lsm)<ext(aoiera5)) stop("lsm smaller than aoi!!!")
    dtmc<-crop(dtmc,aoiera5,snap="out")
    if(class(lsm)[1] != "logical" ) lsm<-crop(lsm,aoiera5,snap="out")
  }
  units(dtmc)<-"m"
  names(dtmc)<-'Elevation'

  # Ensure startdata and enddate to POSIXlt values at start and end of day
  startdate<-as.POSIXlt(paste0(date(startdate)," 00:00:00"),tz=tz(startdate))
  enddate<-as.POSIXlt(paste0(date(enddate)," 23:00:00"),tz=tz(enddate))
  # Layers of data to extract for timeperiod
  lyrs<-which(tme>=startdate & tme<=enddate)

  # Load data as spatrasters
  t2m <- rast(ncfile, subds = "t2m")[[lyrs]] %>% crop(dtmc)
  d2m <- rast(ncfile, subds = "d2m")[[lyrs]] %>% crop(dtmc)
  if('sp' %in% era5vars) pres <- rast(ncfile, subds = "sp")[[lyrs]] %>% crop(dtmc) else if('msl' %in% era5vars){
    psl <- rast(ncfile, subds = "msl")[[lyrs]]
    pres<-psl * (((293-0.0065*dtmc)/293)^5.26) # convert to sea level pressure
  } else stop("Missing necessary pressure variables!!!")
  u10 <- rast(ncfile, subds = "u10")[[lyrs]] %>% crop(dtmc)
  v10 <- rast(ncfile, subds = "v10")[[lyrs]] %>% crop(dtmc)
  tp <- rast(ncfile, subds = "tp")[[lyrs]] %>% crop(dtmc)
  msdwlwrf <- rast(ncfile, subds = "msdwlwrf")[[lyrs]] %>% crop(dtmc)
  if (all(c("fdir", "ssrd") %in% era5vars)) {
    fdir<-rast(ncfile, subds = "fdir")[[lyrs]] %>% crop(dtmc)
    fdir <- fdir/3600
    ssrd<-rast(ncfile, subds = "ssrd")[[lyrs]] %>% crop(dtmc)
    ssrd <- ssrd/3600
  } else if (all(c("msdwswrf", "msdrswrf") %in% era5vars)) {
    fdir <- rast(ncfile, subds = "msdwswrf")[[lyrs]] %>% crop(dtmc)
    ssrd <- rast(ncfile, subds = "msdrswrf")[[lyrs]] %>% crop(dtmc)
  } else stop("Missing necessary SW radiation variables!!!")

  # Convert temp to Celsius
  t2m <- t2m - 273.15

  # Coastal correction of temperature
  if(coastalcorrect){
    lsm_e<-crop(lsm,dtmc)
    a<-as.array(rep(lsm_e,dim(t2m)[3]))
    mu<-(1-a)*dtr_cor_fac+1
    tmn <- .ehr(as.array(.hourtoday(t2m, mean)))
    tdif <- (as.array(t2m) - tmn) * mu
    tc <- tmn + tdif
    tc <- .rast( tc, t2m)
  } else tc<-t2m

  # Reproject variables to original aoi if different from dtmc/era5 data
  if(class(aoi)[1] != "logical"){
    tc <- terra::project(tc, crs(aoi))
    d2m<- terra::project(d2m,crs(aoi))
    pres <- terra::project(pres, crs(aoi))
    u10 <- terra::project(u10, crs(aoi))
    v10 <- terra::project(v10, crs(aoi))
    tp <- terra::project(tp, crs(aoi))
    msdwlwrf <- terra::project(msdwlwrf, crs(aoi))
    fdir <- terra::project(fdir, crs(aoi))
    ssrd <- terra::project(ssrd, crs(aoi))
    dtmc<-project(dtmc,crs(aoi))
  }
  # Calculate output variables
  t0<-now()
  ea <- .satvap(d2m-273.15)
  print(now()-t0)
  #temp <- as.array(tc)
  relhum <- (ea/.satvap(tc)) * 100
  relhum<-ifel(relhum>100,100,relhum)
  pres <- as.array(pres)/1000 ## Surface Pressure!!!!
  swrad <- as.array(ssrd)
  difrad <- swrad - as.array(fdir)
  lwrad <- as.array(msdwlwrf)
  # Calculate windspeed adjusting to output height
  windspeed <- sqrt(as.array(u10)^2 + as.array(v10)^2) * log(67.8 * zo - 5.42)/log(67.8 * 10 - 5.42)
  winddir <- as.array((terra::atan2(u10, v10) * 180/pi + 180)%%360)
  prec <- as.array(tp) * 1000

  # Format outputs
  out <- list(dtm = dtmc, tme = tme, windheight_m = zo, tempheight_m = 2)
  climout <- list(temp = tc, relhum = relhum, pres = pres,
                  swrad = swrad, difrad = difrad, lwrad = lwrad, windspeed = windspeed,
                  winddir = winddir, prec = prec)
  climunits<-c('degC','%','kPa','watt/m^2','watt/m^2','watt/m^2','m/s','deg','mm')
  # Round values
  climround<-c(2,1,1,1,1,1,3,1,4)
  for(n in 1:length(climout)){
    v <- names(climout)[n]
    climout[[v]]<-round(climout[[v]],climround[n])
  }
  # For spatrasters add units and time
  if (toArrays == FALSE){
    climout <- lapply(climout, .rast, tem = dtmc)
    for(n in 1:length(climout)){
      v <- names(climout)[n]
      u<-climunits[n]
      terra::time(climout[[v]])<-tme
      names(climout[[v]])<-tme
      units(climout[[v]])<-u
    }
  }
  out <- c(out, climout)
  return(out)
}
