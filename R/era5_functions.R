#' @title convert era4 ncdf4 file to format required for model
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by mcera5 pkg function request_era5 to the correct formal required for subsequent modelling.
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtmc a SpatRaster object of ERA5 elevations of same or larger extent as aoi and same resolution as ncfile (see details)
#' @param lsm  a SpatRaster object of ERA5 landsea mask of same or larger extent as aoi and same resolution as ncfile
#' @param aoi a  SpatRaster, sf or vect of the area to which outputs are cropped & reprojected - if NA extent/projection same as ncfile
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasters.
#' @param zo= height of wind speed above ground in metres to be output
#' @return a list of the following:
#' \describe{
#'    \item{dtmc}{Digital elevation of downscaled area in metres (as Spatraster)}
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
era5toclimarray <- function(ncfile, dtmc, lsm, aoi=NA, dtr_cor_fac = 1.285, toArrays=TRUE, zo=2)  {
  # Get extent of aoi and check equal or smaller than dtmc and lsm inputs
  if (class(aoi)[1] != "logical"){
    if (!class(aoi)[1] %in% c("SpatRaster", "SpatVector",
                              "sf"))
      stop("Parameter aoi NOT of suitable spatial class ")
    if (class(aoi)[1] == "sf") aoi <- vect(aoi)
    # Check dtmc and lsm extent equals or exceeds aoi
    if(ext(dtmc)<ext(project(aoi,crs(dtmc)))) stop("dtmc smaller than aoi!!!")
    if(ext(lsm)<ext(project(aoi,crs(dtmc)))) stop("dtmc smaller than aoi!!!")
  }
  # If no aoi assume same as dtmc
  if(class(aoi)[1]=="logical") aoi <- dtmc

  # Check vars and get time
  nc <- nc_open(ncfile)
  era5vars <- names(nc$var)
  if("expver" %in% era5vars) tme<-as.POSIXlt(nc$var$expver$dim[[1]]$vals,tz='GMT',origin="1970-01-01")
  if(!"expver" %in% era5vars) tme<-as.POSIXlt(time(rast(ncfile,"t2m")), tz = "UTC")
  nc_close(nc)

  t2m <- rast(ncfile, subds = "t2m")
  d2m <- rast(ncfile, subds = "d2m")
  sp <- rast(ncfile, subds = "sp")
  u10 <- rast(ncfile, subds = "u10")
  v10 <- rast(ncfile, subds = "v10")
  tp <- rast(ncfile, subds = "tp")
  msdwlwrf <- rast(ncfile, subds = "msdwlwrf")
  if (all(c("fdir", "ssrd") %in% era5vars)) {
    fdir <- rast(ncfile, subds = "fdir")/3600
    ssrd <- rast(ncfile, subds = "ssrd")/3600
  } else if (all(c("msdwswrf", "msdrswrf") %in% era5vars)) {
    fdir <- rast(ncfile, subds = "msdwswrf")
    ssrd <- rast(ncfile, subds = "msdrswrf")
  } else stop("Missing necessary SW radiation variables!!!")

  # Aggregate dtmc if not of same resolution as ncfile data
  if(res(dtmc)[1]!=res(t2m)[1]){
    agf <- terra::res(t2m)[1]/terra::res(dtmc)[1]
    dtmc <- terra::aggregate(dtmc, fact = agf, fun = mean, na.rm = T)
  }

  # Reproject and crop all variables to aoi
  t2m <- terra::project(t2m, aoi)
  d2m<- terra::project(d2m, aoi)
  sp <- terra::project(sp, aoi)
  u10 <- terra::project(u10, aoi)
  v10 <- terra::project(v10, aoi)
  tp <- terra::project(tp, aoi)
  msdwlwrf <- terra::project(msdwlwrf, aoi)
  fdir <- terra::project(fdir, aoi)
  ssrd <- terra::project(ssrd, aoi)
  dtmc<-project(dtmc,aoi)
  lsm<-project(lsm,aoi)

  # Convert temp to Celsius
  t2m <- t2m - 273.15

  # Coastal correction of temperature
  tmn <- .ehr(.hourtoday(as.array(t2m), min))
  a<-as.array(rep(lsm,dim(tmn)[3]))
  mu<-(1-a)*dtr_cor_fac+1
  tc <- .rast(((as.array(t2m)) - tmn) * mu + tmn, t2m)

  # Calculate output variables
  ea <- .rast(.satvap(as.array(d2m)-273.15), tc)
  temp <- as.array(tc)
  relhum <- (as.array(ea)/.satvap(temp)) * 100
  pres <- as.array(sp)/1000
  swrad <- as.array(ssrd)
  difrad <- swrad - as.array(fdir)
  lwrad <- as.array(msdwlwrf)
  windspeed <- sqrt(as.array(u10)^2 + as.array(v10)^2) * log(67.8 *
                                                               zo - 5.42)/log(67.8 * 10 - 5.42)
  winddir <- as.array((terra::atan2(u10, v10) * 180/pi + 180)%%360)
  prec <- as.array(tp) * 1000
  out <- list(dtm = dtmc, tme = tme, windheight_m = 10, tempheight_m = 2)
  climout <- list(temp = temp, relhum = relhum, pres = pres,
                  swrad = swrad, difrad = difrad, lwrad = lwrad, windspeed = windspeed,
                  winddir = winddir, prec = prec)

  if (toArrays == FALSE) climout <- lapply(climout, .rast, tem = dtmc)

  out <- c(out, climout)
  return(out)
}
