
#' Download ERA5 reanalysis data
#'
#' @param dir_out - Local output directory
#' @param file_out - Output filename
#' @param year - year of data to retrieve
#' @param area - bounding box with order c(ymax,xmin,ymin,xmax)
#' @param variables - Long names of variables to retrieve
#' @param era5_user - cds user name/email
#' @param era5_key - cds api key
#'
#' @return filepath of output
#' @details Downloads a whole single year of reanalysis-era5-single-levels as netcdf file.
#' Variables chosen and size of area defaults match UK/Eire and variables required for microclimate modelling.
#' Limits on download file may restrict size of area and/or number of variables that can be downloaded to single file.
#' @import curl
#' @import ecmwfr
#' @export
#'
#' @examples
#' \dontrun{
#' dir_out <- tempdir()
#' cds_user<-"your_user_name"
#' cds_key<-"your_key_string"
#' download_era5(dir_out,file_out='test.nc',year=2020,area = c(55,-3,49,3),variables =c("2m temperature"),era5_user=cds_user,era5_key=cds_key)
#' }
download_era5<-function(dir_out,
                        file_out,
                        year,
                        area = c(61, -11.25, 48.75, 2.5),
                        variables = c("10m u-component of wind",
                                      "10m v-component of wind",
                                      "2m dewpoint temperature",
                                      "2m temperature",
                                      "Mean sea level pressure",
                                      "Surface pressure",
                                      "Total precipitation",

                                      "Surface solar radiation downwards",
                                      "Total sky direct solar radiation at surface",

                                      "Mean surface downward short-wave radiation flux",
                                      "Mean surface direct short-wave radiation flux",
                                      "Mean surface downward long-wave radiation flux",
                                      "Mean surface net long-wave radiation flux" ),
                        era5_user,
                        era5_key
){

  wf_set_key(key=era5_key,user= era5_user)

  request <- list(
    dataset_short_name = "reanalysis-era5-single-levels",
    product_type = "reanalysis",
    variable = variables,
    year = year,
    month =c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'),
    day = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',
            '25', '26', '27', '28', '29', '30', '31'),
    time =c('00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'),
    data_format = "netcdf",
    download_format = "unarchived",
    area = area,
    target = file_out
  )
  file<-wf_request(request = request,
                   user = era5_user,
                   transfer = TRUE,
                   path = dir_out)
  return(file)
}

#' Download ancillary era5 data (geopotential and lansea mask)
#'
#' @param dir_out - directory to download files to
#' @param area - lat lon of area
#' @param variables one or both of geopotential and land_sea_mask
#' @param era5_user CDS user name
#' @param era5_key CDS API token
#'
#' @return filepaths of outputs
#' @details Downloads a single instance of ancillary dat for reanalysis-era5-single-levels as netcdf files.
#' @import curl
#' @import ecmwfr
#' @export
#'
#' @examples
#' \dontrun{
#' dir_out <- tempdir()
#' cds_user<-"your_user_name"
#' cds_key<-"your_key_string"
#' download_ancillary_era5(dir_out,area = c(55,-3,49,3),variables =c("geopotential","land_sea_mask" ),era5_user=cds_user,era5_key=cds_key)
#' }
download_ancillary_era5<-function(dir_out,
                                  area = c(61, -11.25, 48.75, 2.5),
                                  variables = c("geopotential",
                                                "land_sea_mask" ),
                                  era5_user,
                                  era5_key
){

  wf_set_key(key=era5_key,user= era5_user)
  for(v in variables){
    request <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type = "reanalysis",
      variable = v,
      year = 2020,
      month =c('01'),
      day = c('01'),
      time =c('00:00'),
      data_format = "netcdf",
      download_format = "unarchived",
      area = area,
      target = paste0("era5_",v,".nc")
    )
    file<-wf_request(request = request,
                     user = era5_user,
                     transfer = TRUE,
                     path = dir_out)
  }
  return(file)
}
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
