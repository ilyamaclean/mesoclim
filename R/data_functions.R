#' Check and summarise coarse resolution climate data inputs
#'
#' @param input_list list of climate and associated variables as output by functions like 'ukcp18toclimarray'
#' @param tstep string decribing whether input data is by hour or day
#'
#' @return where possible returns inputs suitable for downsclaing and/or warning messages if input data has unexpected values
#' Also prints summary statistics and graphs of each variable
#' @import magrittr
#' @import fmsb
#' @export
#'
#' @examples
#' climdata<- read_climdata(system.file('data/era5input.rds',package='mesoclim'))
#' chk_climdata<- checkinputs(climdata, tstep = "hour")
checkinputs <- function(input_list, tstep = c("hour","day")){
  tstep<-match.arg(tstep)

  if('climarray' %in% names(input_list)) input_list[names(input_list$climarray)]<-input_list$climarray[names(input_list$climarray)]

  check.names<-function(nms,char) {
    sel<-which(nms==char)
    if (length(sel) == 0) stop(paste0("Cannot find ",char," in weather"))
  }
  check.vals<-function(x,mn,mx,char,unit) {
    sel<-which(is.na(x))
    if (length(sel)>0) stop(paste0("Missing values in weather$",char))
    sel<-which(x<mn)
    if (length(sel)>0) {
      mx<-min(x)
      stop(paste0(mx," outside range of typical ",char," values. Units should be ",unit))
    }
    sel<-which(x>mx)
    if (length(sel)>0) {
      mx<-max(x)
      stop(paste0(mx," outside range of typical ",char," values. Units should be ",unit))
    }
  }
  check.mean<-function(x,mn,mx,char,unit) {
    me<-mean(x,na.rm=T)
    if (me<mn | me>mx) stop(paste0("Mean ",char," of ",me," implausible. Units should be ",unit))
  }
  up.lim<-function(x,mx,char) {
    mxx<-max(x,na.rm=T)
    x[x>mx]<-mx
    if (mxx>mx) warning(paste0(char," values capped at ",mx))
    x
  }
  check.unpack<-function(r,nme="r") {
    if (class(r)[1] == "PackedSpatRaster" & class(r)[1] == "SpatRaster") {
      stop(paste0(nme," must be a SpatRaster produced using the terra package"))
    }
    if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
    r
  }

  nms<-names(input_list)
  if(tstep=='hour') check.names(nms,"temp")
  if(tstep=='day'){
    check.names(nms,"tmax")
    check.names(nms,"tmin")
  }
  check.names(nms,"relhum")
  check.names(nms,"pres")
  check.names(nms,"prec")
  check.names(nms,"lwrad")
  check.names(nms,"swrad")
  check.names(nms,"windspeed")
  check.names(nms,"winddir")
  check.names(nms,"windheight_m")
  check.names(nms,"tempheight_m")
  check.names(nms,"dtm")
  check.names(nms,"tme")
  # Optional cloud cover
  if(!"cloud" %in% nms) warning('No optional cloud cover variable found - continuing...')

  # Check dtm unpacked and with crs
  input_list$dtm<-check.unpack(input_list$dtm,"dtm")
  if (is.na(crs(input_list$dtm))) stop("dtm must have a coordinate reference system specified")

  # check time and timezone
  sel<-which(is.na(input_list$tme))
  if (length(sel)>0) stop("Cannot recognise all obs_time in weather")

  if(class(input_list$tme)[1]!='POSIXlt') stop("Class of input tme variable is not POSIXlt!!!")

  if (attr(input_list$tme,"tzone") != "UTC") stop("timezone of obs_time in weather must be UTC")

  # Calculate pressure limits based on elevation
  mxelev<-max(.is(input_list$dtm),na.rm=T)
  mnelev<-min(.is(input_list$dtm),na.rm=T)
  mxp<-108.5*((293-0.0065*mnelev)/293)^5.26
  mnp<-87*((293-0.0065*mnelev)/293)^5.26

  # Check weather values
  if(tstep=='hour') check.vals(input_list$temp,-50,65,"temperature","deg C")
  if(tstep=='day'){
    check.vals(input_list$tmax,-50,65,"max temperature","deg C")
    check.vals(input_list$tmin,-50,65,"min temperature","deg C")
  }
  input_list$relhum<-up.lim(input_list$relhum,100,"relative humidity")
  check.vals(input_list$relhum,0,100,"relative humidity","percentage (0-100)")
  check.mean(input_list$relhum,5,100,"relative humidity","percentage (0-100)")
  check.vals(input_list$pres,mnp,mxp,"pressure","kPa ~101.3")
  check.vals(input_list$swrad,0,1350,"shortwave downward radiation","W / m^2")
  check.vals(input_list$lwrad,0,1350,"longwave downward radiation","W / m^2")
  check.vals(input_list$windspeed,0,100,"wind speed","m/s")

  # Wind direction
  mn<-min(input_list$winddir)
  mx<-max(input_list$winddir)
  if (mn<0 | mx>360) {
    input_list$winddir<-input_list$winddir%%360
    warning("wind direction adjusted to range 0-360 using modulo operation")
  }

  # Calculate clear-sky radiation to compare with downward SW down
  #ll<-.latlongfromraster(input_list$dtm)
  #if(tstep=='day') tmean<-(input_list$tmax+input_list$tmin)/2 else tmean<-input_list$temp
  #csr<-.clearskyrad(input_list$tme,ll$lat,ll$long,tmean,input_list$relhum,input_list$pres)
  #csd<-csr+50
  #sel<-which(input_list$swrad>csd)
  #if (length(sel)>0) {
  #  warning("Short wave radiation values significantly higher than expected clear-sky radiation values")
  #}

  # Calculate data entries and time period
  obs_int<-lubridate::int_diff(input_list$tme)
  int <- interval(input_list$tme[1], input_list$tme[length(input_list$tme)])
  data_days<-time_length(int, "days")
  data_months<-round(time_length(int, "months"),2)
  data_yrs<-round(time_length(int, "years"),2)
  print(paste('Weather observations =',length(input_list$tme)))
  print(paste('Timesteps=',mean(obs_int)/3600, 'hrs, max=',max(obs_int)/3600,'hrs, min=', min(obs_int)/3600,'hrs'))
  print(paste('Observations over',data_yrs,'years, or',data_months,'months, or',data_days,'days.'))

  # Check other variables
  if (tstep == "hour") {
    if(length(input_list$tme)%%24 != 0) stop ("Hourly weather data needs to include data for entire days (24 hours)")
  }

  # Print summary stats and figures of input variables?
  statvars<-nms[!nms %in% c('dtm','climarray','tme','windheight_m','tempheight_m')]
  elev<-data.frame(elevation=values(input_list$dtm,mat=FALSE))
  stats_df<-as.data.frame(t(round(sapply(input_list[statvars],summary),3)))[,c('Min.','Mean','Max.')]
  elev_stat<-round(summary(elev$elevation),3)[c('Min.','Mean','Max.')]
  stats_df<-rbind(stats_df, 'elevation'=elev_stat)
  print(stats_df)

  # Plot spatial variation by different time steps using min, max, mean
  if(tstep=='day') vars<-c('relhum','pres','prec','lwrad','swrad','tmax','tmin','windspeed','winddir')
  if(tstep=='hour') vars<-c('relhum','pres','prec','lwrad','swrad','temp','windspeed','winddir')
  if('cloud' %in% nms) vars<-c('cloud',vars)
  if(data_yrs>=6) {
    print('Plotting spatial variation by year: red=max, green=mean, blue=min')
    ncol<-floor(sqrt(length(vars)))
    if(length(vars)%%ncol!=0) nrow<-ncol+1 else nrow<-ncol
    par(mar=c(1,1,1,1))
    par(mfrow=c(nrow,ncol))
    for(v in vars){
      r<-.rast(input_list[[v]],input_list$dtm)
      terra::time(r)<-input_list$tme
      plot_timestats_r(r,var,idx='years')
    }  }
  if(data_yrs<6 && data_months>=6) {
    print('Plotting spatial variation by month: red=max, green=mean, blue=min')
    ncol<-floor(sqrt(length(vars)))
    if(length(vars)%%ncol!=0) nrow<-ncol+1 else nrow<-ncol
    par(mar=c(1,1,1,1))
    par(mfrow=c(nrow,ncol))
    for(v in vars){
      r<-.rast(input_list[[v]],input_list$dtm)
      terra::time(r)<-input_list$tme
      plot_timestats_r(r,v,idx='months')
    }  }
  if(data_months<6) {
    print('Plotting spatial variation by day of year: red=max, green=mean, blue=min')
    ncol<-floor(sqrt(length(vars)))
    if(length(vars)%%ncol!=0) nrow<-ncol+1 else nrow<-ncol
    par(mar=c(1,1,1,1))
    par(mfrow=c(nrow,ncol))
    for(v in vars){
      r<-.rast(input_list[[v]],input_list$dtm)
      terra::time(r)<-input_list$tme
      plot_timestats_r(r,v,idx='doy')
    }  }

  # Plot wind direction
  print('Plotting wind direction figures')
  plot_wind(input_list$winddir,input_list$windspeed)

  return(input_list)
}

#' @title Read climate data
#' @description
#' Reads a climate dataset as output by for example `ukcp18toclimdata()` and written by `write_climdata()`
#' Unwraps any PackedSpatRasters to Spatrasters
#' @param filepath
#'
#' @return R list of climate data
#' @export
#'
#' @examples
#' climdata<-read_climdata(system.file('data/ukcp18rcm.rda',package='mesoclim'))
read_climdata<-function(filepath){
  if(file.exists(filepath)!=TRUE) stop('Filepath provided does NOT exist!!') else{
    climdata<-readRDS(filepath)
    climdata<-lapply(climdata,function(x) if(class(x)[1]=='PackedSpatRaster') terra::unwrap(x) else x)
  }
  return(climdata)
}

#' @title Write climate data
#' @description
#' Writes to a binary file a climate dataset (as output by for example `ukcp18toclimdata()`.
#'Wraps any Spatrasters to PackedSpatRasters before writing.
#' @param climdata a list of climate variables as output by `ukcp18toclimdata()`
#' @param filepath to be written
#'
#' @return writes an external file
#' @export
#'
#' @examples
#' dirct<-tempdir()
#' write_climdata(ukcp18rcm,file.path(dirct,'filename.rds'))
write_climdata<-function(climdata,filepath,overwrite=FALSE){
  if(file.exists(filepath)==TRUE & overwrite==FALSE) stop('Existing file of that name - requires setting overwrite to TRUE !!')
  else{
    climdata<-lapply(climdata,function(x) if(class(x)[1]=='SpatRaster') terra::wrap(x) else x)
    saveRDS(climdata,filepath)
  }
  return()
}
