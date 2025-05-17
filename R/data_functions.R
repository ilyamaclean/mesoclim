#' @title Plot timeseries plots of spatial mean, min and max values by various timesteps
#'
#' @param r stack of spatrasters with time values
#' @param var variable name used in plot titles
#' @param idx time index indicating whether summary plots by year, month, week, day of year etc
#' @param lgd whether to include legend
#'
#' @return NA plots series of timeploat
#' @export
#' @import magrittr
#' @importFrom graphics legend matplot
#' @keywords graph
plot_timestats_r<-function(r,v,idx=c('years', 'months', 'week',  'doy', 'yearmonths', 'yearweeks', '7days','hour'),lgd=FALSE){
  if(idx=='hour'){
    time_mean<- global(r,"mean", na.rm=TRUE)
    time_max<-global(r, "max", na.rm=TRUE)
    time_min<-global(r, "min", na.rm=TRUE)
    plot_df<-as.data.frame(cbind(tstep=terra::time(r),mean=time_mean$mean,max=time_max$max,min=time_min$min))
    plot_df<-plot_df[order(plot_df$tstep),]
  } else{
    time_mean<-tapp(r,index=idx,fun=mean) %>% global("mean", na.rm=TRUE)
    time_max<-tapp(r,index=idx,fun=max) %>% global("max", na.rm=TRUE)
    time_min<-tapp(r,index=idx,fun=min) %>% global("min", na.rm=TRUE)
    plot_df<-as.data.frame(cbind(tstep=as.numeric(sapply(strsplit(rownames(time_mean),'_'),tail,1)),mean=time_mean$mean,max=time_max$max,min=time_min$min))
    plot_df<-plot_df[order(plot_df$tstep),]
  }
  matplot(as_datetime(plot_df$tstep), plot_df[,2:4], type = "l", lty = 1,
          col = c("green", "red", "blue"), xlab = idx, ylab = v, font.main = 1,
          tck = 0.02, cex.main=1, cex.axis=0.7, main = paste(v,'by',idx), cex.main=1)
  if(lgd==TRUE) legend("topright", legend = c("Mean", "Max", "Min"), cex=0.5,
                       col = c("green", "red", "blue"),
                       lty = 1)
}

#' @title direction to cardinal label
#' @noRd
.dir_to_cardinal <- function(x) {
  upper <- seq(from = 11.25, by = 22.5, length.out = 17)
  cardinals <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  ifelse(x>360 | x<0,NA,cardinals[findInterval(x,upper,rightmost.closed = T)+1])
}

#' @title Plot star or radar chart of wind cardinal directions
#' @param winddir wind direction 3D Spatraster
#' @param windspeed windspeed 3D Spatraster of same geometry as `winddir`
#' @return plots star charts showing frequency of wind directions and of
#' average and maximum windspeeds by direction
#' @export
#' @import fmsb
#' @importFrom graphics legend par
#' @keywords graph
#' @examples
#' wdir<-c(357,5,20,45,100,185,265,275,290)
#' wspeed<-seq(3,27,3)
#' plot_wind(wdir,wspeed)
plot_wind<-function(winddir,windspeed){
  # Classify by cardinal group
  dgrp<- .dir_to_cardinal(c(winddir))
  dgrp_mean<-tapply(c(windspeed),dgrp,mean)
  dgrp_max<-tapply(c(windspeed),dgrp,max)
  dgrp_freq<-as.array(table(dgrp))

  # Mean/Max windspeed plot by dir
  cardorder<- c('N', 'NNW', 'NW', 'WNW', 'W', 'WSW', 'SW', 'SSW', 'S', 'SSE', 'SE', 'ESE', 'E', 'ENE', 'NE', 'NNE')
  df<-as.data.frame(matrix(0, 4, 16) )
  names(df)<-cardorder
  rownames(df)<-c('Max','Min','Mean windspeed','Max windspeed')
  df[1,]<-rep(max(dgrp_max),length(cardorder)) # max speed across  all directions
  df[2,]<-rep(0,length(cardorder)) # min speed across all dir
  df[3,names(dgrp_mean)]<-dgrp_mean
  df[4,names(dgrp_max)]<-dgrp_max

  par(mar=c(1,1,1,1))
  par(mfrow=c(1,2))
  fmsb::radarchart(df, title='Mean and max windspeeds', cex=0.75)
  legend(
    x = "bottomleft", legend = c('Mean','Maximum'), horiz = TRUE,
    bty = "n", pch = 20 , col = c("black", "red"),
    text.col = "black", cex = 1, pt.cex = 1
  )

  # Freq plot by dir
  df<-as.data.frame(matrix(0, 3, 16) )
  names(df)<-cardorder
  rownames(df)<-c('Max','Min','Frequency')
  df[1,]<-rep(max(dgrp_freq),length(cardorder)) # max speed across  all directions
  df[2,]<-rep(0,length(cardorder)) # min speed across all dir
  df[3,names(dgrp_freq)]<-dgrp_freq

  fmsb::radarchart(df, title='Wind direction frequency',pcol=c('green'), cex=0.75)

  legend(
    x = "bottomleft", legend = c('Frequency'), horiz = TRUE,
    bty = "n", pch = 20 , col = c("green"),
    text.col = "black", cex = 1, pt.cex = 1
  )
}

#' @title Check and summarise coarse resolution climate data inputs
#' @param input_list list of climate and associated variables as output by functions like 'ukcp18toclimarray'
#' @param tstep string decribing whether input data is by hour or day

#' @return where possible returns inputs suitable for downsclaing and/or warning messages if input data has unexpected values
#' Also prints summary statistics and graphs of each variable
#' @import magrittr
#' @import fmsb
#' @importFrom graphics par
#' @importFrom stats sd
#' @export
#' @keywords preprocess
#' @examples
#' climdata<- read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' chk_climdata<- checkinputs(climdata, tstep = "day")
checkinputs <- function(input_list, tstep = c("hour","day"),plots=TRUE){
  tstep<-match.arg(tstep)

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
  min_timestep<-min(obs_int)/3600
  max_timestep<-max(obs_int)/3600
  if(min_timestep==max_timestep) print(paste('Timesteps=',min_timestep,'hours for whole timeseries.')) else warning(paste("Different timesteps varying from",min_teimstep,"to",max_timestep,"- CHECK DATA !!!"))
  print(paste('Observations from',input_list$tme[1], 'to',input_list$tme[length(input_list$tme)]))
  print(paste('Temperature values are at',input_list$tempheight_m,'metres above ground level.'))
  print(paste('Wind speeds are at',input_list$windheight_m,'metres above ground level.'))
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
  if(plots){
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
        plot_timestats_r(r,v,idx='years')
      }  }
    if(data_months>=6) {
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
  }
  return(input_list)
}

#' @title Read climate data
#' @description
#' Reads a climate dataset as output by for example `ukcp18toclimdata()` and written by `write_climdata()`
#' Unwraps any PackedSpatRasters to Spatrasters
#' @param filepath - either a pre-loaded list or filepath to a list of prepared climate data as written by `write_climdata()`
#'
#' @return R list of climate data
#' @export
#' @keywords preprocess data
#' @examples
#' climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' climdata<-read_climdata(ukcpinput)
read_climdata<-function(filepath){
  if(class(filepath)=="list") climdata<-lapply(filepath,function(x) if(class(x)[1]=='PackedSpatRaster') terra::unwrap(x) else x)
  if(class(filepath)=="character"){
    if(file.exists(filepath)!=TRUE) stop('Filepath provided does NOT exist!!')
    typ<-substr(filepath,(nchar(filepath)-3),nchar(filepath))
    if (typ==".rda"){
      load(filepath)
      filename<-rev(unlist(strsplit(filepath,'/')))[1]
      varname<-substr(filename,1,nchar(filename)-4)
      climdata<-mget(varname)[[1]]
    }  else climdata<-readRDS(filepath)
    climdata<-lapply(climdata,function(x) if(class(x)[1]=='PackedSpatRaster') terra::unwrap(x) else x)
  }
  return(climdata)
}


#' @title Write climate data
#' @description
#' Writes to a binary file a climate dataset (as output by for example `ukcp18toclimdata()`.
#' Wraps any Spatrasters to PackedSpatRasters before writing.
#' @param climdata a list of climate variables as output by `ukcp18toclimdata()`
#' @param filepath to be written
#' @return writes an external file
#' @export
#' @keywords preprocess data
#' @examples
#' climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' dir_temp<-tempdir()
#' write_climdata(climdata,file.path(dir_temp,'filename.rds'))
write_climdata<-function(climdata,filepath,overwrite=FALSE){
  if(file.exists(filepath)==TRUE & overwrite==FALSE) stop('Existing file of that name - requires setting overwrite to TRUE !!')
  typ<-substr(filepath,(nchar(filepath)-3),nchar(filepath))
  if(!typ %in% c(".rda",".rds",".Rds")) stop('File typ not recognized !!')

  climdata<-lapply(climdata,function(x) if(class(x)[1]=='SpatRaster') terra::wrap(x) else x)
  if(typ==".rda") save(climdata,file=filepath) else saveRDS(climdata,filepath)
  return()
}


