#' @title Downscale daily temperatures
#' @description Downscales daily coarse-resolution temperature stats accounting for elevation
#' effects and optionally cold air drainage and coastal effects
#' @param climdata a `climdata` model object containing climate data of the same format as `era5climdata`
#' Expects variables to include a `tmin` and `tmax`
#' @param tmean a coarse resolution array or spatRaster of daily mean temperatures identical to climdata inputs.
#' If NA will calculate from hourluy downscaling of tmin and tmax in `climdata`
#' @param sst a SpatRast of sea-surface temperature data (deg C) timeseries that overlaps climdata$tme
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param uzf optionally, a SpatRast of fine resolution wind speeds as returned by [winddownscale()].
#' Calculated if not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param tempvar string name of element of climdata holding temperature data (ie 'temp','tmax','tmin')
#' @param thgto height above ground of output temperature measurements.
#' @param whgto height above ground of output wind speed measurements.
#' to resolution of `dtmf`.
#' @return a list of three multi-layer SpatRasters of downscaled daily minimum, maximum and mean temperatures (deg C) matching the
#' resolution of dtmf.
#' @details Cold air drainage is calculated by delineating hydrological basins and
#' calculating flow accumulation and the elevation difference from the highest point of the basin.
#' Cold-air drainage is assumed to occur when atmospheric stability is high, namely when
#' the radiation balance is negative and wind speeds are low. Coastal effects are
#' calculated by determining the ratio of land to sea pixels in an upwind direction.
#' The provision of `dtmm` allows this ratio to be derived accoutning for land and
#' sea outside the boundaries of the study area.
#' @rdname tempdaily_downscale
#' @import terra
#' @export
#' @keywords spatial
#' @examples
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' dtmf<-terra::rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))
#' dtmm<-terra::rast(system.file("extdata/dtms/dtmm.tif",package="mesoclim"))
#' basins<-basindelin(dtmf, boundary = 2)
#' wsf<-winddownscale(climdata$windspeed, climdata$winddir, dtmf, dtmm, climdata$dtm, zi=climdata$windheight_m)
#' dailytemps<-tempdaily_downscale(climdata,NA,terra::unwrap(mesoclim::ukcp18sst),dtmf,dtmm,basins,wsf,cad = TRUE,coastal = TRUE,2,2)
#' terra::panel(c(dailytemps$tmin[[13]],dailytemps$tmax[[13]],dailytemps$tmean[[13]]),main=paste(c("Tmin","Tmax","Tmean"),"13/05/2018"))
tempdaily_downscale<-function(climdata,tmean=NA,sst=NA,dtmf,dtmm,basins,uzf,cad=TRUE,coastal=TRUE,thgto=2,whgto=2){
  # Convert variables - unpack any wrapped spatRasters and convert arrays to spatraster
  input_class<-lapply(lapply(climdata,class),"[",1)
  if(any(input_class=="PackedSpatRaster")) climdata[which(input_class=="PackedSpatRaster")]<-lapply(climdata[which(input_class=="PackedSpatRaster")],unwrap)
  if(any(input_class=="array")) climdata[which(input_class=="array")]<-lapply(climdata[which(input_class=="array")],.rast,tem=climdata$dtm)
  if(inherits(sst,"PackedSpatRaster")) sst<-unwrap(sst)

  # Test to see if any valid sst data provided to estimate coastal effects
  if(inherits(sst,"SpatRaster")) if(all(is.na(values(sst[[1]])))) sst<-NA
  if(is.logical(sst)){
    if(coastal==TRUE) message("Ignoring coastal effects as no valid sea surface temperature data provided - assuming inland area!!")
    coastal<-FALSE
  }

  # Check tmin and tmax among climdata variables
  if(!any(c("tmin","tmax") %in% names(climdata))) stop("Cannot find daily min/max temperature variables for tempdaily_downscale!!!")
  if (class(dtmm) == "logical" & coastal) stop("dtmm needed for calculating coastal effects")
  dtmc<-climdata$dtm

  # Get variables
  rh<-climdata$relhum
  pk<-climdata$pres
  tmin<-climdata$tmin
  tmax<-climdata$tmax
  thgti<-climdata$tempheight_m
  whgti<-climdata$windheight_m
  tme<-climdata$tme
  u2<-.windhgt(climdata$windspeed, whgti, thgto) # wind at downscale temperature height above ground (for cad)
  swrad<-climdata$swrad
  lwrad<-climdata$lwrad
  dtmc<-climdata$dtm

  # Get temp independent params
  if (cad) mu<-.cad_multiplier(dtmf, basins, refhgt=thgti)
  if (coastal) {
    # Interpolate sst to all landcells and across full timeseries so no NA
    if (any(global(sst,anyNA))) sstinterp<-.spatinterp(sst) else sstinterp<-sst
    sstinterp<-.tmeinterp(sstinterp,NA,tme)
    if (crs(sst) != crs(dtmf)) sstinterp<-project(sstinterp,crs(dtmf))
    sstf<-.resample(sstinterp,dtmf,method="cubic")
    # Calc windspeed at output height if required
    #if(class(uzf)[1] == "logical") uzf<-winddownscale(climdata$windspeed,climdata$winddir,dtmf,dtmm,dtmc,wca,whgti,thgto)
  }

  # Downscale tmin
  lrc<-lapserate(.is(tmin), .is(rh), .is(pk))
  lrc<-.rast(lrc,dtmc)
  if (crs(lrc) != crs(dtmf)) {
    lrcp<-project(lrc,crs(dtmf))
    lrf<-.resample(lrcp,dtmf)
  } else lrf<-.resample(lrc,dtmf)
  tminf<-.tempelev(tmin,lrf,lrc,dtmf+thgto,ifel(is.na(dtmc),0,dtmc+thgti))
  if (cad) {
    st<-.cad_conditions(tmin,u2,swrad=0,lwrad,dtmc,dtmf,refhgt=thgti) #swrad=0 as we assume tmin in night
    tcad<-.apply_cad(lrf,mu,st)
    tminf<-tminf+tcad
  }
  if (coastal) tminf<-.tempcoastal(tc=tminf,sst=sstf,u2=uzf,wdir=climdata$winddir,dtmf,dtmm,dtmc)

  # Downscale tmax - no cold air drainage
  lrc<-lapserate(.is(tmax), .is(rh), .is(pk))
  lrc<-.rast(lrc,dtmc)
  if (crs(lrc) != crs(dtmf)) {
    lrcp<-project(lrc,crs(dtmf))
    lrf<-.resample(lrcp,dtmf)
  } else lrf<-.resample(lrc,dtmf)
  tmaxf<-.tempelev(tmax,lrf,lrc,dtmf+thgto,ifel(is.na(dtmc),0,dtmc+thgti))
  if (coastal) tmaxf<-.tempcoastal(tc=tmaxf,sst=sstf,u2=uzf,wdir=climdata$winddir,dtmf,dtmm,dtmc)

  # Check tax>tmin (coastal effects can switch)
  diurnal<-(tmaxf-tminf)<0
  if(unique(all(diurnal))!=0){
    temp_tmin<-ifel(diurnal,tmaxf,tminf)
    tmaxf<-ifel(diurnal,tminf,tmaxf)
    tminf<-temp_tmin
    remove(temp_tmin)
    remove(diurnal)
  }

  # Downscale daily mean (calculated from hourly downscaling of climate tmin and tmax)
  # Assumes no cold air drainage
  if(is.logical(tmean)) tmean<-.hourtoday(temp_dailytohourly(tmin, tmax, climdata$tme),mean)
  lrc<-lapserate(.is(tmean), .is(rh), .is(pk))
  lrc<-.rast(lrc,dtmc)
  if (crs(lrc) != crs(dtmf)) {
    lrcp<-project(lrc,crs(dtmf))
    lrf<-.resample(lrcp,dtmf)
  } else lrf<-.resample(lrc,dtmf)
  tmeanf<-.tempelev(tmean,lrf,lrc,dtmf+thgto,ifel(is.na(dtmc),0,dtmc+thgti))
  if (coastal) tmeanf<-.tempcoastal(tc=tmeanf,sst=sstf,u2=uzf,wdir=climdata$winddir,dtmf,dtmm,dtmc)

  #tmeanf<-.hourtoday(temp_dailytohourly(tminf, tmaxf, climdata$tme),mean) # alternative but very slow!!!

  # Format outputs
  terra::time(tminf)<-tme
  terra::time(tmaxf)<-tme
  terra::time(tmeanf)<-tme

  out<-list(tmin=tminf,tmax=tmaxf,tmean=tmeanf)
  return(out)
}
#' @title Downscale hourly temperature
#' @description Downscales coarse-resolution temperature data accounting for elevation
#' effects and optionally cold air drainage and coastal effects
#' @param climdata a `climdata` model object containing climate data of the same format as `era5climdata`.
#' Expects temperature to be in degrees and with the name `temp`.
#' @param sst a SpatRast of sea-surface temperature data (deg C) timeseries that overlaps climdata$tme
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param uzf optionally, a SpatRast of fine resolution wind speeds as returned by [winddownscale()].
#' Calculated if not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param thgto height above ground of output temperature measurements.
#' @param whgto height above ground of output wind speed measurements.
#' @param tempvar string name of element of climdata holding temperature data (default = 'temp')
#' @return a multi-layer SpatRast of downscaled temperatures (deg C) matching the
#' resolution of dtmf.
#' @details Cold air drainage is calculated by delineating hydrological basins and
#' calculating flow accumulation and the elevation difference from the highest point of the basin.
#' Cold-air drainage is assumed to occur when atmospheric stability is high, namely when
#' the radiation balance is negative and wind speeds are low. Coastal effects are
#' calculated by determining the ratio of land to sea pixels in an upwind direction.
#' The provision of `dtmm` allows this ratio to be derived accoutning for land and
#' sea outside the boundaries of the study area.
#' @rdname temphrly_downscale
#' @import terra
#' @export
#' @keywords spatial
#' @examples
#' climdata<- read_climdata(mesoclim::ukcpinput)
#' climhrly<-temporaldownscale(climdaily, adjust = TRUE, clearsky=NA, srte = 0.09, relmin = 10, noraincut = 0)
#' dtmf<-terra::rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))
#' dtmm<-terra::rast(system.file("extdata/dtms/dtmm.tif",package="mesoclim"))
#' basins<-basindelin(dtmf, boundary = 2)
#' wsfhr<-winddownscale(climhrly$windspeed, climhrly$winddir, dtmf, dtmm, climhrly$dtm, zi=climhrly$windheight_m)
#' sst<-terra::unwrap(ukcp18sst)
#' tempf<-temphrly_downscale(climhrly, sst, dtmf, dtmm, basins, uzf = wsfhr,tempvar="temp")
#' terra::panel(tempf[[c(1,6,12,18)]],main=paste('Temperature',c('00:00', '06:00', '12:00', '18:00'),' on 1 May 2018'))
#' matplot(terra::time(tempf),unlist(terra::global(tempf,mean,na.rm=TRUE)), type = "l", lty = 1, main='Mean hourly temperature of whole area')
temphrly_downscale<-function(climhrly, sst, dtmf, dtmm = NA, basins = NA, uzf = NA,
                        cad = TRUE, coastal = TRUE,thgto=2, whgto=2, tempvar='temp') {
  # Convert variables - unpack any wrapped spatRasters and convert arrays to spatraster
  input_class<-lapply(lapply(climhrly,class),"[",1)
  if(any(input_class=="PackedSpatRaster")) climhrly[which(input_class=="PackedSpatRaster")]<-lapply(climhrly[which(input_class=="PackedSpatRaster")],unwrap)
  if(any(input_class=="array")) climhrly[which(input_class=="array")]<-lapply(climhrly[which(input_class=="array")],.rast,tem=climhrly$dtm)
  if(inherits(sst,"PackedSpatRaster")) sst<-unwrap(sst)

  # Test to see if any valid sst data provided to estimate coastal effects
  if(inherits(sst,"SpatRaster")) if(all(is.na(values(sst[[1]])))) sst<-NA
  if(is.logical(sst)){
    if(coastal==TRUE) message("Ignoring coastal effects as no valid sea surface temperature data provided - assuming inland area!!")
    coastal<-FALSE
  }

  dtmc<-climhrly$dtm
  if(is.logical(sst)) coastal<-FALSE else if(all(!is.na(values(sst[[1]])))) coastal<-FALSE
  rh<-climhrly$relhum
  pk<-climhrly$pres
  tc<-climhrly[[tempvar]]
  thgti<-climhrly$tempheight_m
  whgti<-climhrly$windheight_m
  tme<-climhrly$tme
  u2<-.windhgt(climhrly$windspeed, whgti, thgto) # wind at downscale temperature height above ground (for cad)
  swrad<-climhrly$swrad
  lwrad<-climhrly$lwrad
  dtmc<-climhrly$dtm

  # Calculate lapse rates (coarse and fine scale)
  lrc<-lapserate(.is(tc), .is(rh), .is(pk))
  lrc<-.rast(lrc,dtmc)
  if (crs(lrc) != crs(dtmf)) {
    lrcp<-project(lrc,crs(dtmf))
    lrf<-.resample(lrcp,dtmf)
  } else lrf<-.resample(lrc,dtmf)

  # Calculate elevation effects (assume NA in dtmc =sea = height 0) - add input / output heights to elevations
  tcf<-.tempelev(tc,lrf,lrc,dtmf+thgto,ifel(is.na(dtmc),0,dtmc+thgti))

  # Cold air drainage - assume if applied to tmin occurs at nighttime when swrad=0 - is it valid to apply to tmax?
  if (cad) {
    mu<-.cad_multiplier(dtmf, basins, refhgt=thgti)
    st<-.cad_conditions(tc,u2,swrad=0,lwrad,dtmc,dtmf,refhgt=thgti) #swrad=0 as we assume tmin in night
    tcad<-.apply_cad(lrf,mu,st)
    tcf<-tcf+tcad
  }

  # Coastal effects - checks if non NA values in sst AND resamples to res/extent of dtmf
  if (coastal) {
    # Interpolate sst to all landcells and across full timeseries so no NA
    if (any(global(sst,anyNA))) sstinterp<-.spatinterp(sst) else sstinterp<-sst
    sstinterp<-.tmeinterp(sstinterp,NA,tme)
    if (crs(sst) != crs(dtmf)) sstinterp<-project(sstinterp,crs(dtmf))
    sstf<-.resample(sstinterp,dtmf,method="cubic")
    # Calc windspeed at output height if required
    if(class(uzf)[1] == "logical") uzf<-winddownscale(climhrly$windspeed,climhrly$winddir,dtmf,dtmm,dtmc,whgti,thgto)
    tcf<-.tempcoastal(tc=tcf,sst=sstf,u2=uzf,wdir=climhrly$winddir,dtmf,dtmm,dtmc)
  }
  hr_tme<-as.POSIXlt(unlist(lapply(tme,FUN=function(x) x+(60*60*c(0:23)) )))
  terra::time(tcf)<-hr_tme
  return(tcf)
}
#' @title Downscale pressure with elevation effects
#' @description Downscales an array of pressure data applying elevation effects
#' @param pk a SpatRast or array of atmospheric pressures (kPa). If an array dtmc must
#' be provided.
#' @param dtmf a fine-resolution SpatRast of elevations. Temperatures down-scaled
#' to resolution of `dtmf`.
#' @param dtmc optional SpatRast of elevations matching resolution of `pk`. If
#' not supplied, `pk` must be a SpatRast and `dtmc` is derived by resampling `dtmf`
#' to resolution of `pk` and extents must match.
#' @param sealevel optional logical indicating whether `pk` is pressure at sea-level
#' @return a multi-layer SpatRast of elevation-corrected pressures (kPa) matching the
#' resolution of dtmf.
#' @rdname presdownscale
#' @import terra
#' @export
#' @keywords spatial
#' @examples
#' climdata<- read_climdata(mesoclim::ukcpinput)
#' pk <- presdownscale(climdata$pres, terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim')), climdata$dtm)
#' terra::panel(c(pk[[1]],pk[[5]],pk[[10]],pk[[20]]),main=paste0("Atmospheric pressure ",c(1,5,10,15),"/05/2018"))
presdownscale<-function(pk, dtmf, dtmc, sealevel = TRUE) {
  if (class(pk)[1] == "array"){
    pk<-.rast(pk,dtmc)
    tme<-NA } else tme<-time(pk)
  dtmc<-ifel(is.na(dtmc),0,dtmc)
  # Convert to sea level
  if (sealevel == FALSE) {
    pk<-pk/(((293-0.0065*dtmc)/293)^5.26)
  }
  # Convert to fine res elevations
  if (crs(dtmf) != crs(dtmc)) pk<-project(pk,crs(dtmf))
  pkf<-.resample(pk,dtmf)
  pkf<-pkf*((293-0.0065*dtmf)/293)^5.26
  if(!inherits(tme,"logical")) time(pkf)<-tme
  return(pkf)
}

#' @title Downscale longwave radiation
#' @description Downscales an array of longwave downward radiation data with option to
#' simulate terrain shading of skyview
#'
#' @param lwrad - coarse resolution downward longwave radiation (W/m^2)  as spatraster or array.
#' @param tc - coarse resolution temperature (deg C) as spatraster or array.
#' @param tcf - fine resolution spatraster of temperature (deg C) down-scaled to resolution of `dtmf`.
#' @param tme POSIXlt object of times corresponding to radiation values in
#' `lwrad`, `tc` and `tcf`.
#' @param dtmf a fine-resolution SpatRast of elevations.
#' @param dtmc SpatRast of elevations matching the resolution, extend and coordinate
#' reference system of `lwrad` and `tc`.
#' @param terrainshade - boolean TRUE/FALSE
#'
#' @return  a multi-layer SpatRast of downward longwave radiation (W/m^2)
#' matching the resolution of dtmf.
#' @details radiation is downscaled by computing from temperature the difference in upward
#' longwave radiation at fine and coarse resolutions. Corrected LW down is then
#' corrected for terrain shading of skyview if requested.
#' @export
#' @examples
#' climdata<- read_climdata(mesoclim::ukcpinput)
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' tmean_c<-(climdata$tmin+climdata$tmax)/2
#' wsf<-winddownscale(climdata$windspeed, climdata$winddir, dtmf, dtmm, climdata$dtm, zi=climdata$windheight_m)
#' dailytemps<-tempdaily_downscale(climdata,NA,terra::unwrap(mesoclim::ukcp18sst),dtmf,dtmm,NA,wsf,cad = FALSE,coastal = FALSE,2,2)
#' lw<-lwdownscale(climdata$lwrad,tmean_c,dailytemps$tmean, climdata$tme, dtmf=dtmf, dtmc=climdata$dtm)
#' terra::panel(c(lw[[1]],lw[[5]],lw[[10]],lw[[20]]),main=paste0("LW down ",c(1,5,10,15),"/05/2018"))
#'
lwdownscale<-function(lwrad, tc, tcf, tme, dtmf, dtmc, skyview=NA, terrainshade = TRUE) {
  if(class(lwrad)[1]=="array") lwrad<-.rast(lwrad,dtmc)
  lwf<-.resample(lwrad,dtmf, msk=TRUE)
  lwupc<-.lwup(tc)
  if(class(lwupc)[1]=="array") lwupc<-.rast(lwupc,dtmc)
  lwupc<-.resample(lwupc,dtmf)
  lwupf<-.lwup(tcf)
  lwf<-lwf+(lwupc-lwupf)
  if (terrainshade) {
    if(inherits(skyview,'logical')) skyview<-.skyview(dtmf)
    #lwf<-.rast(.is(lwf)*skyview,dtmf)
    lwf<-lwf*skyview
  }
  return(lwf)
}


#' @title Downscale shortwave radiation
#' @description Downscales an array of shortwave radiation data with option to
#' simulate cloud patchiness and calculate diffuse radiation
#' @param swrad an array of shortwave radiation (W/m^2).
#' @param tme POSIXlt object of times corresponding to radiation values in `swrad`.
#' @param dtmf a fine-resolution SpatRast of elevations. Temperatures down-scaled
#' to resolution of `dtmf`.
#' @param dtmc SpatRast of elevations matching the resolution, extend and coordinate
#' reference system of of `swrad`.
#' @param patchsim optional logical indicating whether to simulate cloud cover
#' patchiness during downscaling. More realistically captures variation, but slower.
#' @param nsim optionally the number of independent cloud cover patchiness simulations to perform. Outputs
#' are temporally interpolated (see details).
#' @param hor optional 3D array or spatRaster of horizone angles
#' @param terrainshade optional logical indicating whether to adjust incoming radiation
#' for terrain shading. If TRUE, total downward shortwave radiation is partioned into its
#' direct and diffuse components, the latter adjusted by sky view and the former set
#' to zero if the sun is below the horizon, which is computed explicitly.
#' @return if `terrainshade = FALSE` a multi-layer SpatRast of total shortwave radiation
#' (W/m^2) matching the resolution of dtmf. If if `terrainshade = TRUE` a list of wrapped multi-layer
#' SpatRasts of (1) total direct and difuse shortwave radiation (W/m^2) and (2) diffuse-only radiation  (W/m^2).
#' @details radiation is downscaled by computing average fraction of clear-sky radiation
#' and then adjusting this variable to account for elevation using an emprical
#' adjustment calibrated against 0.05 degree data for the UK. If
#' `patchsim` is set to TRUE cloud cover patchiness is simulated the gstats package.
#' The parameter `nsim` determines the number of independent simulations and hence
#' the time intervals at which these simulations are performed. Simulated anomalies
#' due to local cloud patchiness are then interpolated temporally between these
#' periods. This ensures that, over shorter increments of say and hour, the location
#' clouds within each hour retain a degree of inter-dependence more realistically
#' simulating the trajectory of cloud as they move across the landscape.
#' @rdname swraddownscale
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @import terra
#' @import fields
#' @import gstat
#' @export
#'
#' @keywords spatial
#' @examples
#' climdata<- read_climdata(mesoclim::ukcpinput)
#' sw<-swdownscale(climdata$swrad,dtmf=terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim')), dtmc=climdata$dtm,terrainshade = TRUE)
#' sw<-sw$swf
#' terra::panel(c(sw[[1]],sw[[5]],sw[[10]],sw[[20]]),main=paste0("SW down ",c(1,5,10,20),"/05/2018"))
swdownscale<-function(swrad, tme=NA, dtmf, dtmc, patchsim = FALSE, nsim= dim(swrad)[3],
                      hor=NA, terrainshade = FALSE) {
  if(inherits(swrad,'SpatRaster') & inherits(tme,"logical")) tme<-as.POSIXlt(terra::time(swrad))
  # Work out whether daily or not
  ti<-round(as.numeric(tme[2])-as.numeric(tme[1]))
  # Check no NA in dtmc - convert to 0 elevation
  dtmc<-ifel(is.na(dtmc),0,dtmc)
  # Calculate clear sky fraction
  jd<-juldayvCpp(tme$year+1900, tme$mon+1, tme$mday)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  ll<-.latslonsfromr(dtmc)
  lats<-.rast(ll$lats,dtmc)
  lons<-.rast(ll$lons,dtmc)
  lats<-.is(mask(lats,dtmc))
  lons<-.is(mask(lons,dtmc))
  if (ti < 86400) hourly<-TRUE else hourly<-FALSE
  csr<-array(clearskyradmCpp(jd,lt,as.vector(lats),as.vector(lons),hourly),dim=dim(swrad))

  # Calculate clear-sky fraction and elevation adjust it
  csfc<-.is(swrad)/csr
  csfc[is.na(csfc)]<-0.5
  csfc[csfc>1]<-1
  csfc[csfc<0]<-0
  csfc<-.rast(csfc,dtmc)
  csf<-.cfcelev(csfc,dtmf,dtmc)
  if (patchsim) {
    if (crs(dtmc) != crs(dtmf)) {
      dtmp<-project(dtmc,crs(dtmf))
    } else dtmp<-dtmc
    af<-res(dtmp)[1]/res(dtmf)[1]
    csf<-.simpatch(csf,af,n=nsim,varn="csfrac")
  }
  # Calculate fine-res csf
  ll<-.latslonsfromr(dtmf)
  lats<-.rast(ll$lats,dtmf)
  lons<-.rast(ll$lons,dtmf)
  lats<-.is(mask(lats,dtmf))
  lons<-.is(mask(lons,dtmf))
  if (ti < 86400)  {  # Hourly
    lats<-.rast(ll$lats,dtmf)
    lons<-.rast(ll$lons,dtmf)
    lats<-.is(mask(lats,dtmf))
    lons<-.is(mask(lons,dtmf))
    csrf<-array(clearskyradmCpp(jd,lt,as.vector(lats),as.vector(lons)),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
  }  else {  # daily - simplified assumes constant lat across each row (not always case for OSGB)
    lats<-as.vector(ll$lats[,1])
    csrf<-clearskyradmCpp(jd,rep(12,length(jd)),lats,rep(0,length(lats)),hourly=FALSE)
    # csrf<-array(rep(csrf,each=dim(swrad)[2]),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
    csrf<-array(apply(csrf,2,rep,times=dim(dtmf)[2]),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
  }

  swradf<-.rast(csrf,dtmf)*csf
  if (terrainshade) {
    if (ti < 86400)  {  # Hourly
      swf<-matrix(.is(swradf),ncol=dim(swradf)[3])
      dp<-difpropmCpp(swf,jd,lt,as.vector(lats),as.vector(lons))
      dp<-array(dp,dim=dim(swradf))
      ze<-array(solzenmCpp(jd,lt,as.vector(lats),as.vector(lons)),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
      azi<-.solazi(lt, mean(lats,na.rm=TRUE),mean(lons,na.rm=TRUE),jd)
    } else {  # daily (expand to hour)
      # replicate clear sky fraction 24 times
      csrfh<-.ehr(csrf)
      # Calculate clear sky radiation hourly
      lth<-rep(c(0:23),dim(csrf)[3])
      jdh<-rep(jd,each=24)
      ll<-.latslonsfromr(dtmf)
      lats<-.rast(ll$lats,dtmf)
      lons<-.rast(ll$lons,dtmf)
      lats<-.is(mask(lats,dtmf))
      lons<-.is(mask(lons,dtmf))
      cs<-array(clearskyradmCpp(jdh,lth,as.vector(lats),as.vector(lons)),dim=c(dim(dtmf)[1:2],length(jdh)))
      rad<-cs*csrfh
      swf<-matrix(rad,ncol=length(jdh))
      dp<-difpropmCpp(swf,jdh,lth,as.vector(lats),as.vector(lons))
      dp<-array(dp,dim=dim(csrfh))
      ze<-array(solzenmCpp(jdh,lth,as.vector(lats),as.vector(lons)),dim=c(dim(dtmf)[1:2],dim(csrfh)[3]))
      azi<-.solazi(lth, mean(lats,na.rm=TRUE),mean(lons,na.rm=TRUE),jdh)
    }
    # Calculate horizon angle in 24 directions if required
    if(inherits(hor,'logical')){
      hor<-array(NA,dim=c(dim(dtmf)[1:2],24))
      for (i in 1:24) hor[,,i]<-.horizon(dtmf,(i-1)*15)
    } else hor<-.is(hor)
    # Match horizon angles to azimuths
    nhor<-dim(hor)[3]
    i<-round(azi/(360/nhor))+1; i[i==(nhor+1)]<-1
    hora<-hor[,,i]
    # Calculate terrain shading
    shadowmask<-hora*0+1
    alt<-(90-ze)*pi/180
    shadowmask[hora>tan(alt)]<-0
    shadowmask[(90-ze)<0]<-0
    # Calculate sky view
    #svf<-.rta(.skyview(dtmf),dim(swradf)[3])
    svf<-.rta(.skyview(dtmf),dim(shadowmask)[3])

    # Convert daily SW to hourly
    swradfh<-.ehr(.is(swradf))

    # Adjust radiation to account for sky view factor
    #drf<-.rast(dp*svf*.is(swradf),dtmf) # FAILS HERE FOR DAILY swradf still daily values
    drf<-dp*svf*swradfh
    #swf<-(1-dp)*shadowmask*.is(swradf)+dp*svf*.is(swradf)
    swf<-(1-dp)*shadowmask*swradfh+dp*svf*swradfh
    swf<-.rast(swf,dtmf)
    if (ti == 86400) {  # daily average across days
      swf<-hourtodayCpp(.is(swf),"mean")
      drf<-hourtodayCpp(.is(drf),"mean")
      swf<-.rast(swf,dtmf)
      drf<-.rast(drf,dtmf)
    }
    # Outputs
    terra::time(swf)<-tme
    terra::time(drf)<-tme
    out<-list(swf=swf,drf=drf)
  } else { # If no terrain shading
    out<-list(swf=swradf,drf=NA)
  }
  return(out)
}
#' @title Downscale wind speed accounting for elevation and terrain sheltering effects
#' @description The function `winddownscale` is used to spatially downscale windspeed,
#' performing adjustments for elevation and terrain sheltering
#' a specified wind direction.
#' @param wspeed a coarse-resolution array of wind speeds (m/s)
#' @param wdir a coarse-resolution array of wind directions (decimal degrees from north, from)
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (see details)
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `wspeed`.
#' @param wca - wind shelter coeeficient array as produced by `calculate_windcoeffs` or NA to calculate
#' @param zi height above ground (m) of `wspeed`.
#' @param zo height above ground (m) of output windspeeds (corrected to 0.2 if <0.2?)
#' @return a stacked SpatRast of wind speeds matching the resolution,
#' coordinate reference system and extent of `dtmf`.
#' @details Elevation effects are derived by sampling the dtms at intervals in
#' an upwind direction, determining the elevation difference form each focal cell and
#' performing a standard wind-height adjustment. Terrain sheltering is computed
#' from horizon angles following the method detailed in Maclean et al (2019) Methods
#' Ecol Evol 10: 280-290. By supplying three dtms, the algorithm is able to account for
#' elevation differences outside the boundaries of `dtmf`. The area covered by `dtmm` was
#' extend at least one `dtmc` grid cell beyond `dtmf`. Elevations must be in metres.
#' The coordinate reference system of `dtmf` must be such that x and y are also in metres.
#' `dtmm` and `dtmc` are reprojected to match the coordinate reference system of `dtmf`.
#' To enhance computational efficiency, Wind direction is assumed constant over the extent of
#' `dtmf`.
#' @import terra
#' @export
#' @rdname winddownscale
#' @seealso [windelev()]
#' @keywords spatial
#' @examples
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' wsf <- winddownscale(climdata$windspeed, climdata$winddir, dtmf, dtmm, climdata$dtm, zi=climdata$windheight_m)
#' terra::panel(c(wsf[[1]],wsf[[5]],wsf[[10]],wsf[[20]]),main=paste0("Windspeed down ",c(1,5,10,20),"/05/2018"))
winddownscale <- function(wspeed, wdir, dtmf, dtmm, dtmc, wca=NA, zi=10, zo = 2) {
  # Coarsen dtmm resolution if currently same as dtmf
  if(all(res(dtmm)==res(dtmf))){
    dtmm_res<-round(exp( ( log(terra::res(dtmc)[1]) + log(terra::res(dtmf)[1]) ) / 2 ))
    dtmm<-terra::aggregate(dtmm,dtmm_res / res(dtmf),  na.rm=TRUE)
  }
  # If not supplied, calculate terrain adjustment coefs for 8 directions at wind height zo
  if(class(wca)[1]=='logical') wca2<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo) else wca2<-wca
  # Calculate wind direction of centre of study area
  if(inherits(wdir,"SpatRaster")) wdr<-wdir else wdr<-.rast(wdir,dtmc)
  ll<-.latlongfromrast(wdr)
  xy<-data.frame(x=ll$long,y=ll$lat)
  wdir<-as.numeric(xx<-extract(wdr,xy))[-1]
  if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,mean,na.rm=TRUE)
  # Calculate array of shelter coefficients for every hour
  i<-round(wdir/45)%%8
  windmu<-wca2[,,i+1]
  # Resample wind speeds
  if(inherits(wspeed,"SpatRaster")) wsr<-wspeed else wsr<-.rast(wspeed,dtmc)
  # Adjust windspeed to output height zo
  if (zi != zo) wsr<-.windhgt(wsr, zi, zo)
  #if (uz !=2) wsr<-wsr*4.8699/log(67.8*uz-5.42)
  if (crs(wsr) != crs(dtmf)) {
    ws<-.is(project(wsr,dtmf))
  } else ws<-.is(.resample(wsr,dtmf))
  # Calculate wind speed
  ws<-ws*windmu
  # adjust back to output height if not 2 m NOT NEEDED as already at zo!!
  #if (uz !=2) ws<-ws*log(67.8*uz-5.42)/4.8699
  ws<-.rast(ws,dtmf)
  return(ws)
}
#' @title Downscale relative humidity accounting for downscaled temperature
#' @description The function `relhumdownscale` is used to spatially downscale relative
#' @param rh a coarse-resolution array of relative humidities (percentage)
#' @param tcc a coarse-resolution array of temperatures (deg C)
#' @param tcf a fine-resolution SpatRast of temperatures (deg C) as returned by
#' [temphrly_downscale()] or [tempdaily_downscale()]
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `rh` and `tc`.
#' @param rhmin optional single numeric value indicating the minimum realistically
#' attainable relative humidity within the study region (percentage).
#' @return a stacked SpatRast of relative humidities (percentage) matching the
#' resolution, coordinate reference system and extent of `tcf`.
#' @import terra
#' @export
#' @rdname relhumdownscale
#' @seealso [temphrly_downscale()] and [tempdaily_downscale()]
#' @keywords spatial
#' @examples
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' tmean_c<-mesoclim:::.hourtoday(temp_dailytohourly(climdata$tmin, climdata$tmax, climdata$tme),mean)
#' tcf<-tempdaily_downscale(climdata,sst=NA,dtmf=dtmf,dtmm=NA,cad=FALSE,coastal=FALSE)$tmean
#' rhf<-relhumdownscale(climdata$relhum,tcc=tmean_c,tcf=tcf,dtmc=climdata$dtm,rhmin=20)
#' terra::panel(c(rhf[[1]],rhf[[5]],rhf[[10]],rhf[[20]]),main=paste0("Rel humidity ",c(1,5,10,20),"/05/2018"))
relhumdownscale<-function(rh, tcc, tcf, dtmc, rhmin = 20) {
  eac<-.satvap(tcc)*rh/100
  if(!inherits(eac,"SpatRaster")) eac<-.rast(eac,dtmc)
  if (crs(eac) != crs(tcf)) eac<-project(eac,crs(tcf))
  eaf<-.resample(eac,tcf[[1]])
  rhf<-(eaf/.satvap(tcf))*100
  rhf<-ifel(rhf>100,100,rhf)
  rhf<-ifel(rhf<rhmin,rhmin,rhf)
  return(rhf)
}
#' @title Downscale precipitation accounting for elevation effects
#' @description The function is used to spatially downscale precipitation,
#' performing adjustments for elevation using one of two methods
#' @param prec a coarse-resolution array of precipitation (mm)
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `prec`.
#' @param method One of `Tps` or `Elev` (see details).
#' @param fast optional logical indicating whether to do a faster but less accurate
#' Thin-plate spline downscaling. Ignored if `method = "Elev"`
#' @param noraincut  single numeric value indicating rainfall amounts that should
#' be considered as no rain (see details).
#' @param patchsim optional logical indicating whether to simulate rain patchiness during
#' downscaling. More realistically captures intensity, but slower.
#' @param nsim optionally the number of independent rain patchiness simulations to perform. Outputs
#' are temporally interpolated (see details).
#' @return a stacked SpatRast of precipitation values matching the resolution,
#' coordinate reference system and extent of `dtmf`.
#' @details Precipitation is downscaled by computing the total rainfall amount and
#' fraction of rain days and then adjusting these variables to account for elevation
#' using either a Thin-plate-spine model (`method = "Tps"`) or by performing an
#' empirical adjustment calibrated to UK precipitation data (`method = "Elev"`).
#' If method is set to `Tps` users have the option to specify whether to perform
#' a fast, but slightly less accurate Thin-plate spline using [fields::fastTps()] instead
#' of [fields::Tps()]. One totals have been derived, precipitation for each day is adjusted
#' to ensure that the total amount and number of precipitation days match that expected
#' by a small amount of precipitation on no precipitation days that are the wettest regionally, or
#' by setting precipitation days with little rain to zero as required. To accommodate that some
#' precipitation datasets erroneously include small amounts precipitation (<0.01 mm/day)
#' on no precipitation days, users have the option to set a cut-off via `noraincut`. If
#' `patchsim` is set to TRUE rainfall patchiness caused by e.g. localized rain showed
#' is simulated the gstats package. The parameter `nsim` determines the number of
#' independent simulations and hence the time intervals at which these simulations
#' are performed. Simulated anomalies due to local precipitation events are
#' then interpolated temporally between these periods. This ensures that, over shorter
#' increments of say and hour, the location and magnitude of these local precipitation
#' events within each hour retain a degree of inter-dependence more realistically
#' simulating the trajectory of low pressure systems.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @import terra
#' @import fields
#' @import gstat
#' @export
#' @keywords spatial
#' @examples
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' prcf<-precipdownscale(climdata$prec, dtmf, climdata$dtm, method="Elev", noraincut=0.01)
#' terra::panel(c(prcf[[1]],prcf[[5]],prcf[[12]],prcf[[26]]),main=paste0("Precipitation (mm) ",c(1,5,12,26),"/05/2018"))
precipdownscale <- function(prec, dtmf, dtmc, method = "Tps", fast = TRUE, noraincut = 0, patchsim = FALSE, nsim = dim(prec)[3]){
  if(!inherits(prec,"SpatRaster")) prec<-.rast(prec,dtmc)

  # check how many non NA cells
  if (method != "Tps" & method != "Elev") stop("method must be one of Tps or Elev")
  v<-as.vector(prec[[1]])
  v<-v[is.na(v) == FALSE]
  if (method == "Tps" & length(v) < 500) {
    warning("Not enough non NA cells for sensible thin-plate spline downscale. Changed method to Elev")
    method <- "Elev"
  }
  # Fill gaps in prec
  p2<-resample(aggregate(prec,2,na.rm=TRUE),prec)
  p5<-resample(aggregate(prec,5,na.rm=TRUE),prec)
  a1<-.is(prec)
  a2<-.is(p2)
  a5<-.is(p5)
  s<-which(is.na(a1))
  a1[s]<-a2[s]
  s<-which(is.na(a1))
  a1[s]<-a5[s]
  prec<-.rast(a1,prec)
  if (crs(prec) != crs(dtmf)) prec<-project(prec,crs(dtmf))
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  # Crop prec to area not too far around dtmf
  e<-ext(dtmf)
  ro<-res(prec)[1]
  e$xmin<-e$xmin-10*ro
  e$xmax<-e$xmax+10*ro
  e$ymin<-e$ymin-10*ro
  e$ymax<-e$ymax+10*ro
  prec<-crop(prec,e)
  dtmc<-crop(dtmc,e)
  # Fill in final gaps
  a1<-.is(prec)
  v<-apply(a1,3,median,na.rm=T)
  aa<-.vta(v,prec)
  s<-which(is.na(a1))
  a1[s]<-aa[s]
  prec<-.rast(a1,prec)
  # Calculate total rainfall
  m<-.is(prec)
  m1<-apply(m,c(1,2),sum)
  r1<-.rast(m1,prec)
  # Calculate rain day fraction
  m2<-m
  m2[m>=noraincut]<-1
  m2[m<noraincut]<-0
  m2<-apply(m2,c(1,2),sum)/dim(m2)[3]
  r2<-.rast(m2,prec)
  # Calculate wettest days regionally to enable sensible assignment of rain
  # to no rain days as needed
  rrain<-apply(.is(prec),3,sum,na.rm=TRUE)
  rr2<-as.numeric(.mav(rrain,10))
  s<-which(rrain==0)
  rrain[s]<-rr2[s]
  s<-which(rrain==0)
  rrain[s]<-0.1
  # Calculate resampled rain
  rf3<-.resample(prec,dtmf)
  rf3<-mask(rf3,dtmf)
  if (method == "Tps") {
    # Downscaled rain total
    rf1<-Tpsdownscale(r1, dtmc, dtmf, method = "log", fast)  # Total rain
    # Downscaled rain day fraction
    rf2<-Tpsdownscale(r2, dtmc, dtmf, method = "logit", fast) # Rain day fraction
  } else {
    # Convert NA to zero in dtmc
    dtmc<-ifel(is.na(dtmc),0,dtmc)
    # Downscaled rain total
    dtmcf<-resample(dtmc,dtmf)
    edif<-dtmf-dtmcf
    prat<-9.039606e-03+1.818067e-03*edif-2.923351e-04*dtmf-6.471352e-07*dtmf*edif
    rtmc<-.resample(r1,dtmf)
    rf1<-rtmc*exp(prat)
    # Downscaled rain day fraction
    prat<- -2.787208e-03+-2.787208e-03*edif-2.787208e-03*dtmf-4.503742e-07*dtmf*edif
    rdmc<-.resample(r2,dtmf)
    rf2<-rdmc*exp(prat)
    rf2[rf2<0]<-0
    rf2[rf2>1]<-1
  }
  # Simulate patchiness
  if (patchsim) {
    af<-(res(dtmc)/res(rf3))[1]
    rf3<-.simpatch(rf3,af,nsim,varn="precip")
  }
  # Correct resampled rain for predicted rain total and rain days using Cpp function
  a<-as.array(rf3)
  mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
  rtot<-as.vector(t(rf1))
  rfrac<-as.vector(t(rf2))
  mm<-rainadjustm(mm,rrain,rfrac,rtot)
  a2<-array(mm,dim=dim(a))
  # Convert to raster
  precf<-.rast(a2,dtmf)
  return(precf)
}

#' @title Spatially downscale all climate variables
#' @description Spatially downscales coarse-resolution climate data
#' @param climdata a `climdata` model object containing climate data of the same format as `ukcp18toclimarray`
#' @param sst a SpatRast of sea-surface temperature data (deg C) that overlaps with climdata$tme
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details under [temphrly_downscale()]).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param wca optionally, an array of wind shelter coeffs in 8 directions
#' matching the x,y dimensions of `dtmf` as calculated by `calculate_windcoeffs`. Calculated if
#' not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param thgto height above ground of temperature output.
#' @param whgto height above ground of wind speed output.
#' @param include_tmean if TRUE and daily data will output mean daily temperature based on hourly downsccaling
#' @param rhmin minimum relative humidity (set to avoid relative humidity dropping too low
#' during down-scaling). Default 20 (percent).
#' @param pksealevel optional logical indicating whether input pressure data represent
#' sea-level pressure (default TRUE).
#' @param pathsim optional logical indicating whether to simulate cloud cover and
#' rainfall patchiness see details under [swdownscale()] and [precipdownscale()].
#' @param difr optional logical indicating whether to return diffuse radiation (
#' default TRUE).
#' @param terrainshade optional logical indicating whether to account for terrain shading
#' during shortwave radiation downscale.
#' @param precipmethod One of `Tps` or `Elev` indicating whether to account for
#' elevation effects using a Thin-plate spline model or en emprical adjustment
#' (see details under [precipdownscale()])
#' @param fast optional logical indicating whether to use fast Thin-plate spline
#' down-scaling (see details under [precipdownscale()]).
#' @param noraincut numeric value below which low precipitation amounts are set to
#' zero (see details under [precipdownscale()]).
#' @param toArrays if TRUE climate variables are output as 3D arrays if FALSE as Spatarasters
#' @return a list of multi-layer SpatRast or 3D arrays of downscaled climate variables and metadata as follows:
#' \describe{
#'    \item{dtm}{Digital elevation of downscaled area in metres (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C) - is replaced with `tmax` and `tmin` if input data is daily and `tmean` if `include_tmean`=TRUE.}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2) - not returned unless `terrainshade = TRUE`}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed}{at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @keywords spatial
#' @examples
#'  dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#'  dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#'  climdata<-read_climdata(mesoclim::ukcpinput)
#'  sst<-terra::unwrap(mesoclim::ukcp18sst)
#'  mesodat<-spatialdownscale(climdata, sst, dtmf, dtmm, include_tmean=TRUE,noraincut=0.01)
#'  for(n in 5:length(mesodat)) terra::plot(mesodat[[n]][[12]],main=paste(names(mesodat)[n],mesodat$tme[n]))
spatialdownscale<-function(climdata, sst, dtmf, dtmm = NA, basins = NA, wca=NA, skyview=NA, horizon=NA,
                              cad = TRUE,coastal = TRUE, thgto =2, whgto=2,include_tmean=FALSE,
                               rhmin = 20, pksealevel = TRUE, patchsim = FALSE, terrainshade = TRUE,
                               precipmethod = "Elev",fast = TRUE, noraincut = 0, toArrays=FALSE) {

  # Convert variables - unpack any wrapped spatRasters and convert arrays to spatraster
  input_class<-lapply(lapply(climdata,class),"[",1)
  if(any(input_class=="PackedSpatRaster")) climdata[which(input_class=="PackedSpatRaster")]<-lapply(climdata[which(input_class=="PackedSpatRaster")],unwrap)
  if(any(input_class=="array")) climdata[which(input_class=="array")]<-lapply(climdata[which(input_class=="array")],.rast,tem=climdata$dtm)
  if(inherits(sst,"PackedSpatRaster")) sst<-unwrap(sst)

  # Find out whether daily or hourly
  tme<-as.POSIXlt(climdata$tme,tz="UTC")
  tint<-as.numeric(tme[2])-as.numeric(tme[1])
  hourly<-TRUE
  if (abs(tint-86400) < 5) hourly=FALSE

  # Check likely memory use compared with free RAM - even daily requires converting to hourly for some downscaling
   mmry<-mem_info(dtmf, n=length(tme)*24, print=FALSE)
   if (mmry["needed"]>(0.5*mmry["available"]) & mmry["needed"]<mmry["available"]) warning("High free memory use predicted - consider using spatialdownscale_tile option!!!")
   if (mmry["needed"]>mmry["available"]) warning("Memory demand predicted to exceed available memory - use spatialdownscale_tile option!!!")

  # Extract temp variables - calculate daily means if required
  if (hourly) {
    tc<-climdata$temp
  } else {
    tmin<-climdata$tmin
    tmax<-climdata$tmax
    tmean<-.hourtoday(temp_dailytohourly(tmin, tmax, climdata$tme),mean) # required for relhum downscaling
  }

  # Extract other variables
  dtmc<-climdata$dtm
  rh<-climdata$relhum
  pk<-climdata$pres
  wspeed<-climdata$windspeed
  wdir<-climdata$winddir
  swrad<-climdata$swrad
  lwrad<-climdata$lwrad
  prec<-climdata$prec
  whgti<-climdata$windheight_m
  thgti<-climdata$tempheight_m

  # ================== Downscale variables ===============================  #
  message('Downscaling wind...')
  uzf<-winddownscale(wspeed,wdir,dtmf,dtmm,dtmc,wca,whgti,whgto)

  # winddir
  uu<-wspeed*cos(wdir*pi/180)
  vv<-wspeed*sin(wdir*pi/180)
  if (crs(uu) != crs(dtmf)) uu<-project(uu,crs(dtmf))
  if (crs(vv) != crs(dtmf)) vv<-project(vv,crs(dtmf))
  uu<-.resample(uu,dtmf,method="cubicspline")
  vv<-.resample(vv,dtmf,method="cubicspline")
  wdf<-atan2(vv,uu)*180/pi
  wdf<-.rast(.is(wdf)%%360,dtmf)
  wdf<-mask(wdf,dtmf)

  message('Downscaling temperature...')
  # for daily data will calculate daily means from tmin/tmax temporally downscaled to hourly temp
  if (hourly) {
    tcf<-temphrly_downscale(climdata,sst,dtmf,dtmm,basins,uzf,cad,coastal,thgto,whgto)
  } else {
    dailytemps<-tempdaily_downscale(climdata,tmean,sst,dtmf,dtmm,basins,uzf,cad,coastal,thgto,whgto)
    tminf<-dailytemps$tmin
    tmaxf<-dailytemps$tmax
    tmeanf<-dailytemps$tmean
    remove(dailytemps)
  }

  message('Downscaling relative humidity')
  if (hourly) {
    rhf<-relhumdownscale(rh,tc,tcf,dtmc,rhmin)
  } else {
    rhf<-relhumdownscale(rh,tmean,tmeanf,dtmc,rhmin)
  }

  message('Downscaling pressure...')
  pkf<-presdownscale(pk,dtmf,dtmc,sealevel=pksealevel)

  message('Downscaling SW radiation...')
  if (patchsim) {
    # Calculate mean wind speed
    mws<-mean(as.vector(uzf),na.rm=TRUE)
    # Calculate distance traversed in one timestep
    tstep<-as.numeric(tme[2])-as.numeric(tme[1])
    dtr<-tstep*mws
    # Calculate fraction of study area traversed in one time step
    sze<-max(dim(dtmf)[1]*res(dtmf)[1],dim(dtmf)[2]*res(dtmf)[2])
    nsim<-round((dtr/sze)*dim(swrad)[3])+1
    if (nsim > dim(swrad)[3]) nsim<-dim(swrad)[3]
  } else nsim<-dim(swrad)[3]
  swf<-swdownscale(swrad,tme,dtmf,dtmc,patchsim,nsim,horizon,terrainshade)
  totswrad<-swf$swf
  if (terrainshade) {
    difrad<-swf$drf
    #difrad<-difrad*swf
  }  else difrad = NA

  message('Downscaling LW radiation...')
  if(hourly) lwf<-lwdownscale(lwrad, tc, tcf, tme, dtmf, dtmc, skyview, terrainshade) else lwf<-lwdownscale(lwrad, tmean, tmeanf, tme, dtmf, dtmc, skyview, terrainshade)

  message('Downscaling precipitation...')
  precf<-precipdownscale(prec,dtmf,dtmc,precipmethod,fast,noraincut,patchsim,nsim)

  message('Formatting output...')
  if (hourly) {
    terra::time(tcf)<-climdata$tme
    names(tcf)<-climdata$tme
  } else {
    terra::time(tminf)<-climdata$tme
    terra::time(tmaxf)<-climdata$tme
    terra::time(tmeanf)<-climdata$tme
    names(tminf)<-climdata$tme
    names(tmaxf)<-climdata$tme
    names(tmeanf)<-climdata$tme
  }
  terra::time(rhf)<-climdata$tme
  terra::time(pkf)<-climdata$tme
  terra::time(totswrad)<-climdata$tme
  terra::time(lwf)<-climdata$tme
  terra::time(uzf)<-climdata$tme
  terra::time(wdf)<-climdata$tme
  terra::time(precf)<-climdata$tme
  names(rhf)<-climdata$tme
  names(pkf)<-climdata$tme
  names(totswrad)<-climdata$tme
  names(lwf)<-climdata$tme
  names(uzf)<-climdata$tme
  names(wdf)<-climdata$tme
  names(precf)<-climdata$tme

  out<-list(dtm=dtmf,tme=climdata$tme,windheight_m=whgto,tempheight_m=thgto)
  if (hourly) {
    climout<-list(temp=tcf,relhum=rhf,pres=pkf,swrad=totswrad,lwrad=lwf,windspeed=uzf,winddir=wdf,prec=precf)
    if(toArrays) climout<-lapply(climout,as.array)
    out<-c(out,climout)
  } else {
    if(!include_tmean) climout<-list( tmin=tminf,tmax=tmaxf,relhum=rhf,pres=pkf,swrad=totswrad,lwrad=lwf,windspeed=uzf,winddir=wdf,prec=precf)
    if(include_tmean) climout<-list( tmean=tmeanf,tmin=tminf,tmax=tmaxf,relhum=rhf,pres=pkf,swrad=totswrad,lwrad=lwf,windspeed=uzf,winddir=wdf,prec=precf)
    if(toArrays) climout<-lapply(climout,as.array)
    out<-c(out,climout)
  }
  if (terrainshade){
    terra::time(difrad)<-climdata$tme
    names(difrad)<-climdata$tme
    out$difrad<-difrad
  }
  return(out)
}

