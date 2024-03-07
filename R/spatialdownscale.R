#' @title Downscale temperature
#' @description Downscales coarse-resolution temperature data accounting for elevation
#' effects and optionally cold air drainage and coastal effects
#' @param climdata a `climdata` model object containing climate data of the same format as `era5climdata`
#' @param SST a coarse resolution SpatRast of hourly sea-surface temperature data (deg C) with no NAs, as
#' returned by [SSTinterpolate()]
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param u2 optionally, a SpatRast of high resolution wind speeds as returned by [windownscale()].
#' Calculated if not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param refhgt height above ground of temperature measurements in `climdata`.
#' @param uhgt height above gorund of wind speed measurements in `climdata`
#' to resolution of `dtmf`.
#' @return a multi-layer SpatRast of downscaled temperatures (deg C) matching the
#' resolution of dtmf.
#' @details Cold air drainage is calculated by delineating hydrological basins and
#' calculating flow accumulation and the elevation difference from the highest point of the basin.
#' Cold-air drainage is assumed to occur when atmospheric stability is high, namely when
#' the radiation balance is negative and wind speeds are low. Coastal effects are
#' calculated by determining the ratio of land to sea pixels in an upwind direction.
#' The provision of `dtmm` allows this ratio to be derived accoutning for land and
#' sea outside the boundaries of the study area.
#' @rdname tempdownscale
#' @import terra
#' @export
#' @examples
#' # Takes ~90 seconds to run
#' tmf <- tempdownscale(era5data, rast(era5sst), rast(dtmf), rast(dtmm))
#' plot(tmf[[5]]) # Coldest temperature
tempdownscale<-function(climdata, SST, dtmf, dtmm = NA, basins = NA, u2 = NA,
                        cad = TRUE, coastal = TRUE, refhgt = 2, uhgt = 2) {
  if (class(dtmm) == "logical" & coastal) stop("dtmm needed for calculating coastal effects")
  # Calculate elevation effects
  dtmc<-rast(climdata$dtmc)
  tc<-climdata$climarray$temp
  rh<-climdata$climarray$relhum
  pk<-climdata$climarray$pres
  tcf<-.rast(tc,dtmc)
  if (crs(tcf) != crs(dtmf)) tcf<-project(tcf,crs(dtmf))
  # Calculate resampled coarse-res temp
  tcc<-resample(tcf,dtmf)
  # Calculate elevation effects
  tcf<-.tempelev(tc,dtmf,dtmc,rh,pk)
  if (cad) {
    tcad<-.tempcad(climdata,dtmf,basins,refhgt)
    tcf<-tcf+tcad
  }
  if (coastal) {
    te<-mean(.is(SST))
    if (is.na(te)) {
      SST<-SSTinterpolate(SST,climdata$tme,climdata$tme)
    }
    if (crs(SST) != crs(dtmf)) SST<-project(SST,crs(dtmf))
    SST<-resample(SST,dtmf)
    wspeed<-climdata$climarray$windspeed
    wdir<-climdata$climarray$winddir
    u2<-winddownscale(wspeed,wdir,dtmf,dtmm,dtmc,uhgt)
    tcf<-.tempcoastal(tcf,SST,u2,wdir,dtmf,dtmm,dtmc)
  }
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
#' @examples
#' pk <- presdownscale(era5data$climarray$pres, rast(dtmf), rast(era5data$dtmc))
#' plot(pk[[1]])
presdownscale<-function(pk, dtmf, dtmc, sealevel = TRUE) {
  if (class(pk)[1] == "array") pk<-.rast(pk,dtmc)
  if (class(dtmc)[1] == "logical")  dtmc<-resample(dtmf,pk)
  n<-dim(pk)[3]
  if (sealevel == FALSE) {
    pk<-.is(pk)/(((293-0.0065*.rta(dtmc,n))/293)^5.26)
    pk<-.rast(pk,dtmc)
  }
  if (crs(dtmf) != crs(dtmc)) pk<-project(pk,crs(dtmf))
  pkf<-.is(resample(pk,dtmf))
  pkf<-.rast(pkf*((293-0.0065*.rta(dtmf,n))/293)^5.26,dtmf)
  return(pkf)
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
#' @param terrainshade optional logical indicating whether to adjust incoming radiation
#' for terrain shading. If TRUE, total downward shortwave radiation is partioned into its
#' direct and diffuse components, the latter adjusted by sky view and the former set
#' to zero if the sun is below the horizon, which is computed explicitly.
#' @return if `terrainshade = FALSE` a multi-layer SpatRast of shortwave radiation
#' (W/m^2) matching the resolution of dtmf. If if `terrainshade = TRUE` a list of wrapped multi-layer
#' SpatRasts of (1) shortwave radiation (W/m^2) and (2) the diffuse fraction (range 0-1).
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
#' @import terra, fields, gstat
#' @export
#'
#' @examples
#' # Extract model inputs form inbuilt datasets
#' swrad <- mesoclim::era5data$climarray$swrad
#' tme <- as.POSIXlt(era5data$tme, tz = "UTC")
#' dtmc <- rast(mesoclim::era5data$dtmc)
#' dtmf <- rast(mesoclim::dtmf)
#' # Run with default options
#' swradf<-swdownscale(swrad, tme, dtmf, dtmc)
#' # plot output
#' plot(swradf[[13]], main = as.character(tme)[13])
swdownscale<-function(swrad, tme, dtmf, dtmc, patchsim = FALSE, nsim= dim(swrad)[3],
                      difr = FALSE) {
  # Work out whether daily or not
  ti<-round(as.numeric(tme[2])-as.numeric(tme[1]))
  # Calculate clear sky fraction
  jd<-juldayvCpp(tme$year+1900, tme$mon+1, tme$mday)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  ll<-.latslonsfromr(dtmc)
  if (ti < 86400)  {  # Hourly
    lats<-.rast(ll$lats,dtmc)
    lons<-.rast(ll$lons,dtmc)
    lats<-.is(mask(lats,dtmc))
    lons<-.is(mask(lons,dtmc))
    csr<-array(clearskyradmCpp(jd,lt,as.vector(lats),as.vector(lons)),dim=dim(swrad))
  }  else {  # daily
    lats<-as.vector(ll$lats[,1])
    csr<-clearskyradmCpp(jd,rep(12,length(jd)),lats,rep(0,length(lats)),hourly=FALSE)
    csr<-array(rep(csr,each=dim(swrad)[2]),dim(swrad))
  }
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
  if (ti < 86400)  {  # Hourly
    lats<-.rast(ll$lats,dtmf)
    lons<-.rast(ll$lons,dtmf)
    lats<-.is(mask(lats,dtmf))
    lons<-.is(mask(lons,dtmf))
    csrf<-array(clearskyradmCpp(jd,lt,as.vector(lats),as.vector(lons)),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
  }  else {  # daily
    lats<-as.vector(ll$lats[,1])
    csrf<-clearskyradmCpp(jd,rep(12,length(jd)),lats,rep(0,length(lats)),hourly=FALSE)
    csrf<-array(rep(csrf,each=dim(swrad)[2]),dim=c(dim(dtmf)[1:2],dim(swrad)[3]))
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
    # Calculate hirozon angle in 24 directions
    hor<-array(NA,dim=c(dim(dtmf)[1:2],24))
    for (i in 1:24) hor[,,i]<-.horizon(dtmf,(i-1)*15)
    i<-round(azi/15)+1; i[i==25]<-1
    hora<-hor[,,i]
    # Calculate terrain shading
    shadowmask<-hora*0+1
    alt<-(90-ze)*pi/180
    shadowmask[hora>tan(alt)]<-0
    shadowmask[(90-ze)<0]<-0
    # Calculate sky view
    svf<-.rta(.skyview(dtmf),dim(swradf)[3])
    # Adjust radiation to account for sky view factor
    drf<-.rast(dp*svf*.is(swradf),dtmf)
    swf<-(1-dp)*shadowmask*.is(swradf)+dp*svf*.is(swradf)
    swf<-.rast(swf,dtmf)
    if (ti == 86400) {  # daily average across days
      swf<-hourtodayCpp(.is(swf),"mean")
      drf<-hourtodayCpp(.is(drf),"mean")
      swf<-.rast(swf,dtmf)
      drf<-.rast(swf,dtmf)
    }
    swf<-wrap(swf)
    drf<-wrap(drf)
    out<-list(swf=swf,drf=drf)
  } else {
    out<-swradf
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
#' @param uz height above ground (m) of wind speed measurement in `wspeed`.
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
#' @examples
#' wspeed <- era5data$climarray$windspeed
#' wdir <- era5data$climarray$winddir
#' # Takes a few seconds to run
#' wsfine <- winddownscale(wspeed, wdir, rast(dtmf), rast(dtmm), rast(era5data$dtmc))
#' plot(wsfine[[1]])
winddownscale <- function(wspeed, wdir, dtmf, dtmm, dtmc, uz = 2) {
  # Calculate wind speed in each of 8 directions
  wca<-array(NA,dim=c(dim(dtmf)[1:2],8))
  for (i in 0:7) wca[,,i+1]<-.is(windelev(dtmf,dtmm,dtmc,i*45,uz))
  # smooth results
  wca2<-wca
  for (i in 0:7) wca2[,,i+1]<-0.25*wca[,,(i-1)%%8+1]+0.5*wca[,,i%%8+1]+0.25*wca[,,(i+1)%%8+1]
  # Calculate wind direction of centre of study area
  wdr<-.rast(wdir,dtmc)
  ll<-.latlongfromrast(wdr)
  xy<-data.frame(x=ll$long,y=ll$lat)
  wdir<-as.numeric(xx<-extract(wdr,xy))[-1]
  if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,mean,na.rm=TRUE)
  # Calculate array of shelter coefficients for every hour
  i<-round(wdir/45)%%8
  windmu<-wca2[,,i+1]
  # Resample wind speeds
  wsr<-.rast(wspeed,dtmc)
  # adjust to 2 m if not 2m
  if (uz !=2) wsr<-wsr*4.8699/log(67.8*uz-5.42)
  if (crs(wsr) != crs(dtmf)) {
    ws<-.is(project(wsr,dtmf))
  } else ws<-.is(resample(wsr,dtmf))
  # Calculate wind speed
  ws<-ws*windmu
  # adjust back to uz if not 2 m
  if (uz !=2) ws<-ws*log(67.8*uz-5.42)/4.8699
  ws<-.rast(ws,dtmf)
  return(ws)
}
#' @title Downscale relative humidity accounting for downscaled temperature
#' @description The function `relhumdownscale` is used to spatially downscale relative
#' @param rh a coarse-resolution array of relative humidities (percentage)
#' @param tcc a coarse-resolution array of wind temperatures (deg C)
#' @param tcf a fine-resolution SpatRast of temperatures (deg C) as returned by
#' [tempdownscale()]
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `rh` and `tc`.
#' @param rhmin optional single numeric value indicating the minimum realistically
#' attainable relative humidity within the study region (percentage).
#' @return a stacked SpatRast of relative humidities (percentage) matching the
#' resolution, coordinate reference system and extent of `tcf`.
#' @import terra
#' @export
#' @rdname relhumdownscale
#' @seealso [tempdownscale()]
relhumdownscale<-function(rh, tcc, tcf, dtmc, rhmin = 0) {
  eac<-.satvap(tcc)*rh/100
  eac<-.rast(eac,dtmc)
  if (crs(eac) != crs(tcf)) eac<-project(eac,crs(tcf))
  eaf<-resample(eac,tcf[[1]])
  rhf<-(.is(eaf)/.satvap(.is(tcf)))*100
  rhf<-.rast(rhf,tcf[[1]])
  rhf[rhf>100]<-100
  rhf[rhf<rhmin]<-rhmin
  return(rhf)
}
#' @title Downscale precipitation accounting for elevation effects
#' @description The function `precipdownscale` is used to spatially downscale precipitation,
#' performing adjustments for elevation using one of two methods
#' @param prec a coarse-resolution array of precipitation (mm)
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `prec`.
#' @param method One of `Tps` or `Elev` (see details).
#' @param fast optional logical indicating whether to do a faster but less accurate
#' Thin-plate spline downscaling. Ignored if `method = "Elev"`
#' @param noraincut optional single numeric value indicating rainfall amounts that should
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
#' @import terra, fields, gstat
#' @export
precipdownscale <- function(prec, dtmf, dtmc, method = "Tps", fast = TRUE, noraincut = 0, patchsim = FALSE, nsim = dim(prec)[3]) {
  prec<-.rast(prec,dtmc)
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
  rf3<-resample(prec,dtmf)
  rf3<-mask(rf3,dtmf)
  if (method == "Tps") {
    # Downscaled rain total
    rf1<-Tpsdownscale(r1, dtmc, dtmf, method = "log", fast)  # Total rain
    # Downscaled rain day fraction
    rf2<-Tpsdownscale(r2, dtmc, dtmf, method = "logit", fast) # Rain day fraction
  } else {
    # Downscaled rain total
    dtmcf<-resample(dtmc,dtmf)
    edif<-dtmf-dtmcf
    prat<-9.039606e-03+1.818067e-03*edif-2.923351e-04*dtmf-6.471352e-07*dtmf*edif
    rtmc<-resample(r1,dtmf)
    rf1<-rtmc*exp(prat)
    # Downscaled rain day fraction
    prat<- -2.787208e-03+-2.787208e-03*edif-2.787208e-03*dtmf-4.503742e-07*dtmf*edif
    rdmc<-resample(r2,dtmf)
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
#' @param climdata a `climdata` model object containing climate data of the same format as `era5climdata`
#' @param SST a coarse resolution SpatRast of hourly sea-surface temperature data (deg C) with no NAs, as
#' returned by [SSTinterpolate()]
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details under [tempdownscale()]).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param refhgt height above ground of temperature measurements in `climdata`.
#' @param uhgt height above gorund of wind speed measurements in `climdata`
#' to resolution of `dtmf`.
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
#' @return a list of multi-layer SpatRast of downscaled climate variables as follows:
#' \describe{
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
#'  If input data are hourly `temp` is replaced with `tmax` and `tmin` representing
#'  daily maximum and minimum values. The variable `difrad` is not returned unless
#'  `terrainshade = TRUE`.
spatialdownscale<-function(climdata, SST, dtmf, dtmm = NA, basins = NA, cad = TRUE,
                           coastal = TRUE, refhgt = 2, uhgt = 2, rhmin = 20,
                           pksealevel = TRUE, patchsim = TRUE, terrainshade = TRUE, precipmethod = "Elev",
                           fast = TRUE, noraincut = 0) {
  # Sort out SST
  tme<-as.POSIXlt(climdata$tme,tz="UTC")
  # Find out whether daily or hourly
  tint<-as.numeric(tme[2])-as.numeric(tme[1])
  hourly<-TRUE
  if (abs(tint-86400) < 5) hourly=FALSE
  SST<-SSTinterpolate(SST,tme,tme)
  # Extract variables
  dtmc<-rast(climdata$dtmc)
  climarray<-climdata$climarray
  if (hourly) {
    tc<-climarray$temp
  } else {
    tmin<-climarray$tmin
    tmax<-climarray$tmax
  }
  rh<-climarray$relhum
  pk<-climarray$pres
  wspeed<-climarray$windspeed
  wdir<-climarray$winddir
  swrad<-climarray$swrad
  lwrad<-climarray$lwrad
  prec<-climarray$prec
  # ================== Downscale variables ===============================  #
  # windspeed"
  uzf<-winddownscale(wspeed,wdir,dtmf,dtmm,dtmc,uhgt)
  if (uhgt!=2) {
    u2<-.is(uzf)*4.8699/log(67.8*uhgt-5.42)
    u2<-.rast(u2,dtmf)
  } else u2<-uzf
  # winddir"
  uu<-.rast(wspeed*cos(wdir*pi/180),dtmc)
  vv<-.rast(wspeed*sin(wdir*pi/180),dtmc)
  if (crs(uu) != crs(dtmf)) uu<-project(uu,crs(dtmf))
  if (crs(vv) != crs(dtmf)) vv<-project(vv,crs(dtmf))
  uu<-resample(uu,dtmf)
  vv<-resample(vv,dtmf)
  wdf<-atan2(vv,uu)*180/pi
  wdf<-.rast(.is(wdf)%%360,dtmf)
  wdf<-mask(wdf,dtmf)
  if (hourly) {
    tcf<-tempdownscale(climdata,SST,dtmf,dtmm,basins,u2,cad,coastal,refhgt,uhgt)
  } else {
    climdata$climarray$temp<-tmin
    tminf<-tempdownscale(climdata,SST,dtmf,dtmm,basins,u2,cad,coastal,refhgt,uhgt)
    climdata$climarray$temp<-tmax
    tmaxf<-tempdownscale(climdata,SST,dtmf,dtmm,basins,u2,cad,coastal,refhgt,uhgt)
  }
  # Relative humidity
  if (hourly) {
    rhf<-relhumdownscale(rh,tc,tcf,dtmc,rhmin)
  } else {
    rhf<-relhumdownscale(rh,(tmin+tmax)/2,(tminf+tmaxf)/2,dtmc,rhmin)
  }
  # Pressure
  pkf<-presdownscale(pk,dtmf,dtmc,pksealevel)
  # swrad
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
  # Sw radiation
  swf<-swdownscale(swrad,tme,dtmf,dtmc,patchsim,nsim,difr)
  if (terrainshade) {
    difrad<-rast(swf$dp)
    swf<-rast(swf$swradf)
    difrad<-difrad*swf
  }  else difr = NA
  # lwrad
  lwf<-.rast(lwrad,dtmc)
  if (crs(lwf) != crs(dtmf)) lwf<-project(lwf,crs(dtmf))
  lwf<-resample(lwf,dtmf)
  if (terrainshade) {
    svf<-.rta(.skyview(dtmf),dim(lwf)[3])
    lwf<-.rast(.is(lwf)*svf,dtmf)
  }
  # prech
  precf<-precipdownscale(prec,dtmf,dtmc,precipmethod,fast,noraincut,patchsim,nsim)
  # return values
  if (hourly) {
    tcf<-wrap(tcf)
  } else {
    tminf<-wrap(tminf)
    tmaxf<-wrap(tmaxf)
  }
  rhf<-wrap(rhf)
  pkf<-wrap(pkf)
  swf<-wrap(swf)
  if (difr) difrad<-wrap(difrad)
  lwf<-wrap(lwf)
  uzf<-wrap(uzf)
  wdf<-wrap(wdf)
  if (hourly) {
    out<-return(list(temp=tcf,relhum=rhf,pres=pkf,swrad=swf,difrad=difrad,lwrad=lwf,windspeed=uzf,winddir=wdf,prec=precf))
  } else {
    out<-return(list(tmin=tminf,tmax=tmaxf,relhum=rhf,pres=pkf,swrad=swf,difrad=difrad,lwrad=lwf,windspeed=uzf,winddir=wdf,prec=precf))
  }
  return(out)
}
