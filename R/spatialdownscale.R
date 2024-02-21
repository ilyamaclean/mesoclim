# ============================================================================ #
# ~~~~~~~~~~~~ Temperature downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Downscale temperature with elevation effects
#' @description Downscales an array of temperature data applying elevation effects
#' @param tc a SpatRast or array of temperatures (deg C). If an array dtmc must
#' be provided.
#' @param dtmf a fine-resolution SpatRast of elevations. Temperatures down-scaled
#' to resolution of `dtmf`.
#' @param dtmc optional SpatRast of elevations matching resolution of `tc`. If
#' not supplied, `tc` must be a SpatRast and `dtmc` is derived by resampling `dtmf`
#' to resolution of `tc` and extents must match.
#' @param rh optional SpatRast or array of relative humidities (percentage). If
#' supplied, a humidity-dependent lapse rate is calculated. If not supplied, a
#' fixed lapse rate of -0.05 deg C/m is applied.
#' @param pk optional SpatRast or array of atmospheric pressures (kPa). Needed if
#' lapse-rate humidity dependent.
#' @return a multi-layer SpatRast of elevation-corrected temperatures (deg C) matching
#' the resolution, extent and crs of dtmf.
#' @rdname tempelev
#' @import terra
#' @export
#'
#' # add example
tempelev <- function(tc, dtmf, dtmc = NA, rh = NA, pk = NA) {
  if (class(dtmc)[1] == "logical")  dtmc<-resample(dtmf,tc)
  if (class(tc)[1] == "array") tc<-.rast(tc,dtmc)
  # Calculate lapse rate
  n<-dim(tc)[3]
  if (class(rh) == "logical") {
    lrc<-.rta(0.005*.is(dtmc),n)
    lrf<-.rta(0.005*.is(dtmf),n)
  } else {
    ea<-.satvap(.is(tc))*(.is(rh)/100)
    lrc<-lapserate(.is(tc), ea, .is(pk))
    lrc<-.rast(lrc,dtmc)
    if (crs(lrc) != crs(dtmf)) {
      lrcp<-project(lrc,crs(dtmf))
      lrf<-resample(lrcp,dtmf)
    } else lrf<-resample(lrc,dtmf)
    lrc<-as.array(lrc)*.rta(dtmc,n)
    lrf<-as.array(lrf)*.rta(dtmf,n)
  }
  # Sea-level temperature
  stc<-.rast(.is(tc)+lrc,dtmc)
  if (crs(dtmc) != crs(dtmf)) stc<-project(stc,crs(dtmf))
  stc<-resample(stc,dtmf)
  # Actual temperature
  tcf<-suppressWarnings(stc-.rast(lrf,dtmf))
  return(tcf)
}
#' @title Downscale temperature with cold air drainage effects
#' @description Downscales an array of temperature data adjusting for cold-air drainage
#' @param climdata a model object containing climate data of the same format as `era5climdata`
#' @param dtmf a fine-resolution SpatRast of elevations. Temperatures down-scaled
#' to resolution of `dtmf`.
#' @param basins optionally, a fine-resolution SpatRast of basins matching the
#' coordinate reference system and extent of `dtmf`. Calculated if not supplied.
#' @param dtmc optionally, height of weather measurements in `climdata` (default 2m)
#' @return a SpatRast of temperature differences owing to cold-air drainage (deg C).
#' @rdname tempcad
#' @import terra
#' @export
#' @details Cold air drainage is calaculated by delineating hydrological basins and
#' calculating flow accumulation and the elevation difference from the highest point of the basin.
#' Cold-air drainage is assumed to occur when atmospheric stability is high, namely when
#' the radiation balance is negative and wind speeds are low,
tempcad<-function(climdata, dtmf, basins = NA, refhgt = 2) {
  # Calculate elevation difference between basin height point and pixel
  if (class(basins) == "logical") basins<-basindelin(dtmf,refhgt)
  b<-.is(basins)
  d<-.is(dtmf)
  u<-unique(b)
  u<-u[is.na(u)==FALSE]
  bmx<-b*0
  for (i in 1:length(u)) {
    s<-which(b==u[i])
    mx<-max(d[s],na.rm=TRUE)
    bmx[s]<-mx
  }
  edif<-bmx-.is(dtmf)
  edif<-.rast(edif,dtmf)
  # Calculate lapse-rate multiplication factor
  mu<-edif*.cadpotential(dtmf,basins,refhgt)
  # extract climate variables
  relhum<-climdata$climarray$relhum
  pk<-climdata$climarray$pres
  dtmc<-rast(climdata$dtmc)
  # Calculate lapse rate
  ea<-.satvap(tc)*(relhum/100)
  lr<-lapserate(tc, ea, pk)
  lr<-.rast(lr,dtmc)
  if (crs(lr) != crs(dtmf)) lr<-project(lr,crs(dtmf))
  lr<-resample(lr,dtmf)
  n<-dim(lr)[3]
  cad<-.is(lr)*-.rta(mu,n)
  # determine whether cold-air drainage conditions exist
  d<-0.65*0.12
  zm<-0.1*0.12
  # Extract additional ccimate variables
  u2<-climdata$climarray$windspeed
  swrad<-climdata$climarray$swrad
  lwrad<-climdata$climarray$lwrad
  uf<-(0.4*u2)/log((refhgt-d)/zm)
  H<-(swrad+lwrad-(5.67*10^-8*0.97*(tc+273.15)^4))*0.5
  st<- -(0.4*9.81*(refhgt-d)*H)/(1241*(tc+273.15)*uf^3)
  st<-.rast(st,dtmc)
  if (crs(st) != crs(dtmf)) st<-project(st,crs(dtmf))
  st<-resample(st,dtmf)
  st<-.is(st)
  st[st>1]<-1
  st[st<1]<-0
  ce<-.rast(cad*st,dtmf)
  return(ce)
}
#' @title Downscale temperature with coastal effects
#' @description Downscales an array of temperature data adjusting for coastal effects
#' @param tc an array of hourly temperature data in deg C (either coarse-resolution,
#' matching `dtmc` or fine-scale matching `dtmf`)
#' @param SST a SpatRast of hourly sea-surface temperature data (deg C) with no NAs, as
#' returned by [SSTinterpolate()]
#' @param u2 a SpatRast of high-resolution wind speeds 2 m above ground as returned
#' by [winddownscale()]
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (see details)
#' @param dtmc a coarse-resolution SpatRast of elevations matching
#' the resolution, extent and coordinate reference system of `wspeed`.
#' @return a SpatRast of temperature differences owing to cold-air drainage (deg C).
#' @rdname tempcoastal
#' @import terra
#' @export
#' @examples
#' tc<-era5data$climarray$temp
#' tme<-as.POSIXlt(era5data$tme,tz="UTC")
#' SST<-SSTinterpolate(rast(era5sst),tme,tme)
#' wspeed<-era5data$climarray$windspeed
#' wdir<-era5data$climarray$winddir
#' dtmc<-rast(era5data$dtmc)
#' # Downscale wind - takes a few seconds
#' u2<-winddownscale(wspeed, wdir, rast(dtmf), rast(dtmm), dtmc, uz = 2)
#' # Apply coastal effects (takes a few seconds)
#' tcf<-tempcoastal(tc,SST,u2,wdir,rast(dtmf),rast(dtmm),dtmc)
#' par(mfrow=c(2,1))
#' plot(tcf[[1]])
#' plot(tcf[[12]])
tempcoastal<-function(tc, SST, u2, wdir, dtmf, dtmm, dtmc) {
  # produce land sea mask
  if (crs(dtmm) != crs(dtmf)) dtmm<-project(dtmm,crs(dtmf))
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  if (crs(SST) != crs(dtmf)) SST<-project(SST,crs(dtmf))
  dtmf[is.na(dtmf)]<-0
  dtm<-extend(dtmf,ext(dtmm))
  landsea<-resample(dtmm,dtm)
  m<-.is(dtm)
  m2<-.is(landsea)
  s<-which(is.na(m))
  m[s]<-m2[s]
  landsea<-.rast(m,landsea)
  landsea[landsea==0]<-NA
  lsr<-array(NA,dim=c(dim(dtmf)[1:2],8))
  for (i in 0:7) {
    lsr[,,i+1]<-.is(coastalexposure(landsea, ext(dtmf), i%%8*45))
  }
  # smooth
  lsr2<-lsr
  for (i in 0:7) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%8+1]+0.5*lsr[,,i%%8+1]+0.25*lsr[,,(i+1)%%8+1]
  lsm<-apply(lsr,c(1,2),mean)
  # slot in wind speeds
  wdr<-.rast(wdir,dtmf)
  ll<-.latlongfromrast(wdr)
  xy<-data.frame(x=ll$long,y=ll$lat)
  wdir<-as.numeric(xx<-extract(wdr,xy))[-1]
  if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,mean,na.rm=TRUE)
  # Calculate array land-sea ratios for every hour
  i<-round(wdir/45)%%8
  lsr<-lsr2[,,i+1]
  # Calculate SST weigthing upwind
  b1<-11.003*log(.is(u2))-9.357
  d1<-(1-(1-lsr)+2)/3
  mn<-(2/3)^b1
  rge<-1-(2/3)^b1
  wgt1<-1-(d1^b1-mn)/rge
  # Calculate SST weigthing all directions
  d2<-.mta((1-lsm+2)/3,dim(wgt1)[3])
  b2<-0.6253*log(.is(u2))-3.5185
  mn<-(2/3)^b2
  rge<-1-(2/3)^b2
  wgt2<-(d2^b2-mn)/rge
  swgt<-0.5*wgt1+0.5*wgt2
  if ((dim(tc)[1]*dim(tc)[2]) != (dim(dtmf)[1]*dim(dtmf)[2])) {
    tc<-.rast(tc,dtmc)
    tc<-resample(tc,dtmf)
  }
  if ((dim(SST)[1]*dim(SST)[2]) != (dim(dtmf)[1]*dim(dtmf)[2])) SST<-resample(SST,dtmf)
  tcf<-swgt*.is(SST)+(1-swgt)*.is(tc)
  tcf<-.rast(tcf,dtmf)
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
#'
#' # add example
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
# NB function below still needs to be modified to account for elevation effects
# Worth playing around with thresholds to simulate patchiness a bit better
#' @title Downscale shortwave radiation
#' @description Downscales an array of shortwave radiation data with option to
#' simulate cloud patchiness
#' @param swrad a SpatRast or array of shortwave radiation (W/m^2). If an array dtmc must
#' be provided.
#' @param tme POSIXlt object of times corresponding to radiation values in `swrad`.
#' @param dtmf a fine-resolution SpatRast of elevations. Temperatures down-scaled
#' to resolution of `dtmf`.
#' @param dtmc optional SpatRast of elevations matching resolution of `swrad`. If
#' not supplied, `swrad` must be a SpatRast and `dtmc` is derived by resampling `dtmf`
#' to resolution of `swrad` and extents must match.
#' @param patchiness one of 0 (not simulated), 1 (low patchiness), 2 (medium patchines)
#' or (3) high patchiness
#' @return a multi-layer SpatRast of shortwave radiation (W/m^2) matching the
#' resolution of dtmf.
#' @rdname swraddownscale
#' @import terra, gstat
#' @export
#'
#' # add example
swdownscale<-function(swrad, tme, dtmf, dtmc, patchiness = 0) {
  if (class(dtmc)[1] == "logical")  dtmc<-resample(dtmf,pk)
  # Calculate clear sky fraction
  dtmc<-project(dtmc,"+proj=longlat +datum=WGS84")
  # Compute coarse res clear-sky fraction
  lats <- .latsfromr(dtmc)
  lons <- .lonsfromr(dtmc)
  n<-length(tme)
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  csr<-clearskyrad(.vta(jd,dtmc),.vta(lt,dtmc),.rta(lats,n),.rta(lons,n))
  csf<-.is(swrad)/csr
  csf[is.na(csf)]<-0.5
  csf[is.infinite(csf)]<-0.5
  # Resample clear-sky fraction
  csf<-.rast(csf,dtmc)
  if (crs(dtmf) != crs(dtmc)) csf<-project(csf,crs(dtmf))
  csf<-resample(csf,dtmf)
  if (patchiness > 0) {
    lcsf<-log(csf/(1-csf))
    if (patchiness == 1) mu<-.simpatch(dtmf,n,mn=0.6,mx=1.666667)
    if (patchiness == 2) mu<-.simpatch(dtmf,n,mn=0.4,mx=2.5)
    if (patchiness == 3) mu<-.simpatch(dtmf,n,mn=0.2,mx=5)
    lcsf<-lcsf*mu
    csf<-1/(1+exp(-lcsf))
  }
  # Calculate fine-res clear sky rad
  xy <- data.frame(xyFromCell(dtmf, 1:ncell(dtmf)))
  xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = sf::st_crs(dtmf)$wkt)
  ll <- sf::st_transform(xy, 4326)
  xy<-data.frame(x=sf::st_coordinates(ll)[,1],y=sf::st_coordinates(ll)[,2])
  lats<-matrix(xy$y,ncol=dim(dtmf)[2],byrow=T)
  lons<-matrix(xy$x,ncol=dim(dtmf)[2],byrow=T)
  csr<-clearskyrad(.vta(jd,dtmf),.vta(lt,dtmf),.rta(lats,n),.rta(lons,n))
  csr<-.rast(csr,dtmf)
  rad<-csr*csf
  return(rad)
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
# ============================================================================ #
# ~~~~~~~~~~~~~~~~ Relative humidity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Functions would: (1) Convert coarse-res relative humidity to vapour pressure
# ~~   (temperature needed as input), bilinearly interpolate vapour pressure then
# ~~   back convert to temperature


# ========================== NB - code dump from here ======================= #
# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Precipitation downscale  ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Method probably needs to fit a thin-plate spline to e.g. monthly rainfall
# ~~   to derive a simple multiplication factor. Complexity is handing zero-inflation
# ~~   Probably something clever that can be done with Bartless-Lewis rectangular pulse
# ~~   model. Alternative is to model rain day fraction as function of elevation too
#' Calculates elevation correction factor to apply to fine-resolution rainfall data
#'
#' @description
#' `raincorrect` Calculates an altitudinal correction factor for application to fine-resolution e.g.
#' daily rainfall derived using bilinear interpolation.
#'
#' @param demc a coarse-resolution raster of digital elevation data
#' @param demf a fine-resolution raster of digital elevation data for the region for which correction factors are required
#' @param rainc a coarse-resolution raster of e.g. monthly or annual rainfall covering the same extent
#' as `demc`
#' @param rainf a fine-resolution raster of e.g. monthly or annual rainfall derived by resampling
#' coarse-resolution rainfall data, covering the same extent as `demf`.

#' @return a raster of correction factors
#'
#' @details Fine-scale e.g. daily or hourly ainfall data derived by resampling or interpolating
#' coarse-resolution data fails to adequately capture elevation effects on rainfall. This
#' function fits a thin plate spline model fitted to coarse-resolution rainfall data with
#' elevation as a predictor. The model is then applied to fine-resolution data and the results
#' compared to those obtained using simple raster resampling using bilinear interpolation.
#' The correction factor can then be applied to e.g. daily rainfall.
#'
#' @export
#'
#' @examples
#' rainf <- resample(cornwallrain, dtm100m)
#' plot(mask(rainf, dtm100m)) # rainfall derived using bilinear interpolation
#' cf <- raincorrect(dtm1km, dtm100m, cornwallrain, rainf) # takes ~ 20 seconds to run
#' plot(rainf * cf) # rainfall with correction factor applied
#'
raincorrect <- function(demc, demf, rainc, rainf) {
  xy <- data.frame(xyFromCell(demc, 1:ncell(demc)))
  z <- extract(demc, xy)
  xyz <- cbind(xy, z)
  v <- extract(rainc, xy)
  sel <- which(is.na(v) == F)
  v <- v[is.na(v) == F]
  xyz <- xyz[sel, ]
  tps <- fitTps(xyz, v, m = 2)
  xy <- data.frame(xyFromCell(demf, 1:ncell(demf)))
  z <- extract(demf, xy)
  xyz <- cbind(xy, z)
  v <- extract(demf, xy)
  sel <- which(is.na(v) == F)
  v <- v[is.na(v) == F]
  xyz <- xyz[sel, ]
  xy$z <- NA
  xy$z[sel] <- predict.Tps(tps, xyz)
  rf2 <- rasterFromXYZ(xy)
  cf <- rf2 / rainf
  cf
}
# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Downscale all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Worth writing a wrapper function to combine all of above into a single
# ~~   function
