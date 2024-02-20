#' @title delineate hydrological or cold-air drainage basins
#' @description The function `basindelin` uses a digital elevation dataset to delineate
#' hydrological basins, merging adjoining basis seperated by a low boundary if specified.
#' @param dtm a SpatRast object of elevations
#' @param boundary optional numeric value. If greater than 0, adjoining basins
#' separated by elevation differences < boundary are merged (see details.
#' @return a SpatRast of basins sequentially numbered as integers.
#' @details This function searches for the lowest grid cell in `dtm` and assigns it
#' as basin 1. All immediately adjacent pixels (in 8 directions) not previously assigned
#' are then assigned as being part of this basin if higher than the focal cell. The process is repeated
#' until no higher further cells are found. The next lowest unassigned grid cell is identified
#' and assigned as basin 2 and the process repeated until all grid cells are assigned a basin number.
#' If `boundary > 0`, edge grid cells are identified and the height difference from all
#' surrounding cells calculated. If the height difference is less than `boundary`, basins
#' are merged and the basins renumbered sequentially.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @export
#' @rdname basindelin
#' @examples
#' library(terra)
#' bsn<-basindelin(rast(dtmf))
#' plot(bsn)
basindelin<-function(dtm, boundary = 0) {
  dm<-dim(dtm)
  if (sqrt(dm[1]*dm[2]) > 250) {
    bsn<-.basindelin_big(dtm, boundary)
  } else bsn<-.basindelin(dtm, boundary)
  return(bsn)
}
#' Calculates accumulated flow
#'
#' @description
#' `flowacc` calculates accumulated flow(used to model cold air drainage)
#'
#' @param dtm a SpatRast elevations (m).
#' @param basins optionally a SpatRast of basins numbered as integers (see details).
#'
#' @return a SpatRast of accumulated flow (number of cells)
#' @details Accumulated flow is expressed in terms of number of cells. If `basins`
#' is provided, accumulated flow to any cell within a basin can only occur from
#' other cells within that basin.
#' @import terra
#' @export
#' @rdname flowacc
#' @examples
#' library(terra)
#' fa <- flowacc(rast(dtmf))
#' plot(fa, main = 'Accumulated flow')
flowacc <- function (dtm, basins = NA) {
  dm<-.is(dtm)
  fd<-.flowdir(dm)
  fa<-fd*0+1
  if (class(basins) != "logical") ba<-.is(basins)
  o<-order(dm,decreasing=T,na.last=NA)
  for (i in 1:length(o)) {
    x<-arrayInd(o[i],dim(dm))[1]
    y<-arrayInd(o[i],dim(dm))[2]
    f<-fd[x,y]
    x2<-x+(f-1)%%3-1
    y2<-y+(f-1)%/%3-1
    # If basin file provided only add flow accumulation of from same basin
    if (class(basins) != "logical") {
      b1<-b[x,y]
      b2<-ba[x2,y2]
      if (b1==b2 & x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1
    } else if (x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1
  }
  fa<-.rast(fa,dtm)
  return(fa)
}
#' @title temporally interpolate sea-surface temperature data to hourly
#' @description The function `SSTinterpolate` spatially infills missing spatial spatial
#' sea-surface temperature data and then temporally interpolates data to hourly if not
#' already hourly.
#' @param SST a SpatRast of hourly, daily or e.g. monthly sea-surface tmeperature data
#' @param tmein POSIXlt object of times corresponding to each layer of `SST`
#' @param tmeout POSIXlt object of times for which output sea-surface temperature data are needed (hourly)
#' @return a SpatRast of sea-surface temperature data wiht the same units as `SST` (usually deg C).
#' @details Missing grid cells are interpolated using zoo::na.approx. Output sea-surface
#' temperature data are temporally interpolated to daily values, so the duration of
#' `tmeout` must be at least that of `tmeout`. Hourly values are returned by replicating
#' all values within a day as sea-surface temperatures tend to fluctuate slowly.
#' @import terra, zoo
#' @export
#' @rdname SSTinterpolate
SSTinterpolate<-function(SST, tmein, tmeout) {
  if (length(tmein) != dim(SST)[3]) stop("length of tmein must match dim 3 of SST")
  if (tmeout[1] < tmein[1] | tmeout[length(tmeout)] > tmein[length(tmein)]) {
    stop ("tmeout extends beyond time start and/or end of tmein")
  }
  tostep<-as.numeric(tmeout[2]) - as.numeric(tmeout[2])
  if (tostep != 3600 | tostep != (3600*24)) stop("tmeout must be hourly or daily")
  # ==================================== #
  # Spatially interpolate missing values #
  # ==================================== #
  me<-as.vector(SST)
  n<-which(is.na(me))
  SSTn<-SST
  crs(SSTn)<-crs(SST) # prevents superflious warning
  if (length(n) > 0) {
    for (i in 1:dim(SST)[3]) {
      m<-.is(SST[[i]])
      m1<-m
      m10<-.is(resample(aggregate(dtm,10,na.rm=TRUE),dtm))
      m100<-.is(resample(aggregate(dtm,100,na.rm=TRUE),dtm))
      ma<-array(mean(as.vector(dtm),na.rm=TRUE),dim=dim(dtm)[1:2])
      s<-which(is.na(m1))
      m1[s]<-m10[s]
      s<-which(is.na(m1))
      m1[s]<-m100[s]
      s<-which(is.na(m1))
      m1[s]<-ma[s]
      d<-dim(m1)
      m1[2:(d[1]-1),2:(d[2]-1)]<-m[2:(d[1]-1),2:(d[2]-1)]
      m1a<-na.approx(m1)
      m1b<-t(na.approx(t(m1)))
      m<-(m1a+m1b)/2
      SSTn[[i]]<-suppressWarnings(.rast(m,SST))
    } # end for
  } # end if
  # ==================================== #
  # Temporally interpolate to hourly     #
  # ==================================== #
  if (n == 1) {
    SSTn<-.rta(SSTn,length(tmeout)) # replicate all values if SST is a single layer
  } else { # Temporally interpolate
    tstep<-as.numeric(tmein[2])-as.numeric(tmein[1])
    if (tstep > 86400) {  # checks whether step isn't hourly or daily'
      # Convert tmeout to daily if hourly
      if ((as.numeric(tmeout[2])-as.numeric(tmeout[1])) == 3600) {
        st<-substr(tmeout[1],1,11)
        n<-length(tmeout)/24
        tmed<-as.POSIXlt(c(0:n)*24*3600,origin=st,tz="UTC")
      } else tmed<-tmein
      # Calculate weights matrix
      # ** time difference tmein and tmeout
      wgts<-matrix(NA,nrow=length(tmed),ncol=length(tmein))
      for (i in 1:length(tmed)) wgts[i,]<-abs((as.numeric(tmed[i])-as.numeric(tmein))/3600)
      # Select only the two lowest and set others to NA
      for (i in 1:length(tmed)) {
        v<-wgts[i,]
        co<-v[order(v)][2]
        v[v>co]<-NA
        wgts[i,]<-v
      }
      # Calculate weights
      for (i in 1:length(tmed)) {
        v<-wgts[i,]
        wgts[i,]<-1-v/sum(v,na.rm=T)
      }
      wgts[is.na(wgts)]<-0
      # Convert SST to array
      a<-as.array(SSTn)
      ao<-array(NA,dim=c(dim(a)[1:2],length(tmed)))
      for (i in 1:length(tmed)) {
        ats<-.vta(wgts[i,],SSTn[[1]])*.is(SSTn)
        ao[,,i]<-apply(ats,c(1,2),sum)
      }
      ao<-.ehr(ao) # convert daily to hourly
      # Select only those values within tmeout
      tmih<-seq(as.numeric(tmein[1]),as.numeric(tmein[length(tmein)]),3600)
      tmih<-as.POSIXlt(round(tmih/3600,0)*3600,origin="1970-01-01 00:00",tz="UTC")
      s<-which(tmih>=tmeout[1] & tmoh<=tmeout[length(tmeout)])
      ao<-ao[,,s]
      SSTn<-.rast(ao,SSTn)
    } else { # if tstep<= daily or hourly
      # check whether it is daily oor hourly and stop if not
      if (tstep!=3600 | tstep!=(3600*24)) stop("tmein must have hourly, daily or >daily time increments")
      if (tstep > 3600) { # if daily
        ao<-as.array(SSTn)
        ao<-.ehr(ao)
        SSTn<-.rast(ao,SSTn)
      }
    } # end else tstep daily or hourly
  } # end temporal interpolation
  return(SSTn)
}
#' @title derive wind terrain adjustment coefficient
#' @description The function `windelev` is used to spatially downscale wind, and
#' adjusts wind speed for elevation and applies a terrain shelter coefficient for
#' a specified wind direction.
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (see details)
#' @param dtmc a coarse-resolution SpatRast of elevations usually matching
#' the resolution of climate data used for downscaling (see details)
#' @param wdir wind direction (from, decimal degrees).
#' @param uz height above ground (m) of wind speed measurement
#' @return a SpatRast of wind adjustment coefficients matching the resoltuion,
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
#' @import terra
#' @export
#' @seealso [winddownscale()]
#' @rdname windelev
#' @examples
#' wc <- windelev(rast(dtmf), rast(dtmm), rast(era5data$dtmc), wdir = 0)
#' plot(wc)
windelev <- function(dtmf, dtmm, dtmc, wdir, uz = 2) {
  # Reproject if necessary
  if (crs(dtmm) != crs(dtmf)) dtmm<-project(dtmm,crs(dtmf))
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  # This bit will be wrapped into a function - this for dtmm
  # Calculate wind adjustment 1
  dtmr<-dtmm*0+1
  dtmr[is.na(dtmr)]<-1
  wc1<-.windz(dtmm,dtmc,dtmr,wdir)
  # Calculate wind adjustment 2
  wc1<-resample(wc1,dtmf)
  wc2<-.windz(dtmf,dtmm,wc1,wdir)
  # Average
  wc<-(wc1+wc2)/2
  # Calculate average for coarse grid cell
  wcc<-resample(wc,dtmc,method="near")
  wcc[is.na(wcc)]<-mean(as.vector(wcc),na.rm=TRUE)
  wcc<-resample(wcc,wc)
  wc<-wc/wcc
  # Calculate terrain shelter coefficient
  ws<-.windcoef(dtmm, wdir, hgt = uz)  # coarse
  ws<-resample(ws,dtmf)
  ws2<-.windcoef(dtmf, wdir, hgt = uz) # fine
  ws<-.rast(pmin(.is(ws),.is(ws2)),dtmf)
  wc<-ws*wc
  wc<-suppressWarnings(mask(wc,dtmf))
  return(wc)
}
#' Calculates land to sea ratio in upwind direction
#'
#' @description The function `coastalexposure` is used to calculate an inverse
#' distance^2 weighted ratio of land to sea in a specified upwind direction.
#'
#' @param landsea A SpatRast with NAs (representing sea) or any non-NA value (representing land). The object should have a larger extent than that for which land-sea ratio values are needed, as the calculation requires land / sea coverage to be assessed upwind outside the target area.
#' @param e a terra::ext object indicating the region for which land-sea ratios are required.
#' @param wdir an optional single numeric value specifying the direction (decimal degrees) from which the wind is blowing.
#' @return a SpatRast of representing values ranging between zero
#' (all upwind pixels sea) to one (all upwind pixels land).
#' @details This function calculates a coefficient of the ratio of land to
#' sea pixels in a specified upwind direction, across all elements of a
#' SpatRast, weighted using an inverse distance squared function,
#' such that nearby pixels have a greater influence on the coefficient.
#'
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#'
#' @examples
#' library(terra)
#' e<-ext(rast(dtmf))
#' ce1 <- coastalexposure(rast(landsea), e, 180)
#' ce2 <- coastalexposure(rast(landsea), e, 270)
#' par(mfrow=c(2,1))
#' plot(ce1, main = "Land to sea weighting, southerly wind")
#' plot(ce2, main = "Land to sea weighting, westerly wind")
coastalexposure <- function(landsea, e, wdir) {
  # Calculate sample distances
  e2<-ext(landsea)
  maxdist<-sqrt(xres(landsea)*(e2$xmax-e2$xmin)+yres(landsea)*(e2$ymax-e2$ymin))
  reso<-xres(landsea)
  slr<-landsea*0+1
  slr[is.na(slr)]<-0
  s<-c(0,(8:1000)/8)^2*reso
  s<-s[s<=maxdist]
  if (e2 != e) {
    lss<-crop(slr,e,snap='out')
  } else lss<-slr
  lsm<-.is(lss)
  es<-ext(slr)
  lsr<-invls_calc(lsm,reso,e$xmin,e$ymax,s,wdir,.is(slr),es$xmin,es$xmax,es$ymin,es$ymax)
  lsr<-.rast(lsr,lss)
  return(lsr)
}

# ====================================================================== #
# ~~~~~~~~ Useful functions for processing climate data that we likely
# ~~~~~~~~ want to document.
# ====================================================================== #
#' @title convert era4 ncdf4 file to format required for model
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by [mcera5::request_era5()] to the correct formal required for subsequent modelling.
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtm a SpatRaster object of elevations covering the extent of the study area (see details)
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @return a list of the following:
#' \describe{
#'   \item{tme}{POSIXlt object of times corresponding to climate observations}
#'   \item{climarray}{a list of arrays of hourly weather variables
#'   \item{dtmc}{a coarse resolution digital elevation dataset matching the resolution of input
#'   climate data, but with a coordinate reference system and extent matching `dtm`}
#' }
#' @export
#' @details the model requires that input climate data are projected using a coordinate reference
#' system in which x and y are in metres. Since values returned by [mcera5::request_era5()]
#' are in lat long, the output data are reprojected using the coordinate reference system and
#' extent of dtm (but retain the approximate original grid resolution of the input climate data).
#' Returned climate data match the resolution, corrdinate reference system and extent of `dtmc`.
era5toclimarray <- function(ncfile, dtm, dtr_cor_fac = 1.285)  {
  t2m<-rast(ncfile,subds = "t2m") # Air temperature (K)
  d2m<-rast(ncfile,subds = "d2m") # Dewpoint temperature (K)
  sp<-rast(ncfile,subds = "sp") # Surface pressure (Pa)
  u10<-rast(ncfile,subds = "u10") # U-wind at 10m (m/s)
  v10<-rast(ncfile,subds = "v10") # V-wind at 10m (m/s)
  tp<-rast(ncfile,subds = "tp") # Total precipitation (m)
  msnlwrf<-rast(ncfile,subds = "msnlwrf")   # Mean surface net long-wave radiation flux (W/m^2)
  msdwlwrf<-rast(ncfile,subds = "msdwlwrf") # Mean surface downward long-wave radiation flux (W/m^2)
  fdir<-rast(ncfile,subds = "fdir") #  Total sky direct solar radiation at surface (W/m^2)
  ssrd<-rast(ncfile,subds = "ssrd") # Surface short-wave (solar) radiation downwards (W/m^2)
  lsm<-rast(ncfile,subds = "lsm") # Land sea mask
  # Create coarse-resolution dtm to use as template for resampling
  te<-terra::project(t2m[[1]],crs(dtm))
  agf<-res(te)[1]/res(dtm)[1]
  dtmc<-aggregate(dtm,fact=agf,fun=mean,na.rm=T)
  # Apply coastal correction to temperature data
  tmn<-.ehr(.hourtoday(as.array(t2m)-273.15,min))
  mu<-(1-as.array(lsm))*dtr_cor_fac+1
  tc<-.rast(((as.array(t2m)-273.15)-tmn)*mu+tmn,t2m)
  # Calculate vapour pressure
  ea<-.rast(.satvap(as.array(d2m)-273.15),t2m)
  # Resample all variables to match dtmc
  tc<-terra::project(tc,dtmc)
  ea<-terra::project(ea,dtmc)
  sp<-terra::project(sp,dtmc)
  u10<-terra::project(u10,dtmc)
  v10<-terra::project(v10,dtmc)
  tp<-terra::project(tp,dtmc)
  msnlwrf<-terra::project(msnlwrf,dtmc)
  msdwlwrf<-terra::project(msdwlwrf,dtmc)
  fdir<-terra::project(fdir,dtmc)
  ssrd<-terra::project(ssrd,dtmc)
  # Derive varies
  temp<-as.array(tc) # Temperature (deg c)
  relhum<-(as.array(ea)/.satvap(temp))*100 # Relative humidity (%)
  pres<-as.array(sp)/1000  # SEa-level surface pressure (kPa)
  swrad<-as.array(ssrd)/3600 # Downward shortwave radiation (W/m^2)
  difrad<-swrad-as.array(fdir)/3600 # Downward diffuse radiation (W/m^2)
  lwrad<-as.array(msdwlwrf-msnlwrf) # Downward longwave radiation
  windspeed<-sqrt(as.array(u10)^2+as.array(v10)^2)*log(67.8*2-5.42)/log(67.8*10-5.42) # Wind speed (m/s)
  winddir<-as.array((terra::atan2(u10,v10)*180/pi+180)%%360) # Wind direction (deg from N - from)
  prech<-as.array(tp)*1000
  # Save lists
  climarray<-list(temp=temp,relhum=relhum,pres=pres,swrad=swrad,difrad=difrad,
                  lwrad=lwrad,windspeed=windspeed,winddir=winddir,prech=prech)
  # Generate POSIXlt object of times
  tme<-as.POSIXlt(time(t2m),tz="UTC")
  # Output for returning
  out<-list(tme=tme,climarray=climarray,dtmc=dtmc)
  return(out)
}
# Converts between different humidity types
# h - humidity
# intype - one of relative, absolute, specific or vapour pressure
# outtype - one of relative, absolute, specific or vapour pressure
# returns humidity (Percentage for relatuive, Kg / Kg for specific, kg / m3 for absolute and kPa for bapur pressure)
converthumidity <- function (h, intype = "relative", outtype = "vapour pressure",
                             tc = 11, pk = 101.3) {
  tk <- tc + 273.15
  if (intype != "specific" & intype != "relative" & intype !=
      "absolute" & intype != "vapour pressure") {
    stop("No valid input humidity type specified")
  }
  if (outtype != "specific" & outtype != "relative" & outtype !=
      "absolute" & outtype != "vapour pressure") {
    stop("No valid output humidity type specified")
  }
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- (18.02 / 28.97) * (e0 / pk)
  if (intype == "vapour pressure") {
    hr <- (h/e0) * 100
  }
  if (intype == "specific") {
    hr <- (h/ws) * 100
  }
  if (intype == "absolute") {
    ea <- (tk * h) / 2.16679
    hr <- (ea/e0) * 100
  }
  if (intype == "relative")
    hr <- h
  if (max(hr, na.rm = T) > 100)
    warning(paste("Some relative humidity values > 100%",
                  max(hr, na.rm = T)))
  if (outtype == "specific") h <- (hr / 100) * ws
  if (outtype == "relative") h <- hr
  if (outtype == "absolute") h <- 2.16679 * (ea/tk)
  if (outtype == "vapour pressure") h <- e0 * (hr / 100)
  return(h)
}
#' Calculates clear sky radiation
#' @param jd astronomical Julian day
#' @param lt local time (decimal hours)
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param tc temperature (deg C)
#' @param rh relative humidity (percentage)
#' @param pk atmospheric pressure (kPa)
#' @return expected clear-sky radiation (W/m^2)
#' @export
clearskyrad <- function(jd, lt, lat, long, tc = 15, rh = 80, pk = 101.3) {
  sa<-.solalt(lt,lat,long,jd)*pi/180
  m<-35*sin(sa)*((1224*sin(sa)^2+1)^(-0.5))
  TrTpg<-1.021-0.084*(m*0.00949*pk+0.051)^0.5
  xx<-log(rh/100)+((17.27*tc)/(237.3+tc))
  Td<-(237.3*xx)/(17.27-xx)
  u<-exp(0.1133-log(3.78)+0.0393*Td)
  Tw<-1-0.077*(u*m)^0.3
  Ta<-0.935*m
  od<-TrTpg*Tw*Ta
  Ic<-1352.778*sin(sa)*TrTpg*Tw*Ta
  Ic[is.na(Ic)]<-0
  Ic
}
# Calculates day length
# Inputs:
# julian - astronomical julian day - as returned by julday()
# lat - latitude (decimal degrees)
# Returns daylength in decimal hours (0 if 24 hour darkness, 24 if 24 hour daylight)
daylength <- function(julian, lat) {
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 159.5) / 365.25))
  hc <- -0.01453808/(cos(lat*pi/180)*cos(declin)) -
    tan(lat * pi/180) * tan(declin)
  ha<-suppressWarnings(acos(hc)) * 180 / pi
  m <- 6.24004077 + 0.01720197 * (julian - 2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin (2 * m + 3.5932)
  sr <- (720 - 4* ha - eot) / 60
  ss <- (720 + 4* ha - eot) / 60
  dl <- ss - sr
  sel <- which(hc < -1)
  dl[sel] <- 24
  sel <- which(hc > 1)
  dl[sel] <- 0
  return(dl)
}
#' Calculate Lapse rates
#' ea = temperature (deg C)
#' ea = vapour pressure (kPa)
#' pk = atmospheric pressure (kPa)
lapserate <- function(tc, ea, pk) {
  rv<-0.622*ea/(pk-ea)
  lr<-9.8076*(1+(2501000*rv)/(287*(tc+273.15)))/
    (1003.5+(0.622*2501000^2*rv)/(287*(tc+273.15)^2))
  lr
}
#' Calculates the diffuse fraction from incoming shortwave radiation
#'
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#'
#' @param rad a vector of incoming shortwave radiation values (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}})
#' @param jd the Julian day as returned by [jday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param watts a logical value indicating  whether the units of `rad` are \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}} (TRUE) or \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} (FALSE).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param corr an optional numeric value representing a correction to account for over- or under-estimated diffuse proportions. Values > 1 will apportion a greater ammount of total radiation as diffuse than originally calculated by the formula.
#'
#' @return a vector of diffuse fractions (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}}).
#' @export
#'
#' @details
#' The method assumes the environment is snow free. Both overall cloud cover and heterogeneity in
#' cloud cover affect the diffuse fraction. Breaks in an extensive cloud deck may primarily
#' enhance the beam irradiance, whereas scattered clouds may enhance the diffuse irradiance and
#' leave the beam irradiance unaffected.  In consequence, if hourly data are available, an index
#' is applied to detect the presence of such variable/inhomogeneous clouds, based on variability
#' in radiation for each hour in question and values in the preceding and deciding hour.  If
#' hourly data are unavailable, an average variability is determined from radiation intensity.
#'
#' @examples
#' rad <- c(5:42) / 0.036 # typical values of radiation in W/m^2
#' jd <- jday(2017, 6, 21) # julian day
#' dfr <- difprop(rad, jd, 12, 50, -5)
#' plot(dfr ~ rad, type = "l", lwd = 2, xlab = "Incoming shortwave radiation",
#'      ylab = "Diffuse fraction")
difprop <- function(rad, jd, localtime, lat, long, hourly = FALSE,
                    watts = TRUE, merid = round(long / 15, 0) * 15, dst = 0,
                    corr = 1) {
  if (watts) rad <- rad * 0.0036
  sa <- solalt(localtime, lat, long, jd, merid, dst)
  alt <- sa * (pi / 180)
  k1 <- 0.83 - 0.56 * exp(- 0.06 * sa)
  si <- cos(pi / 2 - alt)
  si[si < 0] <- 0
  k <- rad / (4.87 * si)
  k[!is.finite(k)] <- 0
  k <- ifelse(k > k1, k1, k)
  k[k < 0] <- 0
  rho <- k / k1
  if (hourly) {
    rho <- c(rho[1], rho, rho[length(rho)])
    sigma3  <- 0
    for (i in 1:length(rad)) {
      sigma3[i] <- (((rho[i + 1] - rho[i]) ^ 2 + (rho[i + 1] - rho[i + 2]) ^ 2)
                    / 2) ^ 0.5
    }
  } else {
    sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho ^ 2 - 0.13 *
      exp(-1 * (((rho - 0.931) / 0.134) ^ 2) ^ 0.834)
    sigma3b <- 0.12 + 0.65 * (rho - 1.04)
    sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)
  }
  k2 <- 0.95 * k1
  d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa) / (sa + 3), 1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K ^ 2))
  d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
  alpha <- (1 / sin(alt)) ^ 0.6
  kbmax <- 0.81 ^ alpha
  kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
  dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax) / k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * sa)
  kL <- (k - 0.14) / (kX - 0.14)
  kR <- (k - kX) / 0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL ^ 2 *(1 - kL) * sigma3 ^ 1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 - kR) ^ 2 * sigma3 ^
                    0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[rad == 0] <- 0.5
  d[sa < 0] <- 1
  # apply correction
  dif_val <- rad * d
  dif_val_adj <- dif_val * corr
  d <- dif_val_adj /rad
  d[d > 1] <- 1
  d[d < 0] <- 1
  d[!is.finite(d)] <- 0.5
  d
}
