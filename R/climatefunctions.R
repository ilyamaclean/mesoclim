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
#' @keywords spatial
#' @rdname basindelin
#' @examples
#' bsn<-basindelin(terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim')))
basindelin<-function(dtm, boundary = 0) {
  dm<-dim(dtm)
  if (sqrt(dm[1]*dm[2]) > 250) {
    bsn<-.basindelin_big(dtm, boundary)
  } else bsn<-.basindelin(dtm, boundary)
  return(bsn)
}
#' Calculates accumulated flow
#' @description
#' `flowacc` calculates accumulated flow(used to model cold air drainage)
#' @param dtm a SpatRast elevations (m).
#' @param basins optionally a SpatRast of basins numbered as integers (see details).
#' @return a SpatRast of accumulated flow (number of cells)
#' @details Accumulated flow is expressed in terms of number of cells. If `basins`
#' is provided, accumulated flow to any cell within a basin can only occur from
#' other cells within that basin.
#' @import terra
#' @export
#' @rdname flowacc
#' @keywords spatial
#' @examples
#' fa <- flowacc(terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim')))
#' terra::plot(fa, main = 'Accumulated flow')
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
    if (class(basins) != "logical" & x2>0 & y2>0) {
      b1<-ba[x,y]
      b2<-ba[x2,y2]
      if(!is.na(b1) && !is.na(b2)){ # Check if both b1 and b2 valid basins
        if(b1==b2 & x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1 }
    } else if (x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1
  }
  fa<-.rast(fa,dtm)
  return(fa)
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
#' @keywords spatial
#' @examples
#' dtmf<-rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' climdata<-read_climdata(system.file('data/ukcpinput.rds',package='mesoclim'))
#' wc <- windelev(dtmf, dtmm, climdata$dtm, wdir = 270)
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
  wc1<-.resample(wc1,dtmf)
  wc2<-.windz(dtmf,dtmm,wc1,wdir)
  # Average
  wc<-(wc1+wc2)/2
  # Calculate average for coarse grid cell
  wcc<-resample(wc,dtmc,method="near")
  wcc[is.na(wcc)]<-mean(as.vector(wc),na.rm=TRUE)
  wcc<-.resample(wcc,wc)
  wc<-wc/wcc
  # Calculate terrain shelter coefficient
  ws<-.windcoef(dtmm, wdir, hgt = uz)  # coarse
  ws<-.resample(ws,dtmf)
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
#' @param landsea A SpatRast with NAs (representing sea) or any non-NA value (representing land).
#' The object should have a larger extent than that for which land-sea ratio values are needed,
#' as the calculation requires land / sea coverage to be assessed upwind outside the target area.
#' @param e a terra::ext object indicating the region for which land-sea ratios are required.
#' @param wdir an optional single numeric value specifying the direction (decimal degrees) from which the wind is blowing.
#' @return a SpatRast of representing values ranging between zero
#' (all upwind pixels sea) to one (all upwind pixels land).
#' @details This function calculates a coefficient of the ratio of land to
#' sea pixels in a specified upwind direction, across all elements of a
#' SpatRast, weighted using an inverse distance squared function,
#' such that nearby pixels have a greater influence on the coefficient.
#' @keywords spatial
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#'
#' @examples
#' climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' # landsea == e BUT shouldn't if be >>e ??
#' landsea<- terra::mask(mesoclim:::.resample(dtmm,dtmf),dtmf)
#' ce1 <- coastalexposure(landsea, terra::ext(dtmf), 45)
#' ce2 <- coastalexposure(landsea, terra::ext(dtmf), 270)
#' par(mfrow=c(2,1))
#' #plot(ce1, main = "Land to sea weighting, northeast wind")
#' #plot(ce2, main = "Land to sea weighting, westerly wind")
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
#' Performs thin-plate spline downscaling
#'
#' @description The function `Tpsdownscale` is a thin plate spline model, typically
#' with elevation as a covariate to downscale data.
#'
#' @param r a single layer SpatRast dataset to be downscaled.
#' @param dtmc a coarse resolution SpatRast of elevations matching the resolution, coordinate reference
#' system and extent of `r`.
#' @param dtmf a fine-resolution SpatRast of elevations.
#' @param method one of `normal`, `log` or `logit` (see details)
#' @param fast optional logical indicating whether to use [fields::fastTps()] (faster but
#' less accurate)
#' @return a SpatRast of `r` downscaled to match `dtmf`.
#' @details if `method = "log"` data are log-transformed prior to performing the downscale,
#' and then back-transformed. Use this method if input and output data must always be non-negative.
#' if `method = "logit"` data are logit-transformed prior to performing the downscale,
#' and then back-transformed. Use this method if input and output data must always be
#' bounded by 0 and 1. In both instances, the spacial case where input data are 0 or 1 is handled.
#' If `method = "normal"` no transformation is applied.
#' @import fields
#' @export
#' @keywords spatial
#' @examples
#' climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' rain<-mesoclim:::.rast(climdata$prec,climdata$dtm)[[1]]
#' try(rainf<-Tpsdownscale(rain, climdata$dtm, dtmf, method = "normal", fast = TRUE))
#' plot(rain,main='Input rain')
#' plot(rainf,main='Downscaled rain')
Tpsdownscale<-function(r, dtmc, dtmf, method = "normal", fast = TRUE) {
  # Extract values data.frame
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  if (crs(r) != crs(r)) r<-project(dtmc,crs(dtmf))
  xy<-data.frame(xyFromCell(dtmc, 1:ncell(dtmc)))
  z<-as.vector(extract(dtmc,xy)[,2])
  xyz <- cbind(xy,z)
  v<-as.vector(extract(r,xy)[,2])
  if (method == "log") {
    v2<-suppressWarnings(log(v))
    s<-which(v==0)
    if (length(s) > 0) v2[s]<-min(v2[-s],na.rm=TRUE)
    v<-v2
    v2<-NULL
  }
  if (method == "logit") {
    v2<-suppressWarnings(log(v/(1-v)))
    s<-which(v==0)
    s2<-which(is.infinite(v2)==FALSE)
    if (length(s) > 0) v2[s]<-min(v2[s2],na.rm=T)
    s<-which(v==1)
    if (length(s) > 0) v2[s]<-max(v2[s2],na.rm=T)
    v<-v2
    v2<-NULL
  }
  s<-which(is.na(xyz$z)==FALSE & is.na(v) == FALSE)
  xyz<-xyz[s,]
  v<-v[s]
  # Fit tps model
  if (fast) {
    tps<-suppressWarnings(fastTps(xyz,v,m=2,aRange=res(r)[1]*5))
  }  else {
    tps<-Tps(xyz,v,m=2)
  }
  # Apply Tps model
  xy<-data.frame(xyFromCell(dtmf,1:ncell(dtmf)))
  z<-as.vector(extract(dtmf,xy)[,2])
  xyz<-cbind(xy,z)
  s<-which(is.na(z)==FALSE)
  xy$z<-NA
  xy$z[s] <- predict(tps,xyz[s,])
  if (method == "log") xy$z<-exp(xy$z)
  if (method == "logit") xy$z<-1/(1+exp(-xy$z))
  rfn<-rast(xy,type="xyz")
  return(rfn)
}
#' Calculates the diffuse fraction from incoming shortwave radiation
#'
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#'
#' @param rad a vector of incoming shortwave radiation values (W/m^2)
#' @param julian the Julian day as returned by [julday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
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
#' @keywords climate
#' @examples
#' rad <- c(1:1352) # typical values of radiation in W/m^2
#' jd <- mesoclim:::.jday(as.POSIXlt("2022-06-21")) # julian day
#' dfr <- difprop(rad, jd, 12, 50, -5)
#' plot(dfr ~ rad, type = "l", lwd = 2,
#' xlab = expression(paste("Incoming shortwave radiation (", W*M^-2, ")")),
#' ylab = "Diffuse fraction")
difprop <- function(rad, julian, localtime, lat, long, hourly = FALSE,
                    merid = 0, dst = 0) {
  z <- .solalt(localtime, lat, long, julian, merid, dst)
  k1 <- 0.83 - 0.56 * exp(- 0.06 * (90 - z))
  si <- cos(z * pi / 180)
  si[si < 0] <- 0
  k <- rad / (1352 * si)
  k[is.na(k)] <- 0
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
  d1 <- ifelse(z < 88.6, 0.07 + 0.046 * z / (93 - z), 1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K ^ 2))
  d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
  alpha <- (1 / cos(z * pi / 180)) ^ 0.6
  kbmax <- 0.81 ^ alpha
  kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
  dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax) / k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * (90 - z))
  kL <- (k - 0.14) / (kX - 0.14)
  kR <- (k - kX) / 0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL ^ 2 *(1 - kL) * sigma3 ^ 1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 - kR) ^ 2 * sigma3 ^
                    0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[rad == 0] <- 0.5
  d[z > 90] <- 1
  # apply correction
  dif_val <- rad * d
  d <- dif_val /rad
  d[d > 1] <- 1
  d[d < 0] <- 1
  d[is.na(d)] <- 0.5
  d
}
# ====================================================================== #
# ~~~~~~~~ Useful functions for processing climate data that we likely
# ~~~~~~~~ want to document.
# ====================================================================== #

#' @title Converts between different humidity types
#' @param h- humidity
#' @param intype - one of relative, absolute, specific or vapour pressure
#' @param outtype - one of relative, absolute, specific or vapour pressure
#' @param tc - temperature
#' @param pk - surface pressure
#'
#' @return returns humidity
#' (Percentage for relative,Kg / Kg for specific, kg / m3 for absolute and kPa for vapour pressure)
#' @export
#' @keywords climate
#' @examples
#' rh<-c(25,50,75,90,100)
#' vp<-round(converthumidity(rh),3)
#' print(paste(rh,' relative humidity converts to',vp,'vapour pressure'))
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
#' @title Calculates clear sky radiation
#' @param jd astronomical Julian day
#' @param lt local time (decimal hours)
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param tc temperature (deg C)
#' @param rh relative humidity (percentage)
#' @param pk atmospheric pressure (kPa)
#' @return expected clear-sky radiation (W/m^2)
#' @export
#' @keywords climate
#' @examples
#' jd <- mesoclim:::.jday(as.POSIXlt("2022-06-21")) # julian day
#' tme<-seq(0,23,1)
#' csr<-clearskyrad(jd,tme,50,-2.5)
#' plot(csr ~ tme, type = "l", lwd = 2, xlab = expression(paste("Hour")), ylab = "Clearsky radiation")
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

#' @title Calculates day length
#' @param julian - astronomical julian day - as returned by .jday()
#' @param lat - latitude (decimal degrees)
#' @return Returns daylength in decimal hours (0 if 24 hour darkness, 24 if 24 hour daylight)
#' @export
#' @keywords climate temporal
#' @examples
#' tme<-as.POSIXlt(seq(as.POSIXlt("2022-01-01"),as.POSIXlt("2022-06-30"),60*60*24*8))
#' jd<-sapply(tme,mesoclim:::.jday)
#' dl<-daylength(jd,50)
#' plot(dl ~ jd, type = "l", lwd = 8, xlab = "Day", ylab = "Day length")
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
#' @title Calculate Lapse rates
#' @param ea = temperature (deg C)
#' @param ea = vapour pressure (kPa)
#' @param pk = atmospheric pressure (kPa)
#' @returns lapes rate
#' @export
#' @keywords climate spatial
#' @examples
#' climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
#' ea<-converthumidity(climdata$relhum,tc=climdata$temp , pk=climdata$pres)
#' lr<-lapserate(climdata$temp,ea,climdata$pres)
#' terra::plot(mesoclim:::.rast(lr,climdata$dtm)[[1]])
lapserate <- function(tc, ea, pk) {
  rv<-0.622*ea/(pk-ea)
  lr<-9.8076*(1+(2501000*rv)/(287*(tc+273.15)))/
    (1003.5+(0.622*2501000^2*rv)/(287*(tc+273.15)^2))
  lr
}

