#' @title Derives hourly temperature from daily data
#' @description Generates a plausible vector or SpatRaster of hourly temperatures
#' from daily data.
#' @param tmn a vector or SpatRaster of daily minimum temperatures (deg C)
#' @param tmx a vector or SpatRaster of daily maximum temperatures (deg C)
#' @param tme POSIXlt object of dates. Only required if `tmn` and
#' `tmx` are vectors.
#' @param lat latitude of location (decimal degrees). Only required if `tmn` and
#' `tmx` are vectors
#' @param lon longitude of location (decimal degrees). Only required if `tmn` and
#' `tmx` are vectors
#' @param srte a parameter controlling speed of decay of night time temperatures (see details)
#' @details Daytime temperatures are assumed to follow a sine curve with a peak a
#' short while after solar noon. After dusk, the temperatures are assumed to decay
#' exponentially reaching a minimum at dawn. The day in which tmx and tmn fall is
#' assumed to match UTC days. The parameter `stre` controls the speed of decay of
#' night time temperatures with time. A value of zero ensures values drop to minimum
#' at dawn the following day, but trial and error indicates in most circumstances
#' temperatures decay faster than this. The default value of 0.09 is an optimal
#' value derived using ERA5 data for western Europe, but performs reasonably well
#' globally
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' tc<-mesoclim::era5data$climarray$temp
#' tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2018-05-01", tz = "UTC")
#' # ========================================================================= #
#' # ~~~~~~~~~~~~~~~~~~~~ input provided as vector =========================== #
#' # ========================================================================= #
#' # Derive tmx and tmn from hourly
#' tch <- apply(tc, 3, mean, na.rm = TRUE)
#' tcm <- matrix(tch, ncol = 24, byrow = TRUE)
#' tmn <- apply(tcm, 1, min)
#' tmx <- apply(tcm, 1, max)
#' # Use interpolation function
#' tcp <- hourlytemp(tmn, tmx, tme, 50, -5)
#' tmeh <- as.POSIXct(mesoclim::era5data$tme)
#' # Plot results to compare
#' plot(tch ~ tmeh, type="l", ylim = c(8, 18), xlab = "Date", ylab = "Temperature")
#' par(new = TRUE)
#' plot(tcp ~ tmeh, type="l", ylim = c(8, 18), col = "red", xlab = "", ylab = "")
# as SpatRast input
#' # ========================================================================= #
#' # ~~~~~~~~~~~~~~~~~~~~ input provided as SpatRaster ======================= #
#' # ========================================================================= #
#' dtmc <- rast(mesoclim::era5data$dtmc)
#' # Derive tmx and tmn from hourly
#' tmn <- rast(hourtodayCpp(tc, "min"))
#' tmx <- rast(hourtodayCpp(tc, "max"))
#' Convert to SpatRaster
#' ext(tmn) <- ext(dtmc)
#' ext(tmx) <- ext(dtmc)
#' crs(tmn) <- crs(dtmc)
#' crs(tmx) <- crs(dtmc)
#' # Use interpolation function
#' tp <- hourlytemp(tmn, tmx, tme)
#' # Plot results to compare
#' tc <- rast(tc)
#' ext(tc) <- ext(tp)
#' crs(tc) <- crs(tc)
#' par(mfrow=c(2,1))
#' plot(tc[[i]])
#' plot(tp[[i]])
hourlytemp <- function(tmn, tmx, tme = NA, lat = NA, long = NA, srte = 0.09) {
  if (inherits(tmn, "SpatRaster")) {
    tme<-as.POSIXlt(terra::time(tmn))
  }
  if (inherits(tme, "POSIXlt")) {
    year<-tme$year+1900
    mon<-tme$mon+1
    day<-tme$mday
  } else stop("tme must be of class POSIXlt")
  if (inherits(tmn, "numeric") || inherits(tmn, "integer")) {
    if (inherits(lat, "logical")) stop("lat must be provided unless tmn and tmx are SpatRasts")
    if (inherits(long, "logical")) stop("long must be provided unless tmn and tmx are SpatRasts")
    if (length(tme) != length(tmn)) stop("length of tme must be the same as length of tmn")
    th<-hourlytempv(tmn,tmx,year,mon,day,lat,long,srte)
  } else if (inherits(tmn, "SpatRaster")) {
    ll<-.latslonsfromr(tmn)
    lats<-as.vector(t(ll$lats))
    lons<-as.vector(t(ll$lons))
    r<-tmn[[1]]
    d<-dim(tmn)
    tmn<-matrix(.is(tmn),ncol=d[3])
    tmx<-matrix(.is(tmx),ncol=d[3])
    th<-array(hourlytempm(tmn,tmx,year,mon,day,lats,lons,srte),dim=c(d[1:2],d[3]*24))
    th<-.rast(th,r)
    newtme<-as.POSIXlt(unlist(lapply(tme,FUN=function(x) x+(60*60*c(0:23)) )))
    terra::time(th)<-seconds(newtme) # for some reason must be stored as seconds for date-times!!
  } else stop("tmn and tmx must be vectors or SpatRasts")
  return(th)
}
#' @title Blends Met Office and ERA5 data to produce hourly 1km resolution temperature data
#' @description The function `blendtemp_hadukera5` ERA5 data to 1 km grid resoltuion,
#' calculates the diurnal cycle in each grid cell, and then adjust this by the
#' maximum and minimum daily temperatures in the one km met office data.
#' @param tasmin a stacked SpatRaster of haduk daily minimum temperatures (deg C)
#' @param tasmax a stacked SpatRaster of haduk daily maximum temperatures (deg C)
#' @param era5t2m a stacked SpatRaster of hourly ERA5 temperatures (deg C or K)
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @export
blendtemp_hadukera5<-function(tasmin,tasmax,era5t2m) {
  d1<-dim(tasmin)
  d2<-dim(tasmax)
  d3<-dim(era5t2m)
  if (sum(d1-d2) != 0) stop("dims of tasmin and tasmac must match")
  if (d3[3] != (d2[3]*24)) stop("Hours in era5 must match days in tasmin")
  # Reproject era5
  era5t2m<-project(era5t2m,tasmin)
  # met office dtr
  dtr<-tasmax-tasmin
  # era5 dtr and min
  era5min<-hourtodayCpp(.is(era5t2m),"min")
  era5max<-hourtodayCpp(.is(era5t2m),"max")
  era5dtr<-era5max-era5min
  era5minh<-.ehr(era5min)
  era5dtrh<-.ehr(era5dtr)
  era5frac<-(.is(era5t2m)-era5minh)/era5dtrh
  # met office hourly
  tasminh<-.ehr(.is(tasmin))
  dtrh<-.ehr(.is(dtr))
  tch<-.rast(xx<-era5frac*dtrh+tasminh,tasmin)
  return(tch)
}

# ============================================================================ #
# ~~~~~~~~~~~~~~~~ Relative humidity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ Function should spline interpolate vapour pressure and use temperature cycle
# ~~ to calculate relative humidity
# Derives an array of hourly relative humidity values from daily values.
# Inputs:
# relhum - an array of daily relative humidity values (%)
# tasmin - an array of daily minimum temperature values (deg C)
# tasmax - an array of daily maximum temperature values (deg C)
# temph - an array of hourly temperature values (deg C)
# psl - an array of daily surface-level pressure values(kPa)
# presh - an aray of hourly surface-level pressure values (kPa)
# tme - POSIXlt object of dates corresponding to radsw
# relmin - minimum possible relative humidity value
# adjust - optional logical which if TRUE ensures that, after interpolation, returned
#          hourly values, when averaged to daily, match the input
# Returns an array of hourly relative humidity values (%)
# Details:
# Owing to the strong dependence of relative humidity on temperature, the the
# resulting diurnal patterns, prior to interpolation, relative humidity is first
# converted to specific humidity using tasmin, tasmax and psl. After interpolation,
# the data are back-converted to relative humidity, using temph and presh.
hum_dailytohourly <- function(relhum, tasmin, tasmax, temph, psl, presh, tme, relmin = 2, adjust = TRUE) {
  # Convert to specific humidity
  tc<-(tasmin+tasmax)/2
  hs<-converthumidity(relhum,intype="relative",outtype="specific",tc=tc,pk=psl)
  hr<-daytohour(hs)
  # select days in year only
  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  sel<-which(tmeh$year==tme$year[2])
  hr<-hr[,,sel]
  relh<-suppressWarnings(converthumidity(hr,intype="specific",outtype="relative",tc=temph,pk=presh))

  # make consistent with daily
  if (adjust) {
    if (length(unique(tme$year)) > 1) {
      sel2<-c(2:(length(tme)-1))
    } else sel2<-c(1:length(tme))

    reld <- hourtoday(relh)
    mult <- relhum[,,sel2] / reld
    mult[is.na(mult)] <- 1
    mult <- daytohour(mult, Spline = FALSE)
    relh <- relh * mult
  }
  relh[relh>100]<-100
  relh[relh<relmin]<-relmin
  return(relh)
}
# ============================================================================ #
# ~~~~~~~~~~~~ Atmospheric pressure downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# Derives an array of hourly surface-level pressure values from an array of daily
# surface level pressure values
# Inputs:
# pres - an array of daily mean surface-level pressure values (kPa)
# tme - POSIXlt object of dates corresponding to radsw
# adjust - optional logical which if TRUE ensures that, after interpolation, returned
#          hourly values, when averaged to daily, match the input
# Returns an array of hourly surface-level pressure values (kPa)
# Details:
# note that, as in the era5 dataste, these are surface level pressures in kPa corrsponding to
# the elevation of each geographic location and not sea=level pressures. The
# conversion to millbars is mb = 10 kPa. The conversion to sea-level pressure (Psl) can be
# derived by reversing the equation Psl = pres*((293-0.0065z)/293)**5.26
pres_dailytohourly <- function(pres, tme, adjust = TRUE) {
  mn<-min(pres,na.rm=T)-1
  mx<-max(pres,na.rm=T)+1
  presh <- daytohour(pres)
  presh[presh<mn]<-mn
  presh[presh>mx]<-mx
  # make consistent with daily
  if (adjust) {
    presd <- hourtoday(presh)
    mult <- pres / presd
    mult[is.na(mult)] <- 1
    mult <- daytohour(mult, Spline = FALSE)
    presh <- presh * mult
  }
  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  yr<-tme$year[2]
  sel<-which(tmeh$year==yr)
  presh<-presh[,,sel]
  return(presh)
}
# ============================================================================ #
# ~~~~~~~~~~~~ Downward shortwave radiation downscale ~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# NB - this assumes input radiation is downward flux, not net radiation (as assumed in UKCP)
# to get from net to downward flux we need to recognise that rswnet = (1-alb)*radsw, so
# radsw = rswnet/(1-alb), where alb is white sky albedo. White-sky albedo changes as a function
# of solar angle, but in a manner dependent on ground reflectance, leaf area, leaf inclination
# angles and leaf transmittance and the ratio of diffuse and direct. There are too
# many vegetation parameter unknowns to reverse engineer, so suggest ignoring this.
# discrepancies probably quite minor expect in areas with very low cover and will be handled
# mostly by bias correction anyway


# ~~ * Need to spline interpolate clear-sky fraction (bounding by 1 and 0) and
# ~~   then calculate clear-sky radiation
# ~~ * Need to spline interpolate sky emissvity (bounding by 1 and 0) and
# ~~   then calculate longwave
# Derives an array of hourly radiation values from an array of daily radiation values
# Inputs:
# radsw - an array of daily mean radiation values (W/m**2)
# tme - POSIXlt object of dates corresponding to radsw
# clearky - optionally an array with dimensions matching radsw of daily clearsky
#           radiation as returned by clearskyraddaily(). Calculated if not supplied
# r - a terra::SpatRaster object giving the extent of radsw -
#     used for deriving the lat and long of grid cells
# adjust - optional logical which if TRUE ensures that, after interpolation, returned
#          hourly values, when averaged to daily, match the input
# Returns an array of hourly radiation values (W/m**2)
swrad_dailytohourly <- function(radsw, tme, clearsky = NA, r = r, adjust = TRUE) {
  # If NA, calaculate daily clearsky radiation
  if (class(clearsky) == "logical") clearsky <- .clearskyraddaily(tme,radsw)
  # Calculate clear sky fraction
  radf <- radsw/clearsky
  radf[radf > 1] <- 1
  radf[radf < 0] <- 0
  # Interpolate clear sky fraction to hourly
  radfh <- daytohour(radf)
  radfh[radfh > 1] <- 1
  radfh[radfh < 0] <- 0
  # Calculate hourly clear sky radiation
  lat <- latsfromr(r)
  lon <- lonsfromr(r)
  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  #jd <- julday(tmeh$year + 1900, tmeh$mon + 1, tmeh$mday)
  jd<-juldayvCpp(tme$year+1900, tme$mon+1, tme$mday)
  lt <- tmeh$hour
  lats <- .mta(lat, length(lt))
  lons <- .mta(lon, length(lt))
  jd <- .vta(jd, lat)
  lt <- .vta(lt, lat)
  csh<-clearskyrad(lt, lats, lons, jd)
  csh[is.na(csh)] <- 0
  # Calculate hourly radiation
  radh <- csh * radfh
  # Make consistent with daily
  if (adjust) {
    radd <- hourtoday(radh)
    mult <- radsw / radd
    mult[is.na(mult)] <- 1
    mult <- daytohour(mult, Spline = FALSE)
    radh <- radh * mult
  }
  radh[radh > 1352.778] <- 1352.778
  radh[radh < 0] <- 0
  radh[is.na(radh)] <- 0
  msk<-.mta(radsw[,,1],dim(radh)[3])
  msk[is.na(msk)==F]<-1
  radh<-radh*msk
  yr<-tme$year[2]
  sel<-which(tmeh$year==yr)
  radh<-radh[,,sel]
  return(radh)
}
# ============================================================================ #
# ~~~~~~~~~~~~ Downward longwave radiation downscale ~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# NB - more consistent to code this as downward longwave, but will essentially
# do the calaculations in the function below, but with temperature as an additional
' input'
# Derives an array of hourly effective sky-emissivity values
# Inputs:
# skyem - an array of daily mean sky-emissivity values values (0-1)
# tme - POSIXlt object of dates corresponding to radsw
# adjust - optional logical which if TRUE ensures that, after interpolation, returned
#          hourly values, when averaged to daily, match the input
# Returns an array of hourly sky-emissivity values values (0-1)
# Details:
#  Effective sky emissvity can be used to calaculate downward longwave radiation (Lwd).
#  The formula is Lwd = skyem * Lwu where Lwu is upward longwave radiation given
#  by Lwu=0.97*5.67*10**-8*(tc+273.15). Here tc is average surface temperature (deg C))
#  but an adequate approximation is derived if subtited by air temperature.
skyem_dailytohourly <- function(skyem, tme, adjust = TRUE) {
  skyemh <- daytohour(skyem)
  if (adjust) {
    skyemd <- hourtoday(skyemh)
    mult <- skyem / skyemd
    mult[is.na(mult)] <- 1
    mult <- daytohour(mult, Spline = FALSE)
    skyemh <- skyemh * mult
  }
  skyem[skyem>1] <- 1
  skyem[skyem < 0.2] <- 0.2
  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  yr<-tme$year[2]
  sel<-which(tmeh$year==yr)
  skyemh<-skyemh[,,sel]
  return(skyemh)
}
# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wind speed downscale ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Need to spline interpolate u and v wind vectors. We could simulate
# ~~   inter-hourly variability. Follows a Weiball distribution so quite easy I
# ~~   suspect.
# Derives arrays of hourly wind speed and direction from arrays of daily data
# Inputs:
# ws - an array of daily mean wind speed values (m/s)
# wd - an array of daily mean wind direction values (deg from N)
# tme - POSIXlt object of dates corresponding to radsw
# adjust - optional logical which if TRUE ensures that, after interpolation, returned
#          hourly values, when averaged to daily, match the input
# Returns an a list of two arrays:
# (1) hourly wind speed (m/s)
# (2) hourly wind direction (degrees from north)
# Details:
# For interpolation, u and v wind vectors are derived form wind speed andd direction
# and these are interpolated to hourly, with backward calculations then performed to
# derive wind speed and direction.
wind_dailytohourly <- function(ws, wd, tme, adjust = TRUE) {
  u<- -ws*sin(wd*pi/180)
  v<- -ws*cos(wd*pi/180)
  uh <- daytohour(u)
  vh <- daytohour(v)
  wdh <- (180+atan2(uh,vh)*(180/pi))%%360
  wsh <- sqrt(uh^2+vh^2)
  if (adjust) {
    wsd <- hourtoday(wsh)
    mult <- ws / wsd
    mult[is.na(mult)] <- 1
    mult <- daytohour(mult, Spline = FALSE)
    wsh <- wsh * mult
  }
  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  yr<-tme$year[2]
  sel<-which(tmeh$year==yr)
  wsh<-wsh[,,sel]
  wdh<-wdh[,,sel]
  return(list(wsh=wsh,wdh=wdh))
}
# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Precipitation downscale  ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Need to use Bartlett-Lewis rectangular pulse model + HyetosMinute method.
# ~~   coding is very poor in that package so scope to imporve it.
# sequence of rainfalls
# Z higher level (e.g. storm) sum
.propadjust <- function(rainseq, Z) {
  Xs <- rainseq * (Z / sum(rainseq))
  Xs
}
# runs Bartlett-Lewis until sequence of L wet days is generated
# BLest = Bartlett lewis paramaters
# dailysteps = number of values per day (i.e. 24 for hourly)
# dn = number of days
.level0 <- function(BLest, dailyvals, dn) {
  wet <- function(rainseq) {
    wd <- ifelse(sum(rainseq) > 0, 1, 0)
    wd
  }
  sr <- NA
  w <- NULL
  iter <- 0
  while (length(w) < 1) {
    sim<-BLRPM.sim(BLest[1] / dailyvals, BLest[2] / dailyvals, BLest[3] / dailyvals,
                   BLest[4] / dailyvals, BLest[5] / dailyvals, dn * dailyvals * 100,
                   dn * dailyvals * 100, 1, 0)
    hr <- t(matrix(sim$RR, nrow = dailyvals))
    dr <- apply(hr, 1, wet)
    x <- rle(dr)
    w <- which(x$lengths == dn & x$values == 1)
    if (length(w) > 0) {
      v <- length(rep.int(x$values[1:w[1]], x$lengths[1:w[1]])) - dn + 1
      sr <- hr[v:(v + dn -1),]
    }
    iter <- iter + 1
  }
  return(sr)
}
.level1 <- function(rainseq, BLest, dailyvals, dlim, maxiter) {
  d <- dlim * 2
  iter <- 0
  while (d > dlim & iter < maxiter) {
    dn <- length(rainseq)
    l0 <- .level0(BLest, dailyvals, dn)
    if (length(rainseq) > 1) {
      dr <- apply(l0, 1, sum)
    } else dr <- sum(l0)
    d <- sum(log((rainseq + 0.1) / (dr + 0.1))^2)^0.5
    iter <- iter + 1
  }
  if (d > dlim) l0<-NA
  return(l0)
}
.oneday <- function(dayrain, BLest, dailyvals, dlim, maxiter) {
  l1 <- NA
  while (is.na(sum(l1))) {
    l1 <- .level1(dayrain, BLest, dailyvals, dlim, maxiter)
  }
  l1
}
.level3 <- function(rainseq, BLest, dailyvals, dlim, maxiter) {
  l1 <- matrix(NA, nrow = length(rainseq), ncol = 24)
  r <- matrix(rainseq, nrow = 1)
  while (is.na(sum(l1))) {
    for (i in 1:dim(r)[1]) {
      r1 <- r[i,]
      r1 <- r1[is.na(r1) == F]
      st <- 1
      if (i > 1) {
        xx <- r[1:(i-1),]
        xx <- xx[is.na(xx) == F]
        st <- st + length(xx)
      }
      xx <- r[i,]
      xx <- xx[is.na(xx) == F]
      ed <- st + length(xx) - 1
      if (length(r1) == 1) {
        if (is.na(sum(l1[st:ed,]))) l1[st:ed,] <- .oneday(r1, BLest, dailyvals, dlim, maxiter)
      }  else {
        if (is.na(sum(l1[st:ed,]))) {
          l1[st:ed,] <- .level1(r1, BLest, dailyvals, dlim, maxiter)
          if (is.na(sum(l1[st:ed,]))) {
            ii <- which(r1 == min(r1[1:(length(r1) - 1)]))
            r11 <- r1[1:ii]
            r12 <- r1[(ii + 1):length(r1)]
            rr <- matrix(r[-i,], ncol = ncol(r))
            nc <- max(length(r11),length(r12))
            dm <- dim(rr)[1]
            mxc <- apply(rr, 1, function(x) length(x[is.na(x) ==F]))
            if (dm > 0) nc <- max(nc, max(mxc))
            r <- matrix(NA, nrow = dm + 2, ncol = nc)
            if (dm > 0) r[1:dm,1:mxc] <- rr[1:dm,1:mxc]
            r[dm + 1, 1:length(r11)] <- r11
            r[dm + 2, 1:length(r12)] <- r12
            l1a <- .level1(r11, BLest, dailyvals, dlim, maxiter)
            l1b <- .level1(r12, BLest, dailyvals, dlim, maxiter)
            st2 <- st + length(r11)
            l1[st:(st2 -1), ] <- l1a
            l1[st2:(st2 + length(r12) - 1),] <- l1b
          }
        }
      }
    }
  }
  l1
}
#' Estimate sub-daily rainfall from daily rainfall
#'
#' @description
#' `subdailyrain` estimate sub-daily rainfall using Bartlett-Lewis rectangular pulse rainfall model.
#'
#' @param rain a vector time-series of rainfall
#' @param BLpar a data.frame of Bartlett-Lewis parameters as returned by [findBLpar()].
#' @param dailyvals the number of rainfall values required for each day (i.e. 24 for hourly).
#'
#' @return A matrix with `length(rain)` rows and `dailyvals` columns of sub-daily rainfall.
#'
#' @export
#'
#' @details The function is based on the Bartlett-Lewis Rectangular Pulse model described by
#' Rodriguez-Iturbe (1987 & 1988). The model has six parameters (see [findBLpar()]) and is
#' characterized as a particular form of clustering process in which each cluster of rainfall events
#' (hereafter storms) consists of one or more rainfall cells being generated in the start of the
#' process. The parameters of `BLpar` governs the frequency of storms, the start and end of rainfall
#' events associated with each storms, the intensity of rainfall associated with storms variation
#' in the duration of storms, and can be used to generate data for any time-interval. Since these
#' vary seasonally, or by month, it is wise to generate sb-daily data seperately for each month using
#' different parameter estimates.
#'
#' Singificant element sof the coding have been borrowed from from the HyetosMinute package, and
#' the library must be loaded and attached, i.e. `library(HyetosMinute)' as the function calls C++ code
#' included with the package. The package is not available on CRAN and must be obtained or installed
#' directly from here: http://www.itia.ntua.gr/en/softinfo/3/.
#'
#' @references
#' Rodriguez-Iturbe I, Cox DR & Isham V (1987) Some models for rainfall based on stochastic point
#' processes. Proc. R. Soc. Lond., A 410: 269-288.
#'
#' Rodriguez-Iturbe I, Cox DR & Isham V (1988) A point process model for rainfall: Further
#' developments. Proc. R. Soc. Lond., A 417: 283-298.
#'
#' @examples
#' # =========================================== #
#' # ~~~ Generate hourly data for March 2015 ~~~ #
#' # =========================================== #
#' # ~~~~ Get paramaters for March
#' tme <- as.POSIXlt(dailyrain$obs_time)
#' marchrain <- dailyrain$precipitation[which(tme$mon + 1 == 3)]
#' BLpar <- findBLpar(marchrain) # Takes ~ 30 seconds
#' # ~~~~ Generate hourly data for March 2015
#' sel <- which(tme$mon + 1 == 3 & tme$year + 1900 == 2015)
#' march2015 <- dailyrain$precipitation[sel]
#' hourly <- subdailyrain(march2015, BLpar)
#' # ~~~~ Plots comparing hourly and daily / 24 data
#' o <- as.vector(t(matrix(rep(c(1:31), 24), nrow = 31, ncol = 24)))
#' marchhfd <- march2015[o] / 24
#' hourlyv <-as.vector(t(hourly))
#' dd <- c(1:(31 * 24)) / 24
#' plot(hourlyv ~ dd, type = "l", ylim = c(0, max(hourlyv)),
#'      xlab = "Decimal day", ylab = "Rain (mm / hr)", col = "red")
#' par(new = T)
#' plot(marchhfd ~ dd, type = "l", ylim = c(0, max(hourlyv)),
#'      xlab = "", ylab = "", col = "blue", lwd = 2)

#'
subdailyrain <- function(rain, BLest, dailyvals = 24, dlim = 0.2, maxiter = 1000, splitthreshold = 0.2, trace = TRUE) {
  rain[is.na(rain)] <- 0
  srain <- matrix(0, ncol = dailyvals, nrow = length(rain))
  wd <- ifelse(rain > splitthreshold, 1, 0)
  st <- which(diff(c(0, wd)) == 1)
  ed <- which(diff(c(wd, 0)) == -1)
  if (trace) cat(paste("Number of rainfall events:",length(ed),"\n"))
  for (i in 1:length(ed)) {
    r <- rain[st[i]:ed[i]]
    rseq <- .level3(r, BLest, dailyvals, dlim, maxiter)
    for (j in 1:length(r)) rseq[j,] <- .propadjust(rseq[j,], r[j])
    srain[st[i]:ed[i],] <- rseq
    if (trace) cat(paste("Completed rainfall event:",i,"\n"))
  }
  if (trace) cat("Processing days with rain < splitthreshold \n")
  sel <- which(rain <= splitthreshold & rain > 0)
  if (length(sel) > 0) {
    for (i in 1:length(sel)) {
      r <- rain[sel[i]]
      rseq <- .level3(r, BLest, dailyvals, dlim, maxiter)
      srain[sel[i],] <-  .propadjust(rseq, r)
    }
  }
  srain
}
plotrain <- function(daily, subdaily) {
  dailyvals <- length(subdaily) / length(daily)
  d24 <- as.vector(t(matrix(rep(daily, dailyvals), ncol = dailyvals)))
  d24 <- d24 / dailyvals
  day <- c(0:(length(d24) - 1)) / dailyvals
  sday <-  as.vector(t(subdaily))
  xs <- c(day, max(day), 0)
  ys1 <- c(d24, 0, 0)
  ys2 <- c(sday, 0, 0)
  par(mar=c(5,5,5,5))
  plot(NA, xlim = c(0,max(day)), ylim = c(0,max(sday)),
       ylab = "Rainfall (mm)", xlab = "Day", cex.axis = 2, cex.lab = 2)
  polygon(xs, ys1, col = rgb(1,0,0,0.5))
  polygon(xs, ys2, col = rgb(0,0,1,0.5))
}




# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Downscale all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Worth writing a wrapper function to combine all of above into a single
# ~~   function
