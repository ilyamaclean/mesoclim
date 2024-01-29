# ====================================================================== #
# ~~~~~~~~ Useful functions for processing climate data that we likely
# ~~~~~~~~ want to document.
# ====================================================================== #
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
#' @param tme a POSIXlt object of times
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param tc temperature (deg C)
#' @param rh relative humidity (percentage)
#' @param pk atmospheric pressure (kPa)
#' @return expected clear-sky radiation (W/m^2)
#' @export
clearskyrad <- function(tme, lat, long, tc = 15, rh = 80, pk = 101.3) {
  jd<-.jday(tme)
  lt <- tme$hour+tme$min/60+tme$sec/3600
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
