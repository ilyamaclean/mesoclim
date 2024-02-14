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


# ========================== NB - code dump from here ======================= #
# * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ * #
# ~~ * Functions to calculate coastal / water body effects
#      NB - these are some functions form microclima that will likely be useful.
#      But they definately need to be adapted.
# * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ * #
#' Calculates coastal exposure automatically
#' @import sf
#' @export
.invls.auto <- function(r, steps = 8, use.raster = T, zmin = 0, plot.progress = TRUE, tidyr = FALSE) {
  tidydems <- function(rfine, rc) {
    rfine[is.na(rfine)] <- zmin
    rc <- trim(rc)
    aggf <- floor(mean(res(rc)[1:2]) / mean(res(rfine)))
    if (aggf > 1) {
      rfine2 <- suppressWarnings(aggregate(rfine, aggf, max))
    } else rfine2 <- rfine
    rfine2 <- suppressWarnings(resample(rfine2, rc, method = 'ngb'))
    rfine2 <- crop(rfine2, extent(rfine))
    if (dim(rfine2)[1] * dim(rfine2)[2] == 1) {
      xx <- mean(is_raster(rfine), na.rm = T)
      rfine2 <- if_raster(as.matrix(xx), rfine2)
    }
    a <- array(-999, dim = dim(rfine2)[1:2])
    rc2 <- raster(a)
    extent(rc2) <- extent(rfine2)
    rc2 <- mosaic(rc, rc2, fun = min)
    rc2 <- mosaic(rc2, rfine2, fun = max)
    rc2[rc2 == zmin] <- NA
    rc2
  }
  adjust.lsr <- function(lsr, rs) {
    m <- is_raster(lsr)
    m[m < 0] <- 0
    s <- c(0, (8:10000) / 8) ^ 2 * 30
    sl <- which(s <= rs)
    mval <- length(sl) / 1453
    m2 <- m * (1 - mval) + mval
    m2
  }
  ll <- latlongfromraster(r)
  xs <- seq(0, by = 1.875, length.out = 192)
  ys <- seq(-88.542, by = 1.904129, length.out = 94)
  xc <- xs[which.min(abs(xs - ll$long%%360))]
  yc <- ys[which.min(abs(ys - ll$lat))]
  rll <- raster(matrix(0, nrow = 3, ncol = 3))
  extent(rll) <- extent(xc -  2.8125, xc +  2.8125,
                        yc - 2.856193, yc + 2.856193)
  crs(rll) <- sf::st_crs(4326)$wkt
  ress <- c(30, 90, 500, 1000, 10000)
  ress <- ress[ress > mean(res(r))]
  ress <- c(mean(res(r)), ress)
  ress <- rev(ress)
  # Create a list of dems
  cat("Downloading land sea data \n")
  dem.list <- list()
  rmet <- projectRaster(rll, crs = sf::st_crs(r)$wkt)
  dem <- get_dem(rmet, resolution = ress[1], zmin = zmin)
  dem.list[[1]] <- projectRaster(dem, crs = sf::st_crs(r)$wkt)
  dc <- ceiling(max(dim(r)) / 2)
  rres <- mean(res(r))
  for (i in 2:(length(ress)-1)) {
    d <- 50
    tst <- rres / ress[i] * dc * 2
    if (ress[i] <= mean(res(r)[1:2]))
      d <- dc
    if (tst > 50) d <- tst
    e <- extent(r)
    xy <- c((e@xmin + e@xmax) / 2, (e@ymin + e@ymax) / 2)
    e2 <- extent(c(xy[1] - d * ress[i], xy[1] + d * ress[i],
                   xy[2] - d * ress[i], xy[2] + d * ress[i]))
    rr <- raster()
    extent(rr) <- e2
    crs(rr) <- sf::st_crs(r)$wkt
    dem <- suppressWarnings(get_dem(rr, resolution = ress[i], zmin = zmin))
    dem.list[[i]] <- projectRaster(dem, crs = st_crs(r)$wkt)
  }
  if (tidyr & mean(res(r)) <= 30)
    warning("raster tidying ignored as resolution <= 30")
  if (tidyr & mean(res(r)) > 30) {
    rx <- raster(extent(r))
    res(rx) <- 30
    crs(rx) <- sf::st_crs(r)$wkt
    rfine <- get_dem(rx, resolution = 30, zmin = zmin)
    r <- tidydems(rfine, r)
  }
  dem.list[[i + 1]] <- r
  for (j in 1:i) dem.list[[j]] <- tidydems(r, dem.list[[j]])
  cat("Computing coastal exposure \n")
  lsa.array <- array(NA, dim = c(dim(r)[1:2], steps))
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    cat(paste("Direction:", direction, "degrees"), "\n")
    lsa.list <- list()
    for (i in 1:length(ress)) {
      dem <- dem.list[[i]]
      m <- is_raster(dem)
      m[m == zmin] <- NA
      dem <- if_raster(m, dem)
      lsa <- invls(dem, extent(r), direction)
      lsa[is.na(lsa)] <- 0
      lsm <- is_raster(lsa)
      if (max(lsm) > min(lsm)) {
        if (min(dim(lsm)) < 2) {
          r2 <- aggregate(r, 100)
          xx<- resample(lsa, r2, method = "ngb")
          xx <- resample(xx, r)
          xx <- mask(xx, r)
        } else {
          xx <- resample(lsa, r)
          x <- mask(xx, r)
        }
        mx <- mean(is_raster(xx), na.rm = T)
        if  (is.na(mx)) xx <- xx * 0 + mean(lsm, na.rm = T)
        lsa.list[[i]] <- xx
      } else  {
        lsa.list[[i]] <- if_raster(array(mean(lsm, na.rm = T),
                                         dim = dim(r)[1:2]), r)
      }
    }
    # find min vals
    for (i in 1:length(ress))
      lsa.list[[i]] <- adjust.lsr(lsa.list[[i]], ress[i])
    lsa <- is_raster(lsa.list[[1]])
    for (i in 2:length(ress))
      lsa <- pmin(lsa, is_raster(lsa.list[[i]]))
    lsa <- if_raster(lsa, r)
    if (use.raster) lsa <- mask(lsa, r)
    if (plot.progress) plot(lsa, main = paste("Direction:",direction))
    lsa.array[,,dct + 1] <- is_raster(lsa)
  }
  # Compute land-sea ratio of ncep grid cell
  cat("Computing mean coastal exposure of ncep grid cell \n")
  cncep <-0
  eone <- extent(xc - 1.875 / 2, xc + 1.875 / 2,
                 yc - 1.904129 / 2, yc + 1.904129 / 2)
  rllo <- crop(rll, eone)
  rone <- projectRaster(rllo, crs = sf::st_crs(r)$wkt)
  rone[is.na(rone)] <- zmin
  dem <- dem.list[[1]]
  m <- is_raster(dem)
  m[m == zmin] <- NA
  dem <- if_raster(m, dem)
  e <- extent(r)
  xy <- data.frame(x = (e@xmin + e@xmax) / 2,
                   y = (e@ymin + e@ymax) / 2)
  xy <- sf::st_as_sf(xy, coords = c("x","y"), crs = st_crs(r)$wkt)
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    lsa <- invls(dem, extent(rone), direction)
    f1 <- extract(lsa, xy)
    if (is.na(f1)) f1 <- mean(is_raster(lsa), na.rm = T)
    mncep <- mean(is_raster(lsa), na.rm = T) / f1
    cncep[dct + 1] <- mncep
  }
  cat("Adjusting coastal exposure by ncep means \n")
  for (i in 1:steps) lsa.array[,,i] <- lsa.array[,,i] * cncep[i]
  return(lsa.array)
}
#' Downloads sea-surface temperature data
#' @export
#' @import rnoaa ncdf4
.get_sst <- function(lat, long, tme) {
  sel1 <- which(tme$year == min(tme$year))
  sel2 <- which(tme$year == max(tme$year))
  dms <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  yr <-  unique(tme$year) + 1900
  if (length(yr) > 1) stop("SST intepolation only works on data from one year. Select dates in yearly blocks")
  if(yr%%4 == 0) dms[2] <- 29
  amonth <- min(tme$mon[sel1])
  amonth <- ifelse(amonth == 0, 12, amonth)
  pmonth <- max(tme$mon[sel1]) + 2
  pmonth <- ifelse(pmonth > 12, pmonth -12, pmonth)
  dm <- dms[max(tme$mon[sel2]) + 1]
  omn <- paste0(min(tme$year) + 1900, "-", min(tme$mon[sel1]) + 1, "-01 00:00:00")
  omx <- paste0(max(tme$year) + 1900, "-", max(tme$mon[sel2]) + 1, "-", dm, " 23:59:59")
  tmn <- as.POSIXlt(0, origin = omn, tz = "GMT") - dms[amonth] / 2 * 24 * 3600
  tmx <- as.POSIXlt(1, origin = omx, tz = "GMT") + dms[pmonth] / 2 * 24 * 3600
  tme2 <-  as.POSIXlt(seq(tmn, tmx, by = 'hour'))
  yrs <- unique(tme2$year) + 1900
  tmes <- as.POSIXlt(0, origin = "1900-01-01 00:00", tz = "GMT")
  sstdf <- data.frame(obs_time = tmes, sst = -99)
  for (yr in 1:length(yrs)) {
    dms[2] <- 28
    if (yrs[yr]%%4 == 0) dms[2] <- 29
    if (yrs[yr]%%100 == 0 & yrs[yr]%%400 != 0) dms[2] <- 28
    sel <- which(tme2$year + 1900 == yrs[yr])
    mths <- unique(tme2$mon[sel]) + 1
    for (mt in 1:length(mths)) {
      nc <- ersst(yrs[yr], mths[mt])
      sst <- ncvar_get(nc, 'sst')
      sst <- t(sst)
      sst <- apply(sst, 2, rev)
      sst <- na.approx(sst, na.rm = F)
      sstr <- raster(sst)
      extent(sstr) <- extent(-1, 359, -89, 89)
      long2 <- long%%360
      long2 <- ifelse(long2 > 359, long2 - 360, long2)
      s <- extract(sstr, data.frame(long2, lat))
      orn <- paste0(yrs[yr], "-", mths[mt], "-01 00:00")
      tmem <- as.POSIXlt(0 + 3600 * 24 * dms[mths[mt]] / 2, origin = orn, tz = "GMT")
      sstdo <- data.frame(obs_time = tmem, sst = s)
      sstdf <- rbind(sstdf, sstdo)
    }
  }
  sstdf <- sstdf[-1,]
  ssth <- spline(sstdf$sst~as.POSIXct(sstdf$obs_time), n = length(tme2))$y
  tmeh <- spline(sstdf$sst~as.POSIXct(sstdf$obs_time), n = length(tme2))$x
  tmeh <- as.POSIXlt(tmeh, origin = "1970-01-01 00:00", tz = "GMT")
  tme <- as.POSIXlt(tme + 0, tz = "GMT")
  sel <- which(tmeh >=min(tme) & tmeh <= max(tme))
  return(ssth[sel])
}
#' Calculates coastal effects from NCEP data
#'
#' @param landsea a raster object with NAs (representing sea) or any non-NA value and a projection system defined.
#' @param ncephourly a dataframe of hourly climate data as returned by [hourlyNCEP()].
#' @param steps an optional integer. Coastal effects are calculated in specified directions upwind. Steps defines the total number of directions used. If the default 8 is specified, coastal effects are calculated at 45º intervals.
#' @param use.raster an optional logical value indicating whether to mask the output values by `landsea`.
#' @param zmin optional assumed sea-level height. Values below this are set to zmin
#' @param plot.progress logical value indicating whether to produce plots to track progress.
#' @param tidyr logical value indicating whether to download 30m resolution digital elevation data
#' and use these data to improve the assignment of land and sea pixels in `landsea`.
#' @details coastalNCEP downloads digital elevation and varying resolutions to calculate
#' coastal effects by applying a variant of [invls()] over a greater extent then that specified by landsea, but the resulting outputs are
#' cropped to the same extent as landsea. Land temperatyres as a function of coastal exposure
#' and sea-surface temperature data are then calculated. Sea surface temperature data are
#' automatically downloaded from the NOAA.
#'
#' @return a three-dimension array of hourly temperature values for each pixel of `landsea`
#' @export
#' @import raster sp zoo rnoaa ncdf4 sf
#' @examples
#' library(raster)
#' # Download NCEP data
#' ll <- latlongfromraster(dtm100m)
#' tme <-as.POSIXlt(c(0:31) * 3600 * 24, origin = "2015-03-15", tz = "GMT")
#' ncephourly<-hourlyNCEP(NA, ll$lat, ll$long, tme)
#' aout <- coastalNCEP(dtm100m, ncephourly)
#' # Calculate mean temperature and convert to raster
#' mtemp <- if_raster(apply(aout, c(1, 2), mean), dtm100m)
#' plot(mtemp, main = "Mean temperature")
coastalNCEP <- function(landsea, ncephourly, steps = 8, use.raster = T, zmin = 0, plot.progress = T, tidyr = F) {
  bound <- function(x, mn = 0.1, mx = 1) {
    x[x < mn] <- mn
    x[x > mx] <- mx
    x
  }
  lsr1 <- .invls.auto(landsea, steps, use.raster, zmin, plot.progress, tidyr)
  lsr2 <- lsr1
  for (i in 1:steps) {
    mn <- i - 1
    mn <- ifelse(mn == 0, 8, mn)
    mx <- i + 1
    mx <- ifelse(mx >8, mx - 8, mx)
    lsr2[,,i] <- 0.25 * lsr1[,,mn]  + 0.5 * lsr1[,,i] + 0.25 * lsr1[,,mx]
  }
  lsrm <- apply(lsr1, c(1,2), mean)
  ll <- latlongfromraster(landsea)
  tme <- as.POSIXlt(ncephourly$obs_time)
  cat("Downloading sea-surface temperature data \n")
  sst <- .get_sst(ll$lat, ll$long, tme)
  wd <- round(ncephourly$winddir%%360 / (360 / 8), 0) + 1
  wd <- ifelse(wd > 8, wd - 8, wd)
  lws <- log(ncephourly$windspeed)
  b1 <- 11.003 * lws - 9.357
  b1[b1 < 0] <- 0
  b2 <- 0.6253 * lws - 3.5185
  dT <- sst - ncephourly$temperature
  p1 <- 1.12399 - 0.39985 * dT - 0.73361 * lws
  p2 <- -0.011142 + 0.01048 * dT + 0.043311 * lws
  d2 <- bound((1 - lsrm + 2) / 3)
  aout <- array(NA, dim = c(dim(landsea)[1:2], length(tme)))
  cat("Applying coastal effects \n")
  for (i in 1:length(tme)) {
    d1 <- bound((1 - lsr2[,,wd[i]] + 2) / 3)
    xx <- p1[i] * d1^b1[i] + p2[i] * d2^b2[i]
    xx[xx > 6] <- 6
    xx[xx < -6] <- -6
    pdT <- dT[i] + xx
    aout[,,i] <- (pdT - sst[i]) * (-1)
    if (plot.progress == T & i%%500 == 0) {
      plot(if_raster(aout[,,i], landsea), main = tme[i])
    }
  }
  return(aout)
}
# ~~ * Functions to calculate cold-air drainage effects
#' Cold air drainage direct from emissivity
# Below taken form microclima
.cadconditions2 <- function (em, wind, startjul, lat, long, starttime = 0, hourint = 1,
                             windthresh = 4.5, emthresh = 0.725, con = TRUE) {
  jd <- floor(c(1:length(em)) * hourint/24 - hourint/24 + startjul +
                starttime/24)

  st <- suntimes(jd, lat, long, merid = 0)
  hrs <- (c(1:length(em)) * hourint - hourint + starttime)%%24
  dn <- ifelse(hrs > (st$sunrise + 3) & hrs < st$sunset, 1, 0)  # sunrise before sunset
  dn2 <- ifelse(hrs > (st$sunrise + 3) | hrs < st$sunset, 1, 0)  # sunrise after sunset
  sel <- which(st$sunrise > st$sunset)
  dn[sel] <- dn2[sel]
  cad <- ifelse(dn < 1 & wind < windthresh & em < emthresh,
                1, 0)
  if (hourint == 1 & con) {
    cad[2] <- ifelse(cad[1] & cad[2] == 1, 1, 0)
    cad <- c(cad[1:2], cad[3:length(cad)] * cad[2:(length(cad) -
                                                     1)] * cad[1:(length(cad) - 2)])
  }
  cad
}
#' Internal function used to calculate flow direction
#' @export
.flowdir <- function(md) {
  fd <- md * 0
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  v <- c(1:length(md))
  v <- v[is.na(md) == F]
  x <- arrayInd(v, dim(md))[, 1]
  y <- arrayInd(v, dim(md))[, 2]
  for (i in 1:length(x)) {
    md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T))), 0)
  }
  fd
}
#' Delineates hydrological basins for large datasets
#'
#' @description
#' `basindelin_big` is for use with large digital elevation datasets, to
#' delineate hydrological or cold-air drainage basins.
#'
#' @param dem a raster object of elevations.
#' @param dirout an optional character vector containing a single path directory for temporarily storing tiles. Deleted after use. Tilde expansion (see [path.expand()]) is done.
#' @param trace a logical value indicating whether to plot and report on progress.
#'
#' @return a raster object with individual basins numbered sequentially as integers.
#' @import raster
#' @export
#' @seealso [basindelin()] for working with smaller datasets.
#'
#' @details
#' The function `basindelin_big` divides the large dataset into tiles and then
#' uses [basindelin()] to delineate basins for each tile before mosaicing back
#' together and merging basins along tile edges if not separated by a boundary
#' greater than 0. If `dirout` is unspecified, then a directory `basinsout` is
#' temporarily created within the working directory. If `trace` is TRUE (the
#' default) then progress is tracked during three stages: (1) the basins
#' of each tile are plotted, (2) basins after mosaicing, but prior
#' to merging are plotted and (3) on each merge iteration, the number of basins
#' to merge is printed and processed basin is plotted.
#'
#' @examples
#' library(raster)
#' basins <- basindelin_big(dtm1m)
#' plot(basins, main = "Basins")
basindelin_big <- function(dem, dirout = NA, trace = TRUE) {
  dem <- trim(dem)
  dmsx <- ceiling(dim(dem)[2] / 200) - 1
  dmsy <- ceiling(dim(dem)[1] / 200) - 1
  if (dmsx < 1 & dmsy < 1) {
    basins <- basindelin(dem)
  }
  else {
    xres <- xres(dem)
    yres <- yres(dem)
    ed <- extent(dem)
    if (is.na(dirout)) dirout <- "basinsout/"
    dir.create(dirout)
    fol <- ""
    ii <- 1
    for (i in 0:dmsx) {
      for (j in 0:dmsy) {
        xmn <- ed@xmin + i * 200 * xres
        xmx <- min((xmn + 200 * xres), ed@xmax)
        ymn <- ed@ymin + j * 200 * yres
        ymx <- min((ymn + 200 * yres), ed@ymax)
        e <- extent(c(xmn, xmx, ymn, ymx))
        r <- crop(dem, e)
        v <- getValues(r)
        if (is.na(mean(v, na.rm = TRUE)) == FALSE) {
          b <- basindelin(r)
          if (trace) {
            progress <- paste0(round(ii / ((dmsx + 1) * (dmsy + 1)) * 100, 1),
                               "%")
            plot(b, main = paste0("basin slice progress: ", progress))
          }
          fo <- paste0(dirout, "b", i, "_", j, ".tif")
          fol <- c(fol, fo)
          writeRaster(b, filename = fo, overwrite = TRUE)
        }
        ii <- ii + 1
      }
    }
    fol <- fol[fol != ""]
    basins <- raster(fol[1])
    ta <- max(getValues(basins), na.rm = TRUE)
    if (length(fol) > 1) {
      for (i in 2:length(fol)) {
        r <- raster(fol[i]) + ta
        basins <- mosaic(basins, r, fun = mean)
        ta <- max(getValues(basins), na.rm = TRUE)
      }
    }
    if (trace) {
      plot(basins, main = "Basin mosaic complete")
    }
    unlink(dirout)
    iter <- 1
    test <- 0
    basins[is.na(basins)] <- 9999
    dem[is.na(dem)] <- 9999
    while (test != 1) {
      basins <- basinsort(dem, basins)
      bmm <- getValues(basins, format = "matrix")
      lst <- as.list(c(1:max(getValues(basins), na.rm = TRUE)))
      ii <- 1
      if (dmsx > 0) {
        for (i in 1:dmsx) {
          xmn <- ed@xmin + i * 200 * xres - 1
          xmx <- min((xmn + 2 * xres), ed@xmax)
          e <- extent(c(xmn, xmx, ed@ymin, ed@ymax))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- .onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      if (dmsy > 0) {
        for (j in 1:dmsy) {
          ymn <- ed@ymin + j * 200 * yres - 1
          ymx <- ymn + 2 * yres
          e <- extent(c(ed@xmin, ed@xmax, ymn, ymx))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- .onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      lst <- lst[1:ii]
      ub <- unique(unlist(lst))
      ub <- ub[order(ub)]
      ub2 <- ub
      for (i in 1:length(ub)) {
        for (j in 1:ii) {
          v <- lst[[j]]
          tst <- which(v == ub[i])
          if (length(tst) > 0) ub2[i] <- min(ub2[i], min(v))
        }
      }
      sel <- which(ub != ub2)
      if (trace) {
        print(paste0("Basins to merge: ", length(sel)))
      }
      if (length(sel) == 0) {
        test <- 1
        if (trace) plot(basins, main = "Merge complete")
      }
      else {
        mbout <- bmm
        for (i in 1:length(sel)) {
          mbout[bmm == ub[sel[i]]] <- ub2[sel[i]]
        }
        basins <- raster(mbout, template = basins)
        if (trace) {
          plot(basins, main = paste0("Basin merge iteration: ", iter))
        }
        iter <- iter + 1
      }
    }
  }
  basins[dem == 9999] <- NA
  if (trace) plot (basins, main = "Merge complete")
  basins <- basinsort(dem, basins)
  basins
}
#' Merges adjoining basins
#'
#' @description
#' `basinmerge` merges adjoining basins if the height differences between the
#' bottom of the basin and the pour point is less than than that
#' specified by `boundary`.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#' @param boundary a single numeric value. Basins seperated by boundaries below this height are merged (should have same units as `dtm`).
#'
#' @return a raster object, two-dimensional array or matrix with basins numbered as integers.
#' @import raster
#' @export
#'
#' @details
#' If `dem` is a raster object, then a raster object is returned.
#' If the differences in height between the pour-point and bottom of the basin is
#' less than that specified by `boundary` the basin is merged with basin to which
#' water or air would pour.
#'
#' @examples
#' library(raster)
#' basins2 <- basinmerge(dtm100m, basins100m, 1)
#' par(mfrow=c(1, 2))
#' plot(basins100m, main = "Basins")
#' plot(basins2, main = "Merged basins")
basinmerge <- function(dem, basins, boundary) {
  if (all.equal(dim(basins)[1:2], dim(dem)[1:2]) == FALSE)
    stop ("basins and dem have different dimensions")
  r <- basins
  dem <- is_raster(dem)
  basins <- is_raster(basins)
  test <- F
  while (test == F) {
    mb2 <- basins
    bvars <- .basinchars(dem, basins)
    bkeep <- .bvarsrem(bvars, boundary)
    if (dim(bkeep)[1] == 0) test <- T
    for (b in 1:(dim(bkeep)[1])) {
      sel <- which(bkeep$pourpointbasin == bkeep$basin[b])
      if (length(sel) > 0) {
        bkeep$pourpointbasin[sel] <- bkeep$pourpointbasin[b]
      }
      sel <- which(basins == bkeep$basin[b])
      mb2[sel] <- bkeep$pourpointbasin[b]
    }
    u <- unique(as.vector(mb2))
    sel <- which(is.na(u) == F)
    u <- u[sel]
    for (i in 1:length(u)) {
      sel <- which(mb2 == u[i])
      basins[sel] <- i
    }
  }
  basins <- basinsort(dem, basins)
  if_raster(basins, r)
}
#' Calculates accumulated flow
#'
#' @description
#' `flowacc` is used by [pcad()] to calculate accumulated flow to each cold air drainage
#' basin
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#'
#' @return a raster object, two-dimensional array or matrix of accumulated flow.
#' @details Accumulated flow is expressed in terms of number of cells.
#' @export
#'
#' @examples
#' library(raster)
#' fa <- flowacc(dtm100m)
#' plot(fa, main = 'Accumulated flow')
flowacc <- function (dem)
{
  dm <- is_raster(dem)
  fd <- .flowdir(dm)
  fa <- fd * 0 + 1
  o <- order(dm, decreasing = T, na.last = NA)
  for (i in 1:length(o)) {
    x <- arrayInd(o[i], dim(dm))[1]
    y <- arrayInd(o[i], dim(dm))[2]
    f <- fd[x, y]
    x2 <- x + (f - 1)%%3 - 1
    y2 <- y + (f - 1)%/%3 - 1
    if (x2 > 0 & x2 < dim(dm)[1] & y2 > 0 & y2 < dim(dm)[2])
      fa[x2, y2] <- fa[x, y] + 1
  }
  if_raster(fa, dem)
}



# ~~ * Function to flexibly apply 1 or more of above
# The function below is taken from microclima with climate as a data,frame
# and also applies cold air drainage effects
# Calculates elevation effects
#' @import sf
#' @export
.eleveffects <- function(hourlydata, dem, lat, long, windthresh = 4.5,
                         emthresh = 0.78,  weather.elev, cad.effects) {
  xy <- data.frame(x = long, y = lat)
  if (weather.elev == 'ncep') {
    elevncep <- raster::extract(demworld, xy)
  } else if (weather.elev == 'era5') {
    lar <- round(lat*4,0)/4
    lor <- round(long*4,0)/4
    e <- extent(lor-0.875,lor+0.875,lar-0.875,lar+0.875)
    e5d <- get_dem(r = NA, lat, long, resolution = 10000, xdims = 10, ydims = 10)
    e5d <- projectRaster(e5d, crs = sf::st_crs(4326)$wkt)
    rte <- raster(e)
    res(rte) <- 0.25
    e5d <- resample(e5d,rte)
    elevncep <- raster:: extract(e5d, xy)
  } else elevncep <- as.numeric(weather.elev)
  if (is.na(elevncep)) {
    warnings("elevation of input weather data NA. Setting to zero")
    elevncep <- 0
  }
  xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 4326)
  xy <- sf::st_transform(xy, sf::st_crs(dem)$wkt)
  elev <- raster::extract(dem, xy)
  if (is.na(elev)) elev <- 0
  lr <- lapserate(hourlydata$temperature, hourlydata$humidity, hourlydata$pressure)
  elevt <- lr * (elev - elevncep) + hourlydata$temperature
  tme <- as.POSIXlt(hourlydata$obs_time)
  ## cad effects
  if (cad.effects) {
    jds <- julday(tme$year[1] + 1900, tme$mday[1] + 1, tme$mday[1])
    cad <- .cadconditions2(hourlydata$emissivity, hourlydata$windspeed,
                           jds, lat, long, starttime = tme$hour[1], hourint = 1,
                           windthresh = windthresh, emthresh = emthresh)
    pxls <- dim(dem)[1] * dim(dem)[2]
    if (pxls > 300 * 300) {
      basins <- basindelin_big(dem)
    } else {
      basins <- basindelin(dem)
    }
    basins <- basinmerge(dem, basins, 2)
    basins <- basinsort(dem, basins)
    fa <- flowacc(dem)
    pfa <- is_raster(fa) * 0
    td <- is_raster(fa) * 0
    bm <- is_raster(basins)
    dm <- is_raster(dem)
    for (b in 1:max(bm, na.rm = TRUE)) {
      sel <- which(bm == b)
      fao <- log(is_raster(fa)[sel])
      pfa[sel] <- fao/max(fao, na.rm = TRUE)
      ed <- max(dm[sel], na.rm = TRUE) - dm[sel]
      td[sel] <- ed
    }
    cdif <- pfa * td
    cdif <- if_raster(cdif, dem)
    cdif <- extract(cdif, xy)
    if (is.na(cdif)) cdif <- 0
    cadt <- cdif * lr * cad
  } else {
    basins <- dem*0+1
    fa <- basins
    cadt <- rep(0,length(tme))
  }
  tout <- data.frame(tref = hourlydata$temperature,
                     elev = elev, elevncep = elevncep,
                     telev = elevt, tcad = cadt, lapserate = lr)
  return(list(tout = tout, basins = basins, flowacc = fa))
}
# ============================================================================= #
# ********** NB another code dump for coastal / water body effects
.mta <- function(m,n) {
  a<-array(rep(m,n),dim=c(dim(m),n))
  a
}
#' Convert vector to array
.vta <- function(v,m) {
  va<-rep(v,each=dim(m)[1]*dim(m)[2])
  a<-array(va,dim=c(dim(m),length(v)))
  a
}
# Increase size of array
.expanda <- function(a,n) {
  x<-dim(a)[3]
  i<-floor((x/n)*c(0:(n-1))+1)
  a<-a[,,i]
  a
}
.rastt<-function(x,r) {
  xr<-terra::rast(x)
  terra::ext(xr)<-terra::ext(r)
  terra::crs(xr)<-terra::crs(r)
  xr
}
.bound <- function(x, mn = 0.1, mx = 1) {
  x[x < mn] <- mn
  x[x > mx] <- mx
  x
}
.interrast<-function(a,rco,rfn,maskr=TRUE,outarray=TRUE) {
  rs<-rast(a)
  if (class(rco)[1] == "RasterLayer") rco<-rast(rco)
  if (class(rfn)[1] == "RasterLayer") rfn<-rast(rfn)
  ext(rs)<-ext(rco)
  crs(rs)<-crs(rco)
  rfo<-terra::resample(rs,rfn, method="bilinear")
  if (maskr) rfo<-mask(rfo,rfn)
  if (outarray) rfo<-as.array(rfo)
  rfo
}
.invls <- function(r, steps = 8, use.raster = T, zmin = 0, plot.progress = TRUE, tidyr = FALSE) {
  tidydems <- function(rfine, rc) {
    rfine[is.na(rfine)] <- zmin
    rc <- trim(rc)
    aggf <- floor(mean(res(rc)[1:2]) / mean(res(rfine)))
    if (aggf > 1) {
      rfine2 <- suppressWarnings(aggregate(rfine, aggf, max))
    } else rfine2 <- rfine
    rfine2 <- suppressWarnings(resample(rfine2, rc, method = 'ngb'))
    rfine2 <- crop(rfine2, extent(rfine))
    if (dim(rfine2)[1] * dim(rfine2)[2] == 1) {
      xx <- mean(is_raster(rfine), na.rm = T)
      rfine2 <- if_raster(as.matrix(xx), rfine2)
    }
    a <- array(-999, dim = dim(rfine2)[1:2])
    rc2 <- raster(a)
    extent(rc2) <- extent(rfine2)
    rc2 <- mosaic(rc, rc2, fun = min)
    rc2 <- mosaic(rc2, rfine2, fun = max)
    rc2[rc2 == zmin] <- NA
    rc2
  }
  adjust.lsr <- function(lsr, rs) {
    m <- is_raster(lsr)
    m[m < 0] <- 0
    s <- c(0, (8:10000) / 8) ^ 2 * 30
    sl <- which(s <= rs)
    mval <- length(sl) / 1453
    m2 <- m * (1 - mval) + mval
    m2
  }
  ll <- latlongfromraster(r)
  xs <- seq(0, by = 1.875, length.out = 192)
  ys <- seq(-88.542, by = 1.904129, length.out = 94)
  xc <- xs[which.min(abs(xs - ll$long%%360))]
  yc <- ys[which.min(abs(ys - ll$lat))]
  rll <- raster(matrix(0, nrow = 3, ncol = 3))
  extent(rll) <- extent(xc -  2.8125, xc +  2.8125,
                        yc - 2.856193, yc + 2.856193)
  crs(rll) <- sf::st_crs(4326)$wkt
  ress <- c(30, 90, 500, 1000, 10000)
  ress <- ress[ress > mean(res(r))]
  ress <- c(mean(res(r)), ress)
  ress <- rev(ress)
  # Create a list of dems
  if (plot.progress) cat("Downloading land sea data \n")
  dem.list <- list()
  rmet <- projectRaster(rll, crs = sf::st_crs(r)$wkt)
  dem <- suppressMessages(get_dem(rmet, resolution = ress[1], zmin = zmin))
  dem.list[[1]] <- projectRaster(dem, crs = sf::st_crs(r)$wkt)
  dc <- ceiling(max(dim(r)) / 2)
  rres <- mean(res(r))
  for (i in 2:(length(ress)-1)) {
    d <- 50
    tst <- rres / ress[i] * dc * 2
    if (ress[i] <= mean(res(r)[1:2]))
      d <- dc
    if (tst > 50) d <- tst
    e <- extent(r)
    xy <- c((e@xmin + e@xmax) / 2, (e@ymin + e@ymax) / 2)
    e2 <- extent(c(xy[1] - d * ress[i], xy[1] + d * ress[i],
                   xy[2] - d * ress[i], xy[2] + d * ress[i]))
    rr <- raster()
    extent(rr) <- e2
    crs(rr) <- sf::st_crs(r)$wkt
    dem <- suppressWarnings(suppressMessages(get_dem(rr, resolution = ress[i], zmin = zmin)))
    dem.list[[i]] <- projectRaster(dem, crs = st_crs(r)$wkt)
  }
  if (tidyr & mean(res(r)) <= 30)
    warning("raster tidying ignored as resolution <= 30")
  if (tidyr & mean(res(r)) > 30) {
    rx <- raster(extent(r))
    res(rx) <- 30
    crs(rx) <- sf::st_crs(r)$wkt
    rfine <- suppressMessages(get_dem(rx, resolution = 30, zmin = zmin))
    r <- tidydems(rfine, r)
  }
  dem.list[[i + 1]] <- r
  for (j in 1:i) dem.list[[j]] <- tidydems(r, dem.list[[j]])
  if (plot.progress) cat("Computing coastal exposure \n")
  lsa.array <- array(NA, dim = c(dim(r)[1:2], steps))
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    if (plot.progress) cat(paste("Direction:", direction, "degrees"), "\n")
    lsa.list <- list()
    for (i in 1:length(ress)) {
      dem <- dem.list[[i]]
      m <- is_raster(dem)
      m[m == zmin] <- NA
      dem <- if_raster(m, dem)
      lsa <- invls(dem, extent(r), direction)
      lsa[is.na(lsa)] <- 0
      lsm <- is_raster(lsa)
      if (max(lsm) > min(lsm)) {
        if (min(dim(lsm)) < 2) {
          r2 <- aggregate(r, 100)
          xx<- resample(lsa, r2, method = "ngb")
          xx <- resample(xx, r)
          xx <- mask(xx, r)
        } else {
          xx <- resample(lsa, r)
          x <- mask(xx, r)
        }
        mx <- mean(is_raster(xx), na.rm = T)
        if  (is.na(mx)) xx <- xx * 0 + mean(lsm, na.rm = T)
        lsa.list[[i]] <- xx
      } else  {
        lsa.list[[i]] <- if_raster(array(mean(lsm, na.rm = T),
                                         dim = dim(r)[1:2]), r)
      }
    }
    # find min vals
    for (i in 1:length(ress))
      lsa.list[[i]] <- adjust.lsr(lsa.list[[i]], ress[i])
    lsa <- is_raster(lsa.list[[1]])
    for (i in 2:length(ress))
      lsa <- pmin(lsa, is_raster(lsa.list[[i]]))
    lsa <- if_raster(lsa, r)
    if (use.raster) lsa <- mask(lsa, r)
    if (plot.progress) plot(lsa, main = paste("Direction:",direction))
    lsa.array[,,dct + 1] <- is_raster(lsa)
  }
  return(lsa.array)
}
.coastalTpsa <- function(temp, wst, lswa, ud, demc, demf, tasl=FALSE, plot.progress = T, ploti=25) {
  if (class(demc)[1] == "RasterLayer") demc<-rast(demc)
  if (class(demf)[1] == "RasterLayer") demf<-rast(demf)
  # Calculate mean coastal exposure
  lsm<-.rastt(apply(lswa,c(1,2),mean),demf)
  # Remove altitudinal effect
  if (tasl == FALSE) {
    ca<-.mta(as.matrix(demc,wide=T),dim(temp)[3])
    temp<-temp+0.005*ca
  }
  # coarsen arrays
  lsrc<-.interrast(lswa,demf,demc,maskr=FALSE)
  lsrc<-.rastt(.expanda(lsrc,dim(temp)[3]),demc)
  lsmc<-resample(lsm,demc)
  # Temperature difference raster
  dTa<-wst-temp
  dTr<-.rastt(dTa,demc)
  a<-array(NA,dim=dim(lswa))
  options(warn=-1)
  mns<-pmin(apply(wst,3,min,na.rm=T),apply(temp,3,min,na.rm=T))
  mxs<-pmax(apply(wst,3,max,na.rm=T),apply(temp,3,max,na.rm=T))
  wstr<-.rastt(wst,demc)
  wstf<-terra::resample(wstr,.rastt(lsm,demf),method="bilinear")
  wd <- round(ud%%360 / (360 / 8), 0) + 1
  wd <- ifelse(wd > 8, wd - 8, wd)
  lsr<-.rastt(lswa[,,wd],demf)
  a<-array(NA,dim=dim(lsr))
  for (i in 1:dim(dTa)[3]) {
    # Fit thin plate spline model
    xy <- data.frame(xyFromCell(lsrc[[i]],1:ncell(demc)))
    xyz<-cbind(xy,as.vector(lsrc[[i]]),as.vector(lsmc))
    v<-as.vector(dTr[[i]])
    names(xyz)<-c("x","y","z1","z2")
    sel<-which(is.na(xyz$z1)==F & is.na(xyz$z2)==F & is.na(v)==F)
    xyz<-xyz[sel,]
    v<-v[sel]
    tps <- fields::Tps(xyz, v, m = 3)
    # Apply thin plate spline model
    xy <- data.frame(xyFromCell(lsr[[i]], 1:ncell(lsr[[i]])))
    xyz<-cbind(xy,as.vector(lsr[[i]]),as.vector(lsm))
    names(xyz)<-c("x","y","z1","z2")
    sel<-which(is.na(xyz$z1)==F & is.na(xyz$z2)==F)
    xyz<-xyz[sel,]
    xy$z <- NA
    xy$z[sel] <- fields::predict.Krig(tps, xyz)
    xy$z<-as.vector(wstf[[i]])-xy$z
    xy$z<-.bound(xy$z,mns[i],mxs[i])
    if (tasl == FALSE) xy$z<-xy$z-0.005*as.vector(demf)
    r <- rast(xy)
    if (plot.progress) {
      if (i%%ploti==0) plot(r,main=i)
    }
    a[,,i]<-as.matrix(r,wide=TRUE)
  }
  options(warn=0)
  return(a)
}
.coastalEmpa<-function(ud,uz,wst,temp,lswa) {
  lsrm <- apply(lswa, c(1,2), mean)
  wd <- round(ud%%360 / (360 / 8), 0) + 1
  wd <- ifelse(wd > 8, wd - 8, wd)
  lws <- log(uz)
  b1 <- 11.003 * lws - 9.357
  b1[b1 < 0] <- 0
  b2 <- 0.6253 * lws - 3.5185
  dT <- wst - temp
  p1 <- 1.12399 - 0.39985 * dT - 0.73361 * lws
  p2 <- -0.011142 + 0.01048 * dT + 0.043311 * lws
  d2 <- .bound((1 - .mta(lsrm,dim(temp)[3])+2)/3)
  lsr3<-lswa[,,wd]
  d1 <- .bound((1-lsr3+2)/3)
  xx<-p1*d1^b1+p2*d2^b2
  xx[xx>6]<-6
  xx[xx < -6] <- -6
  pdT<-dT+xx
  aout<-(pdT-wst)*(-1)
  return(aout)
}
.intwst<-function(wst,demc) {
  mns<-apply(wst,3,mean,na.rm=T)
  mns<-.vta(mns,as.array(demc))
  imx<-floor(log(min(dim(wst)[1:2]),2))
  imx<-log(imx,2)
  wr<-.rastt(wst,demc)
  for (i in 1:imx) {
    af<-2^i
    wc<-aggregate(wr,af,fun='mean',na.rm=T)
    wc<-resample(wc,demc)
    wc<-as.array(wc)
    sel<-which(is.na(wst))
    wst[sel]<-wc[sel]
  }
  mns<-apply(wst,c(1,2),mean,na.rm=T)
  mns<-.vta(mns,as.array(demc))
  sel<-which(is.na(wst))
  wst[sel]<-mns[sel]
  wst<-.rastt(wst,demc)
  wst<-aggregate(wst,2)
  wst<-resample(wst,demc)
  return(wst)
}
#' Compute coastal effects
#'
#' The function `waterbodyeffects` applies coastal effects using one of three user-specified
#' methods
#'
#' @param tme an POSIXlt object of times corresponding to values in `temperature`.
#' @param temperature a vector or array of hourly temperatures (deg C). Must be an array of `method > 1` (see details).
#' @param windspeed a vector or array of hourly wind speeds in m/s (see details).
#' @param winddir a vector or array of hourly wind directions (decimal degrees) in m/s (see details).
#' @param wst a vector, matrix or array of water body surface temperatures in deg C or NA (see details).
#' @param demf a RasterLayer or SpatRaster object of elevations (m) matching the spatial resolution for which outputs are needed
#' @param demc a RasterLayer or SpatRaster object of elevations (m) matching the spatial resoltuion of input temperature data. Ignored if method = 1.
#' @param method Optional integer in range 1-3 indicating which method is used to derive water body effects (see details)
#' @param tasl optional logical indicating whether values in `temperature` correpsond to sea-level or
#' include an altitudinal effects (see details).
#' @param steps an optional integer. Coastal effects are calculated in specified directions upwind. Steps defines the total number of directions used. If the default 8 is specified, coastal effects are calculated at 45º intervals.
#' @param use.raster an optional logical value indicating whether to mask the output values by `demf`.
#' @param zmin optional assumed sea-level height. Values below this are set to zmin
#' @param plot.progress logical value indicating whether to produce plots to track progress.
#' @param tidyr logical value indicating whether to download 30m resolution digital elevation data
#' and use these data to improve the assignment of land and sea pixels when calaculating coastal exposure
#'
#' @return a three dimensional array of temperatures with the first two dimensions matching x and y in `demf` and
#' the third dimension corresponding to times in `tme`
#'
#' @details if `wst = NA` water-bodies are assumed to be sea and sea-surface
#' temperature data are downloaded from the NOAA. If `method = 1`, Water body-land temperature differences are applied
#' using an empirical method calibrated against data obtained in Cornwall, UK.
#' If `method > 1` Thin-plate spline interpolation is used to derive coastal
#' effects with coastal exposure (upwind and mean in all directions)
#' as covariates). `Temperature` should be a three-dimensional array with x and y
#' dimensions matching those of `demc` Models are fitted with coarse-resolution
#' estimates of coastal exposure and applied at high-resolution using `demf`. To
#' ensure coastal effects capture waterbodies beyond the study area, digital
#' elevation data are automatically downloaded and any elevations < `zmin` are
#' assumed to correspond to waterbodies. Future versions of this function will
#' allow users to supply their own wide extend water-body data. The Thin-plate
#' spline methods are quite slow. If `method = 2` Thin-plate spline models are
#' fitted and applied at hours corresponding to daily maximum and minimum
#' temperatures only and the diurnal cycle in `temps` is applied to derive values
#' for hours in between. If `method = 3`, Thin-plate spline models are applied
#' hourly.
#'
#' The Thin-plate spline models derive estimates of land-water-body
#' temperature differences. `wst` can be supplied as a vector of hourly values
#' (assumes no spatial variability), a matrix of values (assumes no temporal
#' variability) or as a 3D array of in which case the x and y dims must match
#' those of `temperature`. However if a 3D array is supplied there is no need for
#' the intervals to be hourly: hourly values are derived by spline interpolation.
#' As land-water body temperature differences are needed, if `wst` contains values
#' for water-body pixels only, values for land-areas are derived by interpolation.
#'
#' If supplied as a 3D array the x and y dimensions of `temperature`, `windspeed` and
#' `winddir` and `sst` should match `demc`, but can also match `demf`. If `method > 1`
#' and `tasl = FALSE` elevational effects are removed by applying a fixed environmental
#' lapse rate of -5 degrees per 1000m prior to interpolation and the effect is
#' added back-in after downscaling.
waterbodyeffects <- function(tme, temperature, windspeed, winddir, wst, demf, demc = NA,
                             method = 1, tasl = FALSE, steps = 8, use.raster = TRUE,
                             zmin = 0, plot.progress = TRUE, tidyr = FALSE) {
  # Internal functions
  getmode <- function(v) {
    uniqv <- unique(round(v,0))
    uniqv <- uniqv[is.na(uniqv)==F]
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  splinf<-function(v,n) spline(v,n=n)$y
  # dimensions
  ll <- latlongfromraster(demf)
  mf<-is_raster(demf)
  mc<-is_raster(demc)
  n<-length(tme)
  # Coastal exposure
  if (plot.progress) cat("Calculating coastal exposure\n")
  lsr1 <- suppressWarnings(.invls(demf, steps, use.raster, zmin, plot.progress, tidyr))
  lsr2 <- lsr1
  for (i in 1:steps) {
    mn <- i - 1
    mn <- ifelse(mn == 0, steps, mn)
    mx <- i + 1
    mx <- ifelse(mx > steps, mx - steps, mx)
    lsr2[,,i] <- 0.25 * lsr1[,,mn]  + 0.5 * lsr1[,,i] + 0.25 * lsr1[,,mx]
  }
  if (class(demc)[1] == "RasterLayer") demc<-rast(demc)
  if (class(demf)[1] == "RasterLayer") demf<-rast(demf)
  if (plot.progress) cat("Formatting input data\n")
  # Water surface temperature data
  if (class(wst)[1] == "logical") {
    wst <- .get_sst(ll$lat, ll$long, tme)
    if (method==1) {
      wst<-.vta(wst,mf)
    } else wst<-.vta(wst,mc)
  } else if (class(wst)[1] == "numeric") {
    if (n > length(wst)) wst<-splinf(wst,n)
    if (method==1) {
      wst<-.vta(wst,mf)
    } else wst<-.vta(wst,mc)
  } else if (class(wst)[1] == "matrix") {
    if (method==1) {
      if (dim(wst)[1] == dim(demc)[1] & dim(wst)[2] == dim(demc)[2]) {
        wst<-.rastt(wst,demc)
        wst<-resample(wst,demf)
      }
      wst<-.mta(as.matrix(wst,wide=T),n)
      if (dim(wst)[1] != dim(mf)[1]) stop("wst must have same x and y dimensions as demf or demc")
      if (dim(wst)[2] != dim(mf)[2]) stop("wst must have same x and y dimensions as demf or demc")
    } else {
      if (dim(wst)[1] == dim(demf)[1] & dim(wst)[2] == dim(demf)[2]) {
        wst<-.rastt(wst,demf)
        wst<-resample(wst,demc)
      }
      wst<-.mta(as.matrix(wst,wide=T),n)
      if (dim(wst)[1] != dim(mc)[1]) stop("wst must have same x and y dimensions as demf or demc")
      if (dim(wst)[2] != dim(mc)[2]) stop("wst must have same x and y dimensions as demf or demc")
    }
  } else if (class(wst)[1] == "array") {
    if (dim(wst)[3]<n) {
      wst<-apply(wst,c(1,2),splinf,length(tme))
      wst<-aperm(wst,c(2,3,1))
    }
    if (method == 1) {
      if (dim(wst)[1] == dim(mc)[1] & dim(wst)[2] == dim(mc)[2]) {
        wst<-.interrast(wst,demc,demf)
      }
      if (dim(wst)[1] != dim(mf)[1]) stop("wst must have same x and y dimensions as demf or demc")
      if (dim(wst)[2] != dim(mf)[2]) stop("wst must have same x and y dimensions as demf or demc")
    } else {
      if (dim(wst)[1] == dim(mf)[1] & dim(wst)[2] == dim(mf)[2]) {
        wst<-.interrast(wst,demf,demc)
      }
      if (dim(wst)[1] != dim(mc)[1]) stop("wst must have same x and y dimensions as demf or demc")
      if (dim(wst)[2] != dim(mc)[2]) stop("wst must have same x and y dimensions as demf or demc")
    }
  } else stop("wst must be a vector, array or specified as NA")
  # Interpolate NAS
  sel<-which(is.na(wst))
  if (length(sel)>1) {
    if (method==1) {
      wst<-.rastt(wst,demf)
      wst<-resample(wst,demc)
      wst<-.intwst(wst,demc)
      wst<-resample(wst,demf)
      wst<-as.array(wst)
    } else {
      wst<-.intwst(wst,demc)
      wst<-as.array(wst)
    }
  }
  # Climate data
  # ** Temperature
  if (class(temperature)[1]== "numeric") {
    if (method == 1) {
      if (length(temperature)!=n) stop ("number of temperature records must match length of tme")
      temp<-.vta(temperature,mf)
    } else stop("for method > 1, temperature data must be an array")
  } else if (class(temperature)[1] == "array") {
    if (dim(temperature)[1] != dim(mc)[1]) stop("temperature must have same x and y dimensions as demc")
    if (dim(temperature)[2] != dim(mc)[2]) stop("temperature must have same x and y dimensions as demc")
    if (dim(temperature)[3] != n) stop ("number of temperature records must match length of tme")
    if (method==1) {
      temp<-.interrast(temperature,demc,demf)
    } else temp<-temperature
  } else stop("temperature data in wrong format")
  # ** Wind speed
  if (method == 1) {
    if (class(windspeed)[1]== "numeric" | class(windspeed)[1] == "integer") {
      if (length(windspeed)!=n) stop ("number of wind speed records must match length of tme")
      uz<-.vta(windspeed,mf)
    } else if (class(windspeed)[1] == "array") {
      if (dim(windspeed)[1] != dim(mc)[1]) stop("wind speed must have same x and y dimensions as demc")
      if (dim(windspeed)[2] != dim(mc)[2]) stop("wind speed must have same x and y dimensions as demc")
      if (dim(windspeed)[3] != n) stop ("number of wind speed records must match length of tme")
      uz<-.interrast(windspeed,demc,demf)
    } else stop("wind speed data in wrong format")
  }
  # ** Wind direction
  if (class(winddir)[1]== "numeric" | class(winddir)[1] == "integer") {
    ud<-winddir
  } else if (class(winddir)[1] == "array") {
    ud<-apply(winddir,3,getmode)
  } else stop("wind direction data in wrong format")
  if (length(ud)!=n) stop ("number of wind direction records must match length of tme")
  if (plot.progress) cat("Calculating coastal effects\n")
  if (method==1) {
    aout<-.coastalEmpa(ud,uz,wst,temp,lsr2)
  }
  if (method==2) {
    # determine which values are daily min and max
    th<-apply(temp,3,mean,na.rm=T)
    td<-matrix(th,ncol=24,byrow=T)
    ta<-(c(1:dim(td)[1])-1)*24
    smn<-apply(td,1,which.min)+ta
    smx<-apply(td,1,which.max)+ta
    # Compute array for minimum values
    amn<-.coastalTpsa(temp[,,smn],wst[,,smn],lsr2,ud[smn],demc,demf,tasl,plot.progress,10)
    amx<-.coastalTpsa(temp[,,smx],wst[,,smx],lsr2,ud[smx],demc,demf,tasl,plot.progress,10)
    # Calculate frac of dtr above min
    tmn<-rep(apply(td,1,min),each=24)
    tmx<-rep(apply(td,1,max),each=24)
    dtr<-tmx-tmn
    tfr<-.vta((th-tmn)/dtr,mf)
    # Apply tfrac
    amn<-.expanda(amn,n)
    amx<-.expanda(amx,n)
    adtr<-amx-amn
    aout<-(adtr*tfr)+amn
  }
  if (method==3) {
    aout<-.coastalTpsa(temp,wst,lsr2,ud,demc,demf,tasl,plot.progress,25)
  }
  return(aout)
}



# ============================================================================ #
# ~~~~~~~~~~~~~~~~ Relative humidity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Functions would: (1) Convert coarse-res relative humidity to vapour pressure
# ~~   (temperature needed as input), bilinearly interpolate vapour pressure then
# ~~   back convert to temperature

# ============================================================================ #
# ~~~~~~~~~~~~ Atmospheric pressure downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Function would take dtm as input and apply elevation adjustment
# ============================================================================ #

# ============================================================================ #
# ~~~~~~~~~~~~ Downward shortwave radiation downscale ~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Function to apply elevation effect. Suggested method would be to (i) calculate
# ~~   temporally average clear-sky radiation fraction, use thin-plate spline with
# ~~   elevation as a covariate to  derive a single correction coeffcient,  (ii)
# ~~   calculate clear-sky fraction in every hour, (iii) apply correction coefficient,
# ~~   (iv) back-calculate downward shortwave. Will be far too slow to apply Tps for every hour
# ~~ * Function to (optionally) simulate cloud patchiness using method currently in UKCP
# ~~   package
# ~~ * Function to estimate diffuse fraction probably best put here


# ============================================================================ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wind speed downscale ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ~~ * Method probably needs to consider downwind elevation and apply a standard
# ~~   logairthmic height profile. Suggest sampling at progressively large spacings
# ~~  and averaging wind profile



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
#' @importFrom rgcvpack fitTps
#' @importFrom rgcvpack predict.Tps
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
