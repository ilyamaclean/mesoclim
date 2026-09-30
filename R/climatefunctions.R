#' @title Calculate wind coefficients - for use in wind downscaling
#' @param dtmc - coarse resolution raster overing wider extent than dtmf
#' @param dtmm - medium resolution raster covering wider extent than dtmf
#' @param dtmf - fine scale raster of downscaling area
#' @param zo - wind height in metres of coeeficients
#' @param toArray - if FALSE will produce a spatraster stack otherwise as 3D array
#' @return wca2 - array of wind coeeficients for each of 8 wind directions
#' @export
calculate_windcoeffs<-function(dtmc,dtmm,dtmf,zi=10,zo=2,toArray=TRUE){
  if(all(terra::res(dtmm)==terra::res(dtmf))){
    dtmm_res<-round(exp( ( log(terra::res(dtmc)[1]) + log(terra::res(dtmf)[1]) ) / 2 ))
    dtmw<-terra::aggregate(dtmm,dtmm_res / res(dtmf),  na.rm=TRUE)
  } else if(any(terra::res(dtmm)<terra::res(dtmf))) stop("dtmm in calculate windcoeffs must be same or coarser resolution than dtmf!!") else dtmw<-dtmm
  # Calculate terrain adjustment coefs in each of 8 directions for output wind height zo
  wca<-array(NA,dim=c(dim(dtmf)[1:2],8))
  for (i in 0:7) wca[,,i+1]<-.is(windelev(dtmf,dtmw,dtmc,i*45,zi,zo))
  # smooth results
  wca2<-wca
  for (i in 0:7) wca2[,,i+1]<-0.25*wca[,,(i-1)%%8+1]+0.5*wca[,,i%%8+1]+0.25*wca[,,(i+1)%%8+1]
  if(!toArray){
    wca2<-.rast(wca2,dtmf)
    names(wca2)<-paste("Shelter coef dir",seq(1:nlyr(wca2)))
  }
  return(wca2)
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
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' wc <- windelev(dtmf, dtmm, climdata$dtm, wdir = 270)
#' terra::plot(wc)
windelev <- function(dtmf, dtmm, dtmc, wdir, zi = 10, zo=2) {
  # Reproject if necessary
  if (crs(dtmm) != crs(dtmf)) dtmm<-project(dtmm,crs(dtmf))
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  # This bit will be wrapped into a function - this for dtmm
  # Calculate wind adjustment 1
  dtmr<-dtmm*0+1
  dtmr[is.na(dtmr)]<-1
  wc1<-.windz(dtmm,dtmc,dtmr,wdir,zi,zo)
  # Calculate wind adjustment 2
  wc1<-.resample(wc1,dtmf)
  wc2<-.windz(dtmf,dtmm,wc1,wdir,zi,zo)
  # Average
  wc<-(wc1+wc2)/2
  # Calculate average for coarse grid cell
  wcc<-resample(wc,dtmc,method="near")
  wcc[is.na(wcc)]<-mean(as.vector(wc),na.rm=TRUE)
  wcc<-.resample(wcc,wc)
  wc<-wc/wcc
  # Calculate terrain shelter coefficient
  ws<-.windcoef(dtmm, wdir, hgt = zo)  # coarse
  ws<-.resample(ws,dtmf)
  ws2<-.windcoef(dtmf, wdir, hgt = zo) # fine
  ws<-.rast(pmin(.is(ws),.is(ws2)),dtmf)
  wc<-ws*wc
  wc<-suppressWarnings(mask(wc,dtmf))
  return(wc)
}

#' @title delineate hydrological or cold-air drainage basins
#' @description The function `basindelin` uses a digital elevation dataset to delineate
#' hydrological basins, merging adjoining basins separated by a low boundary if specified.
#' @param dtm a SpatRast object of elevations
#' @param boundary optional numeric value. If greater than 0, adjoining basins whose
#' lowest crossing point (pour point) is within `boundary` metres of the lower basin's
#' floor elevation are merged. Merges are applied transitively via union-find in a
#' single C++ pass. Default 0 (no merging).
#' @param method Character. Controls how a strictly higher neighbouring cell is
#' absorbed into a basin. `"any"` (default): a higher cell joins a basin as soon as
#' any already-claimed cell of that basin is adjacent to it — the original behaviour.
#' `"steepest"`: a higher cell only joins the basin that owns its single steepest
#' downhill neighbour, tying basin membership to local flow direction. Cells at
#' exactly the same elevation always merge regardless of method. `"steepest"` is more
#' physically meaningful for cold-air drainage on complex terrain; `"any"` preserves
#' backwards compatibility with earlier versions.
#' @return a SpatRast of basins sequentially numbered as integers.
#' @details Basin delineation uses a min-heap priority queue (O(N log N)), seeding each
#' new basin at the globally lowest unclaimed cell and growing it to exhaustion before
#' the next basin is considered. If `boundary > 0`, adjacent basins are merged using a
#' pour-point algorithm: for each pair of basins sharing a boundary, the lowest possible
#' crossing (pour point) is found; basins are merged when that crossing is within
#' `boundary` metres of the lower basin's own floor elevation. All merges are applied
#' transitively in a single pass.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @export
#' @keywords spatial
#' @rdname basindelin
#' @examples
#' bsn<-basindelin(terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim')))
basindelin<-function(dtm, boundary = 0, method = "any") {
  method <- match.arg(method, c("any","steepest"))
  dm<-dim(dtm)
  if (sqrt(dm[1]*dm[2]) > 250) {
    bsn<-.basindelin_big(dtm, boundary, method=method)
  } else bsn<-.basindelin(dtm, boundary, method=method)
  return(bsn)
}
#' @title Calculates accumulated flow
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
#' r<-terra::rast(matrix(c(25,24,21,24,25,15,14,11,14,15,5,4,1,4,5,15,14,11,14,15,25,24,21,24,25),nrow=5, ncol=5))
#' terra::plot(log(flowacc(r)))
flowacc<-function (dtm, basins = NA) {
  dm <- .is(dtm)
  fd <- .flowdir(dm)
  fa <- fd * 0 + 1
  if (class(basins) != "logical")
    ba <- .is(basins)
  o <- order(dm, decreasing = T, na.last = NA)
  for (i in 1:(length(o)-1)) {
    x<-arrayInd(o[i],dim(dm))[2]
    y<-arrayInd(o[i],dim(dm))[1]
    f<-fd[y,x]
    y2<-y+(f-1)%%3-1
    x2<-x+(f-1)%/%3-1
    if (class(basins) != "logical" & x2 > 0 & y2 > 0 & x2 <=
        dim(dm)[2] & y2 <= dim(dm)[1]) {
      b1 <- ba[y, x]
      b2 <- ba[y2, x2]
      if (!is.na(b1) && !is.na(b2)) {
        if (b1 == b2 & x2 > 0 & x2 < dim(dm)[2] & y2 >
            0 & y2 < dim(dm)[1])
          fa[y2,x2]<-fa[y2,x2]+fa[y,x]
      }
    }
    else if (x2 > 0 & x2 < dim(dm)[2] & y2 > 0 & y2 < dim(dm)[1])
      fa[y2,x2]<-fa[y2,x2]+fa[y,x]
  }
  fa <- .rast(fa, dtm)
  return(fa)
}

#' @title Calculates land fraction in the upwind direction
#'
#' @description The function `coastalexposure` calculates, for each land cell, the
#' proportion of land along a line in the upwind direction, with nearby land and sea
#' carrying more weight than distant land and sea.
#'
#' @param landsea a SpatRast with a projected coordinate reference system, with NA
#' representing sea and any non-NA value representing land. Output covers its whole extent.
#' @param wdir direction (decimal degrees clockwise from north) from which the wind is
#' blowing, or `"all"` for the mean of 8 directions at 45 degree intervals.
#' @param coarse optionally, a SpatRast or list of SpatRasts (NA = sea) with the same
#' coordinate reference system as `landsea`, typically at coarser resolution and covering
#' a wider area, used for samples beyond the extent of `landsea` (see details). Rasters in a
#' different coordinate reference system are projected to that of `landsea`.
#' @param n positive numeric controlling how strongly nearby land and sea is weighted
#' relative to distant land and sea (see details). Default 2.
#' @param jitter logical. If TRUE, the azimuth of each sample is perturbed by up to
#' +/-10 degrees, which removes stripes caused by samples aligning with grid rows or columns.
#'
#' @details Ported from `coastalexposure()` in the terravars package, but returning land
#' rather than sea fraction. Land and sea are sampled at distances of k^n x resolution
#' (k in steps of 1/8) along the upwind line, out to the diagonal of the largest supplied
#' raster, and the result is the mean of the samples. There is no explicit distance
#' weight: nearby cells carry more weight because samples are more closely spaced near
#' the focal cell. Each sample is read from the finest of `landsea` and `coarse` that
#' covers it, so a wide coarse raster extends the search far beyond `landsea`. Samples
#' outside all supplied rasters are ignored. When `n < 2` the search distance is capped
#' at twice the diagonal of `landsea`, with a warning if this shortens it. Jitter offsets
#' are deterministic, so results are reproducible.
#'
#' @return a SpatRast of the proportion of land upwind, from 0 (all sea) to 1 (all
#' land). Sea cells are NA.
#' @keywords spatial
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#' @examples
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' ce1 <- coastalexposure(dtmf, 45, coarse = dtmm)
#' ce2 <- coastalexposure(dtmf, 270, coarse = dtmm)
#' terra::plot(c(ce1, ce2), main = c("Land fraction, northeast wind", "Land fraction, westerly wind"))
coastalexposure <- function(landsea, wdir, coarse = NULL, n = 2, jitter = TRUE) {
  if (inherits(landsea, "PackedSpatRaster")) landsea <- terra::unwrap(landsea)
  stopifnot(inherits(landsea, "SpatRaster"))
  if (terra::is.lonlat(terra::crs(landsea))) stop("landsea must have a projected coordinate reference system!!")
  stopifnot(is.logical(jitter), length(jitter) == 1L, !is.na(jitter))
  stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n > 0)
  if (is.character(wdir)) {
    wdir <- match.arg(wdir, "all")
    dirs <- seq(0, 315, by = 45)
  } else {
    stopifnot(is.numeric(wdir), length(wdir) == 1L)
    dirs <- wdir
  }
  if (is.null(coarse)) {
    coarse <- list()
  } else if (inherits(coarse, c("SpatRaster", "PackedSpatRaster"))) {
    coarse <- list(coarse)
  }
  coarse <- lapply(coarse, function(r) if (inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r)
  if (!all(vapply(coarse, inherits, logical(1), "SpatRaster"))) stop("coarse must be a SpatRaster or a list of SpatRasters!!")

  # Levels sorted finest first: C++ uses the first level covering each sample point
  coarse <- lapply(coarse, function(r) if (terra::same.crs(r, landsea)) r else terra::project(r, terra::crs(landsea), method = "near"))
  levels <- c(list(landsea), coarse)
  resos <- vapply(levels, function(r) terra::res(r)[1], numeric(1))
  ord <- order(resos)
  levels <- levels[ord]
  resos <- resos[ord]
  to_binary <- function(r) {
    b <- r * 0 + 1
    b[is.na(b)] <- 0
    b
  }
  mats <- lapply(levels, function(r) terra::as.matrix(to_binary(r), wide = TRUE))
  exts <- lapply(levels, terra::ext)
  reso <- resos[1]

  # Search to the diagonal of the largest level
  maxdist_of <- function(r) {
    ex <- terra::ext(r)
    sqrt((ex$xmax - ex$xmin)^2 + (ex$ymax - ex$ymin)^2)
  }
  maxdist <- max(vapply(levels, maxdist_of, numeric(1)))
  if (n < 2) {
    landsea_cap <- 2 * maxdist_of(landsea)
    if (maxdist > landsea_cap) {
      warning(sprintf("n = %.3g is below 2: capping search distance at %.0f m (2 x landsea diagonal) instead of %.0f m!!", n, landsea_cap, maxdist), call. = FALSE)
      maxdist <- landsea_cap
    }
  }
  kmax <- max(8, ceiling((maxdist / reso)^(1 / n) * 8))
  s <- (c(0, (8:kmax) / 8))^n * reso
  s <- s[s <= maxdist]

  lss <- to_binary(landsea)
  lsm <- terra::as.matrix(lss, wide = TRUE)
  e <- terra::ext(landsea)
  jitter_deg <- if (jitter) 10 else 0
  rasters <- lapply(dirs, function(d) {
    lsr <- coastal_exposure_cpp(lsm, reso, e$xmin, e$ymax, s, d, mats, resos,
                                vapply(exts, function(x) x$xmin, numeric(1)), vapply(exts, function(x) x$xmax, numeric(1)),
                                vapply(exts, function(x) x$ymin, numeric(1)), vapply(exts, function(x) x$ymax, numeric(1)),
                                jitter_deg)
    .rast(lsr, lss)
  })
  if (length(dirs) > 1L) {
    out <- terra::mean(terra::rast(rasters))
    names(out) <- "coastalexposure_all"
  } else {
    out <- rasters[[1]]
    names(out) <- paste0("coastalexposure_", dirs)
  }
  out
}

#' @title Calculate coastal exposure for coastal temperature effects
#'
#' @description The function `calculate_coastalexposure` calculates the land fraction
#' upwind of each cell for a set of wind directions. The result depends only on the
#' land/sea geometry, so it can be calculated once for an area and passed via the `cex`
#' parameter to [spatialdownscale()], [spatialdownscale_tiles()], [tempdaily_downscale()]
#' and [temphrly_downscale()].
#'
#' @param dtmf a fine-resolution SpatRast of elevations (NA = sea) defining the output grid.
#' @param dtmm optionally, a SpatRast of elevations (NA = sea) covering a wider area than
#' `dtmf`. Used at its own resolution and extent for samples beyond `dtmf`. NA to use
#' `dtmf` only.
#' @param coarse optionally, further SpatRasts (NA = sea) covering wider areas, e.g. a
#' national land/sea mask (see [coastalexposure()]).
#' @param ndir number of wind directions, evenly spaced from 0 degrees. Default 32.
#' @param smooth size in cells of the moving window used to smooth the directional layers. Default 5.
#' @param n,jitter passed to [coastalexposure()].
#' @param filename optional file to which output is written (as 32-bit floats),
#' recommended for large areas.
#'
#' @details For each direction, land fraction is calculated with [coastalexposure()],
#' blended with the two neighbouring directions (weights 0.25, 0.5, 0.25) and smoothed
#' spatially with a `smooth` x `smooth` moving-window mean.
#'
#' @return a SpatRast matching `dtmf` with `ndir` + 1 layers: the smoothed land
#' fraction for each direction (named `dir_<azimuth>`), and `all`, the unsmoothed mean
#' across directions.
#' @keywords spatial
#' @export
#' @examples
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' cex<-calculate_coastalexposure(dtmf, dtmm)
#' terra::plot(cex[[c("dir_90","dir_270","all")]])
calculate_coastalexposure<-function(dtmf, dtmm=NA, coarse=NULL, ndir=32, smooth=5, n=2, jitter=TRUE, filename=""){
  if(inherits(dtmf,"PackedSpatRaster")) dtmf<-unwrap(dtmf)
  if(inherits(dtmm,"PackedSpatRaster")) dtmm<-unwrap(dtmm)
  if(inherits(coarse,c("SpatRaster","PackedSpatRaster"))) coarse<-list(coarse)
  if(inherits(dtmm,"SpatRaster")) coarse<-c(list(dtmm),coarse)
  dirs<-(0:(ndir-1))*360/ndir
  lsr<-.is(rast(lapply(dirs,function(d) coastalexposure(dtmf,d,coarse=coarse,n=n,jitter=jitter))))
  lsr2<-lsr
  for (i in 0:(ndir-1)) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%ndir+1]+0.5*lsr[,,i+1]+0.25*lsr[,,(i+1)%%ndir+1]
  lsr2<-focal(.rast(lsr2,dtmf),w=smooth,fun="mean",na.policy="omit",na.rm=TRUE)
  cex<-c(lsr2,.rast(apply(lsr,c(1,2),mean),dtmf))
  names(cex)<-c(paste0("dir_",dirs),"all")
  if(filename!="") cex<-writeRaster(cex,filename,datatype="FLT4S",overwrite=TRUE)
  return(cex)
}
#' @title Performs thin-plate spline downscaling
#' @description The function `Tpsdownscale` is a thin plate spline model, typically
#' with elevation as a covariate to downscale data.
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
#'  \dontrun{
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' rain<-climdata$prec
#' try(rainf<-Tpsdownscale(rain, climdata$dtm, dtmf, method = "normal", fast = TRUE))
#' terra::plot(rain,main='Input rain')
#' terra::plot(rainf,main='Downscaled rain')
#' }
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
    tps<-suppressWarnings(fields::fastTps(xyz,v,m=2,aRange=res(r)[1]*5))
  }  else {
    tps<-fields::Tps(xyz,v,m=2)
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
#' @title Calculates the diffuse fraction from incoming shortwave radiation
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#' @param rad a vector of incoming shortwave radiation values (W/m^2)
#' @param julian the Julian day as returned by [julday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @return a vector of diffuse fractions (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}}).
#' @export
#' @details
#' The method assumes the environment is snow free. Both overall cloud cover and heterogeneity in
#' cloud cover affect the diffuse fraction. Breaks in an extensive cloud deck may primarily
#' enhance the beam irradiance, whereas scattered clouds may enhance the diffuse irradiance and
#' leave the beam irradiance unaffected.  In consequence, if hourly data are available, an index
#' is applied to detect the presence of such variable/inhomogeneous clouds, based on variability
#' in radiation for each hour in question and values in the preceding and deciding hour.  If
#' hourly data are unavailable, an average variability is determined from radiation intensity.
#' @keywords climate
#' @examples
#' rad <- c(1:1352) # typical values of radiation in W/m^2
#' jd <- 2459752 #"2022-06-21" as julian day
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
#' @param pk - surface pressure in kPa
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
  e0 <- 0.6108 * exp((17.27 * tc)/(tc + 237.3))
  ws <- (18.02 / 28.97) * (e0 / pk)
  # ws <- (18.02 / 28.97) *pK
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
  if (intype == "relative")  hr <- h

  if (max(.is(hr), na.rm = T) > 100)
    warning(paste("Some relative humidity values > 100%",
                  max(hr, na.rm = T)))
  if (outtype == "specific") h <- (hr / 100) * ws
  if (outtype == "relative") h <- hr
  if (outtype == "absolute"){
    ea<-e0*(hr/100)
    h <- 2.16679 * (ea/tk)
  }
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
#' jd <- 2459752 #"2022-06-21" as julian day
#' jd<-2459215 # 01/01/2021
#' tme<-seq(0,23,1)
#' csr<-clearskyrad(jd,tme,60,0)
#' plot(csr ~ tme, type = "l", lwd = 2, xlab = expression(paste("Hour")), ylab = "Clearsky radiation")
#' print(sum(csr))
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
#'  \dontrun{
#' tme<-as.POSIXlt(seq(as.POSIXlt("2022-01-01"),as.POSIXlt("2022-06-30"),60*60*24*8))
#' jd<-sapply(tme,mesoclim:::.jday)
#' dl<-daylength(jd,50)
#' plot(dl ~ jd, type = "l", lwd = 8, xlab = "Day", ylab = "Day length")
#' }
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
#' @title Calculate Lapse rates via vapour pressure
#' @param tc = temperature (deg C) as vector, array or spatraster
#' @param ea = vapour pressure (kPa) as vector, array or spatraster
#' @param pk = atmospheric pressure (kPa) as vector, array or spatraster
#' @returns lapse rates in same format as tc
#' @details Alternatively a value of 0.005 can be used if lacking pk or relhum data
#' @export
#' @keywords climate spatial
#' @examples
#' climdata<-read_climdata(mesoclim::ukcpinput)
#' ea<-converthumidity(climdata$relhum,tc=climdata$temp , pk=climdata$pres)
#' lr<-lapserate(climdata$temp,ea,climdata$pres)
#' terra::plot(mesoclim:::.rast(lr,climdata$dtm)[[1]])
lapserate <- function(tc, rh=NA, pk=NA) {
  if(class(tc)[1]=="Spatraster"){
    asArrays<-FALSE
    r<-tc[[1]]
    tc<-.is(tc)
    rh<-.is(rh)
    pk<-.is(pk)
  } else asArrays<-TRUE
  ea<-.satvap(tc)*(rh/100)
  rv<-0.622*ea/(pk-ea)
  lr<-9.8076*(1+(2501000*rv)/(287*(tc+273.15)))/
    (1003.5+(0.622*2501000^2*rv)/(287*(tc+273.15)^2))
  if(!asArrays) lr<-.rast(lr,r)
  lr
}

#' Convert sea to atmospheric pressure
#'
#' @param psl - numeric or spatraster sea level pressure
#' @param dtm - numeric or spatraster elevation
#'
#' @return numeric or spatraster of sea level pressures matching dtm
#' @export
#'
#' @examples
#' print(paste("Sea level pressure of 100 kPa atmospheric pressure at 500m elevation =",round(sea_to_atmos_pressure(100,500),1),"kPa"))
sea_to_atmos_pressure<-function(psl,dtm){
  if(inherits(psl,"SpatRaster")){
    toArrays<-FALSE
    tem<-psl[[1]]
  } else toArrays<-TRUE
  psl<-.is(psl)
  dtm<-.is(dtm)
  dtm<-ifelse(is.na(dtm),0,dtm)
  if(inherits(psl,"numeric")) arraylength<-length(psl) else arraylength<-dim(psl)[3]
  pres<-psl * (((293-0.0065*.mta(dtm,arraylength))/293)^5.26)
  if(!toArrays) pres<-.rast(pres,tem)
  if(inherits(psl,"numeric")) pres<-as.vector(pres)
  return(pres)
}

#' Convert atmospheric to sea level pressure
#'
#' @param pres - vector or 3Darray or spatraster of atmospheric pressure - if array expects format[x,y,time]
#' @param dtm - vector or matrix or spatraster of elevations - if matrix or spatraster expected to match pres
#'
#' @return numeric or spatraster of atmospheric pressures matching dtm
#' @export
#'
#' @examples
#' atmos_to_sea_pressure(ukcpinput$pres,ukcpinput$dtm)
atmos_to_sea_pressure<-function(pres,dtm){
  if(inherits(pres,"SpatRaster")){
    toArrays<-FALSE
    tem<-pres[[1]]
  } else toArrays<-TRUE
  pres<-.is(pres)
  dtm<-.is(dtm)
  dtm<-ifelse(is.na(dtm),0,dtm)
  if(inherits(pres,"numeric")) arraylength<-length(pres) else arraylength<-dim(pres)[3]
  psl<-pres / (((293-0.0065*.mta(dtm,arraylength))/293)^5.26)
  if(!toArrays) psl<-.rast(psl,tem)
  if(inherits(pres,"numeric")) psl<-as.vector(psl)
  return(psl)
}

#' Calculate horizon for different solar azimuths and total skyview
#' @details Skyview places equal importance on each sector of sky
#' @param dtm digital terrain model SpatRaster
#' @param steps number of angular sectors used for skyview/horizon calculation (default 36)
#' @param toArrays logical; if TRUE returns arrays instead of SpatRasters
#' @param skyview_only logical; if TRUE returns only the skyview SpatRaster (or array),
#'   skipping construction of the per-direction horizon output. Faster when horizon
#'   angles are not needed.
#'
#' @return When \code{skyview_only = FALSE} (default): a list with elements
#'   \code{skyview} and \code{horizon} (SpatRasters or arrays depending on
#'   \code{toArrays}). When \code{skyview_only = TRUE}: the skyview SpatRaster
#'   (or array) directly.
#' @export
#'
#' @examples
#' dtmf<-terra::rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))
#' results<-calculate_terrain_shading(dtmf)
#' #plot(results$skyview)
#' #plot(results$horizon[[c(1,6,12,18)]])
calculate_terrain_shading<-function(dtm, steps=36, toArrays=FALSE, skyview_only=FALSE){
  r<-dtm
  dtm<-ifel(is.na(dtm),0,dtm)
  # Accumulate horizon angles; only allocate full hor array when needed
  sv<- array(0, dim(dtm)[1:2])
  if(!skyview_only) hor<-array(NA,dim=c(dim(dtm)[1:2],steps))
  for (i in 1:steps){
    h<-.horizon(dtm,(i-1)*(360/steps))
    sv<-sv+atan(h)
    if(!skyview_only) hor[,,i]<-h
  }
  sv<-sv/steps
  sv<-tan(sv)
  sv<-0.5*cos(2*sv)+0.5

  if(!toArrays){
    sv<-.rast(sv,dtm)
    sv<-mask(sv,r)
  }
  if(skyview_only) return(sv)

  if(!toArrays){
    hor<-.rast(hor,dtm)
    hor<-mask(hor,r)
  }
  return(list(skyview=sv,horizon=hor))
}


