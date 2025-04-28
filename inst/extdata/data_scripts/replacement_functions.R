# radsw<-daily100m$swrad
# tme<-as.POSIXlt(time(daily100m$swrad))
# clearsky = NA
# r = daily100m$dtm
 # adjust = TRUE

##### DEVELOPING FUNCTIONS

#' @title Daily to hourly downward longwave radiation
#' @description Derives an array of hourly effective sky-emissivity values.
#'
#' @param skyem - an array of daily mean sky-emissivity values values (0-1)
#' @param tme - POSIXlt object of dates corresponding to radsw
#' @param adjust  - optional logical which if TRUE ensures that, after interpolation, returned
#'          hourly values, when averaged to daily, match the input
#'
#' @return an array of hourly sky-emissivity values  (0-1)
#' @export
#' @details NB - more consistent to code this as downward longwave, but will essentially
#' do the calaculations in the function below, but with temperature as an additional
#'  Effective sky emissvity can be used to calaculate downward longwave radiation (Lwd).
#'  The formula is Lwd = skyem * Lwu where Lwu is upward longwave radiation given
#'  by Lwu=0.97*5.67*10**-8*(tc+273.15). Here tc is average surface temperature (deg C))
#'  but an adequate approximation is derived if subtited by air temperature.
#' @keywords temporal
#' @examples
#' lwhr<- lw_dailytohourly(lw=daily100m$lwrad, dtm=daily100m$dtm, hrtemps=hrtemps, hrrh=hrrh, hrpres=hrpres, tme=daily100m$tme,  adjust = FALSE)
lw_dailytohourly <- function(lw=daily100m$lwrad, dtm=daily100m$dtm, hrtemps=hrtemps, hrrh=hrrh, hrpres=hrpres, tme,  adjust = FALSE) {
  if (inherits(lw, "PackedSpatRaster")) lw<-unwrap(pres)
  if(inherits(lw, "SpatRaster")){
    toArrays<-FALSE
    tem<-lw[[1]]
  } else toArrays<-TRUE

  # Calculate LW up
  lwup<-.lwup(.is(hrtemps)) #= .lwup function 0.97*5.67*10**-8*(tc+273.15)# where tc = average surface temperature approximated by tair (hrly or daily)

  # Calculate sky emissivity
  ea<-converthumidity(.is(hrrh),intype = "relative", outtype = "vapour pressure", tc = .is(hrtemps), pk = .is(hrpres))
  tdp<-.dewpoint(ea,.is(hrtemps))
  #skyem = 5.31 x 10^-13 * T_air^2 # swindon
  skyem<- 0.787 + 0.764 * log((tdp+273.15) / 273.15) # Clear sky em from clark and allen - Tdp= dewpoint temperature

  # Adjust skyem for cloud cover N if supplied?
  # skyem <- skyem * (1 + (0.0224*N) - (0.0035*N^2) + (0.00028*N^3))

  # Calculate LW down from skyem and LW up
  lwdhr <- skyem * lwup

  ## ADJUST removes special dimension eg to temperature etc??
  # plot(.rast(lwdhr,daily100m$lwrad[[1]])[[1:12]] )
  #lwhr <- .daytohour(.is(lw))
  if (adjust) {
    lwdd <- .hourtoday(lwdhr)
    mult <- .is(lw) / lwdd
    mult[is.na(mult)] <- 1
    mult <- .daytohour(mult, Spline = FALSE)
    lwdhr <- lwdhr * mult
  }

  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  yr<-tme$year[2]
  sel<-which(tmeh$year==yr)
  lwdhr<-lwdhr[,,sel]

  if(!toArrays){
    lwdhr<-.rast(lwdhr,tem)
    time(lwdhr) <- tmeh

  }
  return(lwdhr)
}

#### NEW FUNCTIONS for SW temporal downscale - works when not already in lat lon crs
# To work for r when not already lat/lon projection
.clearskyraddaily_v2 <- function(tme, r) {
  #1440 = minutes in day
  # dmean - sums for each day
  dmean<-function(x) {
    x<-matrix(x, ncol = 1440, byrow=T)
    apply(x, 1, mean, na.rm=T)
  }
  # Get lat for each cell in r
  lats<-.latslonsfromr(r)$lats # reprojects to 4326
  #e <- ext(r)
  #lats <- seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]) # ASSUMES crs in lat lon??
  jd <- rep(juldayvCpp(tme$year + 1900, tme$mon + 1, tme$mday),  each = 1440)
  lt <- rep(c(0:1439)/60,length(tme))

  # Create matrices of times and latitudes
  n1 <- length(lats) # num of cells/cols
  n2 <- length(jd) # num of minutes
  lats <- matrix(rep(as.vector(lats), n2), ncol = n2)
  jd<-matrix(rep(jd, each = n1), ncol = n2)
  lt<-matrix(rep(lt, each = n1), ncol = n2)

  # Calculate clear sky radiation and convert to daily
  csr <- clearskyrad(lt, lats, long=0, jd)
  csd <- apply(csr, 1, dmean)
  csda <- array(rep(csd, dim(r)[2]), dim=c(length(tme), dim(r)[1:2]))
  csda <- aperm(csda, c(2,3,1))
  csda
}

# To work for r when not already lat/lon projection - see following functions as well
# outputs spatrast or array determined to match inputs
#swradhr<-swrad_dailytohourly_v2(radsw=daily100m$swrad, tme=as.POSIXlt(terra::time(unwrap(daily100m$swrad))), clearsky = NA, r = daily100m$dtm, adjust = TRUE, toArray=FALSE)

swrad_dailytohourly_v2 <- function(radsw, tme, clearsky = NA, r = r, adjust = TRUE, toArray=TRUE) {
  if (inherits(radsw, "PackedSpatRaster")) radsw<-unwrap(radsw)
  # Check geo info in either radsw or as r; convert radsw to array
  if(class(radsw)[1]=='SpatRaster') r<-radsw[[1]] else{
    if(class(r)[1]!='SpatRaster') stop('Lacking geo information. Need radsw or r to be class SpatRaster!')
  }
  # If NA, calaculate daily clearsky radiation
  if (class(clearsky) == "logical") clearsky <- .clearskyraddaily_v2(tme,r)
  # Calculate clear sky fraction
  radsw<-.is(radsw)

  radf <- radsw/clearsky
  radf[radf > 1] <- 1
  radf[radf < 0] <- 0

  # Interpolate clear sky fraction to hourly
  radfh <- .daytohour(radf)
  radfh[radfh > 1] <- 1
  radfh[radfh < 0] <- 0

  # Calculate hourly clear sky radiation
  lat<-.latslonsfromr(r)$lats # reprojects to 4326
  lon<-.latslonsfromr(r)$lons # reprojects to 4326

  tmeh <- as.POSIXlt(seq(tme[1],tme[length(tme)]+23*3600, by = 3600))
  jd<-juldayvCpp(tmeh$year+1900, tmeh$mon+1, tmeh$mday)
  lt <- tmeh$hour
  lons<-as.vector(t(lon))
  lats<-as.vector(t(lat))
  csh<-clearskyradmCpp(jd,lt,lats, lons)
  csh[is.na(csh)] <- 0
  # make spatraster
  csh.r<-rep(r,length(lt))
  values(csh.r)<-csh

  # Calculate hourly radiation
  radh <- .is(csh.r) * radfh
  # Make consistent with daily
  if (adjust) {
    radd <- .hourtoday(radh)
    mult <- radsw / radd
    mult[is.na(mult)] <- 1
    mult <- .daytohour(mult, Spline = FALSE)
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
  if(!toArray){ # convert to spatraster
    radh<-.rast(radh,r)
    time(radh) <- seq.POSIXt(tme[1],  by = "hour", length.out=length(tme)*24)
  }
  return(radh)
}



.clearskyrad <- function(tme, lat, long, tc = 15, rh = 80, pk = 101.3)

#' @title Produces a matrix of latitudes form a terra::SpatRaster object
#' Inputs:
#' @param r - a terra::SpatRaster object
#' @returns a matrix of latidues
#'@noRd
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' @title Produces a matrix of longitudes form a terra::SpatRaster object
#' Inputs:
#' @param r - a terra::SpatRaster object
#' @returns a matrix of latidues
#'@noRd
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' @title get lats and lons from raster
#' @noRd
.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}
