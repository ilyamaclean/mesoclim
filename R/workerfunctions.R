# ******** NB - suggest putting any functions we don't want to make
# available as stand-alone documented functions in here. Precede by dot as R
# studio won't then expect oxygen2 when documenting package
# ** These are functions that are definately used, checked and working

#' @title Check if input is a SpatRaster or PackedSpatRaster and convert to matrix (if single layer) or array (multiple layers)
#' @param r - possible raster
#' @keywords internal
#' @noRd
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-terra::rast(r)
  if (class(r)[1] == "SpatRaster") {
    if (dim(r)[3] == 1) {
      m<-as.matrix(r,wide=TRUE)
    } else m<-as.array(r)
  } else {
    m<-r
  }
  return(m)
}
#' @title Convert matrix or rast to array
#' @param r - matrix or raster
#' @param n - third dim
#' @details
#' Converts to array of given dimensions n - replicating 2D matrix or single layer rast n times.
#' ALWAYS converts to array (not matrix) even if n=1
#' @keywords internal
#' @noRd
.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}

#' @title Convert vector to array defined by spatraster or matrix
#' @keywords internal
#' @noRd
.vta <- function(v,r) {
  m<-.is(r)
  va<-rep(v,each=dim(m)[1]*dim(m)[2])
  a<-array(va,dim=c(dim(m),length(v)))
  a
}

#' @title Create SpatRaster object using a template
#' @import terra
#' @details expects m as 2D matrix/array or 3D array
#' @keywords internal
#' @noRd
.rast <- function(m,tem) {
  # Throw warning if dim of m do not match dim of tem
  if(any(dim(m)[1:2]!=dim(tem)[1:2])){
    warning("In .rast dimensions of matrix/array do not match dimensions of rast template!!!")
    #print(dim(m))
    #print(dim(tem))
  }
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' @title expand 3D daily array to hourly array where 3rd dim = days
#' @keywords internal
#' @noRd
.ehr<-function(a) {
  n<-dim(a)[1]*dim(a)[2]
  o1<-rep(c(1:n),24*dim(a)[3])
  o2<-rep(c(1:dim(a)[3]),each=24*n)-1
  o2<-o2*max(o1,na.rm=T)
  o<-o1+o2
  ah<-rep(a,24)
  ah<-ah[o]
  ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
  ah
}
#' @title Calculate moving average
#' @noRd
.mav <- function(x, n = 5) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}

#' @title Produces a matrix of latitudes/y coordinates form a terra::SpatRaster object
#' @param r - a terra::SpatRaster object
#' @returns a matrix of latidues
#'@noRd
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' @title Produces a matrix of longitudes/x coordinatew form a terra::SpatRaster object
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
#' @title Gets lats and lons for each cell in a spatraster
#' @details Extract x/y coordinates in native rast crs before reprojecting coordinates to EPSG 4326
#' This ensure that the dimensions macth those of the input rast - would not be case if first projected then extracted.
#' @returns named list of a 2D matrix of lats and lons
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

#' @title version of terra resample that equates to NA.RM = TRUE
#' @param r1 - is resampled to same geometry as r2
#' @param r2 - target geometry
#' @param msk=TRUE if output to be masked out where r2 cells = NA
#' @param method for resample and project can be set
#' @export
.resample <- function(r1,r2, msk=FALSE, method='bilinear'){
  if (terra::crs(r1) != terra::crs(r2)) r1<-terra::project(r1, terra::crs(r2), method)
  af<-terra::res(r2)[1] /terra::res(r1)[1]
  if (round(af,10) > 1) {			         # If resolution different aggregate
    ro<-terra::aggregate(r1, af, na.rm=TRUE)
    if (ext(ro) != terra::ext(r2)){           # if extents different then also resample
      ro<-terra::resample(ro, r2, method)
    }
  } else ro<-terra::resample(r1, r2, method)  # if res same then just resample

  if(msk) ro<-terra::mask(ro, r2)

  #ensure time and names of layers in new raster match r1
  names(ro)<-names(r1)
  terra::time(ro)<-terra::time(r1)
  return(ro)
}

#' Interpolates sea surface temperature data to coastal regions
#' Resamples sst data to same extent, resolution and projection as aoi
#' Interpolates from mean of adjacent sea cells
#' @param sst.r -  sea surface temperature data
#' @param aoi.r -  dtm of area of interest where sea = NA
#'
#' @return sst spatrasters of same projection, res and extent as aoi with matching coastline
#' @export
#'
#' @examples
#' sst.r<-rast(system.file("extdata/sst/NWSClim_NWSPPE_r001i1p00000_2018_gridT.nc",pkg="mesoclim"),"SST")[[1]]
#' aoi.r<-rast(system.file("extdata/dtms/dtmm.tif",pkg="mesoclim"))
#' plot(.sea_to_coast(sst.r,aoi.r))
.sea_to_coast<-function(sst.r,aoi.r, ext_cells=8){
  # Extend area to make sure sea surface temperature data matches coastal areas
  aoibuf.r<-terra::extend(aoi.r,c(ext_cells,ext_cells),fill=1)
  # Re-project/sample sst to extended aoi
  newsst.r<-project(sst.r,aoibuf.r)
  # Interpolate coastal sea cells
  target<-sum(c(crop(newsst.r[[1]],aoi.r),aoi.r),na.rm=T)
  n<-1
  while(anyNA(values(target)) & n<=5){
    plot(target,main=n)
    newsst.r<-focal(newsst.r, w=9, fun=mean, na.policy="only", na.rm=T)

    target<-sum(c(crop(newsst.r[[1]],aoi.r),aoi.r),na.rm=T)
    n<-n+1
    if(n>5) warning("CHECK interpolation of coastal cells in create_ukcpsst_data function!!!")
  }
  newsst.r<-mask(crop(newsst.r,aoi.r),aoi.r,inverse=TRUE)
  return(newsst.r)
}


# ============================================================================ #
# Replacements of SSTinterpolate
# spatial_interA is a direct copy of SSTinterpolate code
# time_interp  uses terra::approximate but valid for any combi of timeseries
# ============================================================================ #
#' @title Spatially interpolate missing values in Spatrasters
#' @param r 2 or 3D Spatraster
#' @return Spatraster of same dim as r but without NA
#' @export
#' @import terra
#' @keywords internal
#' @examples
#' r <- terra::rast(ncols=5, nrows=5, vals=rep(c(1,2,NA,4,5),5))
#' rout<-mesoclim:::.spatinterp(r)
#' terra::plot(c(r,rout))
.spatinterp<-function(r){
  tme<-terra::time(r)
  me<-as.vector(r)
  n<-which(is.na(me))
  rout<-r
  crs(rout)<-crs(r) # prevents superfluous warnings
  if (length(n) > 0) {
    dmx<-max(dim(r)[1:2])
    m1<-.is(r)
    m2<-.is(resample(aggregate(r,2,na.rm=T),r))
    s<-which(is.na(m1))
    m1[s]<-m2[s]
    if (dmx >= 5) {
      m2<-.is(resample(aggregate(r,5,na.rm=T),r))
      m1[s]<-m2[s]
    }
    if (dmx >= 10) {
      m2<-.is(resample(aggregate(r,10,na.rm=T),r))
      s<-which(is.na(m1))
      m1[s]<-m2[s]
    }
    if (dmx >= 100) {
      m2<-.is(resample(aggregate(r,100,na.rm=T),r))
      s<-which(is.na(m1))
      m1[s]<-m2[s]
    }
    s<-which(is.na(m1))
    if (length(s) > 0) {
      me<-unlist(global(r,mean,na.rm=TRUE))
      m2<-.vta(me,r[[1]])
      m1[s]<-m2[s]
    }
    rout<-.rast(m1,r)
  }
  terra::time(rout)<-tme
  return(rout)
}
#' @title Temporal interpolation of Spatraster timeseries
#' @param rin = 3D Spatraster
#' @param tmein = NA or Date or POSIX vector of datetimes of r - if NA assumes time(r) provides tmein values
#' @param tmeout = Date or POSIX vector of datetimes required
#' @return Interpolated (linear) Spatraster stack of nlyrs=length(tmeout)
#' @export
#' @import terra
#' @import lubridate
#' @keywords internal
#' @examples
#' r <- terra::rast(ncols=5, nrows=5, vals=rep(1,25))
#' r<-c(r,r*2,r*3,r*4,r*5)
#' tmein<-as.Date(seq(17532,17652,30))
#' terra::time(r)<-tmein
#' tmeout<-as.Date(seq(17532,17651,1))
#' rout<-mesoclim:::.tmeinterp(r,NA,tmeout)
#' terra::plot(rout[[c(1,15,31,95,119,120)]])
.tmeinterp<-function(r,tmein,tmeout){
  # Check formats of inputs
  if(class(r)[1]== "PackedSpatRaster") r<-rast(r)
  if(class(r)[1]!= "SpatRaster") stop("Input r must be a SpatRaster!!!")
  if(class(tmein)[1]=='Date') tmein<-as.POSIXlt(tmein)
  if(class(tmeout)[1]=='Date') tmeout<-as.POSIXlt(tmeout)
  if(class(tmein)[1]=='logical') tmein<-as.POSIXlt(terra::time(r))
  if(class(tmein)[1]!="POSIXlt" | class(tmeout)[1]!="POSIXlt") stop("Time parameters must be POSIXlt!!!")
  if(nlyr(r)!=length(tmein)) stop('Number of layers in r does not match length of input times!!!')

  # Check tmeout within tmein
  if(tail(tmeout,1)<tmein[1] | tmeout[1]>tail(tmein,1)) stop('Output times do not overlap with input times!!!')
  if(tmeout[1]<tmein[1]|tail(tmeout,1)>tail(tmein,1)) warning('Output times extend beyond input times!!!')

  # If only single tmein - simply replicate across tmeout
  if(length(tmein)==1){
    warning('Single input timeseries - returning constant output!!!')
    rout<-rep(r,length(tmeout))
  }

  if(length(tmein)>1){
    # Check relative size of timesteps in tmein and tmeout
    instep<-as.numeric(tmein[2])-as.numeric(tmein[1])
    outstep<-as.numeric(tmeout[2])-as.numeric(tmeout[1])
    if(instep/outstep>744) warning('Input time step MUCH lower than output time step - this may be slow!!!')

    # Extend tmeout to include tmein at or after last tmeout and at or before tmeout[1]
    if(!tmeout[1] %in% tmein) tstart<-tmein[tail(which(tmein<tmeout[1]),n=1)] else tstart<-tmeout[1]
    if(!tail(tmeout,1) %in% tmein) tend<-tmein[which(tmein>tail(tmeout,1))[1]] else tend<-tail(tmeout,1)
    tme<-as.POSIXlt(seq(tstart,tend,outstep))

    # Create output rast
    rout<-rast(setValues(r[[1]],NA),nlyr=length(tme))
    terra::time(rout)<-tme
    names(rout)<-tme

    # Insert and fill missing date layers of output r - could use zoo::na.approx or na.spline
    rout[[c(which(tme %in% tmein))]]<-r[[which(tmein %in% tme)]]
    rout<-terra::approximate(rout,method="linear")

    # Return only values in tmeout
    sel<-which(as.POSIXlt(terra::time(rout)) %in% tmeout)

  }
  return(rout[[sel]])
}

# ============================================================================ #
# ~~~~~~~~~ Climate processing worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #

#' @title Calculate saturated vapour pressure
#' @noRd
.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
}
#' @title Calculates the astronomical Julian day
#'@noRd
.jday <- function(tme) {
  yr<-tme$year+1900
  mth<-tme$mon+1
  dd<-tme$mday+(tme$hour+(tme$min+tme$sec/60)/60)/24
  madj<-mth+(mth<3)*12
  yadj<-yr+(mth<3)*-1
  jd<-trunc(365.25*(yadj+4716))+trunc(30.6001*(madj+1))+dd-1524.5
  b<-(2-trunc(yadj/100)+trunc(trunc(yadj/100)/4))
  jd<-jd+(jd>2299160)*b
  jd
}
#' @title Calculates solar time
#'@noRd
.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-localtime+(4*(long-merid)+eot)/60-dst
  st
}
#' @title Calculates the solar altitude
#'@noRd
.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  sa<-(180*atan(sh/sqrt(1-sh^2)))/pi
  sa
}
#' @title Calculates solar azimuth
#'@noRd
.solazi <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  hh<-0.5*pi-acos(sh)
  sazi<-cos(d)*sin(tt)/cos(hh)
  cazi<-(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))/
    sqrt((cos(d)*sin(tt))^2+(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))^2)
  sqt <- 1-sazi^2
  sqt[sqt<0]<-0
  solz<-180+(180*atan(sazi/sqrt(sqt)))/pi
  solz[cazi<0 & sazi<0]<-180-solz[cazi<0 & sazi<0]
  solz[cazi<0 & sazi>=0]<-540-solz[cazi<0 & sazi>=0]
  solz
}
#' @title Calculate skyview
#'
#' @param dtm digital terrain spatRaster - usually of at downscale resolution and extent
#' @param steps - number of segments to use in calculation
#' @return spatRaster of skyview suitable for use as a parameter to [`swdownscale`] and [`lwdownscale`] functions.
#' @noRd
.skyview<-function(dtm,steps=36) {
  r<-dtm
  dtm[is.na(dtm)]<-0
  ha <- array(0, dim(dtm)[1:2])
  for (s in 1:steps) { # uses horizon angle in calc but places equal importance on each sector of sky
    ha<-ha+atan(.horizon(dtm,s*360/steps))
  }
  ha<-ha/steps
  ha<-tan(ha)
  svf<-0.5*cos(2*ha)+0.5
  svf<-.rast(svf,dtm)
  svf<-mask(svf,r)
  return(svf)
}

#' @title Simulate cloud or rain patchiness
#' @import gstat
#' @import terra
#' @noRd
.simpatch<-function(varf,af,n=dim(varf)[3],varn="precip") {
  if (varn == "precip") {
    varf<-log(varf+1)
    sdr<- -0.14934+0.1777*log(af)
    if (sdr> 1) sdr<-1
  }
  if (varn == "csfrac") {
    vf<-.is(varf)
    s0<-which(vf==0)
    s1<-which(vf==1)
    vf<-suppressWarnings(log(vf/(1-vf)))
    vf[is.infinite(vf)]<-NA
    vf[s0]<-min(vf,na.rm=TRUE)
    vf[s1]<-max(vf,na.rm=TRUE)
    varf<-.rast(vf,varf)
    sdr<- 0.012950+0.175718*log(af)
    if (sdr> 1) sdr<-1
  }
  # Generate n simulations
  xy <- expand.grid(1:dim(varf)[2], 1:dim(varf)[1])
  names(xy) <- c('x','y')
  g1<-gstat(formula=z~1, locations=~x+y, dummy=T, beta=0,
            model=vgm(psill=1,
                      range=af/2,
                      nugget=0,
                      model='Sph'), nmax = 40)
  ni<-dim(varf)[3]
  hps<-round(ni/(n-1),0)
  if (n != ni) {
    days<-ni/hps
  } else days <- ni+1
  yy1 <- predict(g1,newdata=xy,nsim=days-1,debug.level=0)
  anom<-.rast(yy1,varf[[1]])
  anom<-resample(aggregate(anom,2),anom)
  anom<-mask(anom,varf[[1]])
  ni<-dim(varf)[3]
  # if n not equal to ni generate ni anomoly layers by averaging
  # allows for temporal autocorrelation in simulations
  if (n != ni) {
    ad<-.is(anom)
    wgts<-c(0:(hps-1))/hps
    am<-array(NA,dim=dim(varf))
    ij<-1
    for (i in 1:days) {
      for (j in 1:hps) {
        am[,,ij]<-(1-wgts[j])*ad[,,i]+wgts[j]*ad[,,i+1]
        ij<-ij+1
      }
    }
    anom<-.rast(am,anom[[1]])
  }
  # adjust to ensure anomaly centered on zero and with correct sd
  me<-apply(.is(anom),3,mean,na.rm=T)
  sds<-apply(.is(anom),3,sd,na.rm=T)
  sda<-apply(.is(varf),3,sd,na.rm=T)
  mu <- (sda*sdr)/sds
  a<-.is(anom)
  a<-(a-.vta(me,varf[[1]]))*.vta(mu,varf[[1]])
  anom<-.rast(a,anom[[1]])
  varnf<-varf+anom
  if (varn == "precip") {
    varnf<-exp(varnf)-1
    varnf[varnf<0]<-0
  }
  if (varn == "csfrac") {
    varnf<-1/(1+exp(-varnf))
  }
  return(varnf)
}
#' @title Appply elevation adjustment to cloud fractional cover
#' @import terra
#' @noRd
.cfcelev<-function(csfc,dtmf,dtmc) {
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  if (crs(csfc) != crs(dtmf)) csfc<-project(csfc,crs(dtmf))
  # Calculate mean
  ca<-.is(csfc)
  cmean<-apply(ca,c(1,2),mean,na.rm=TRUE)
  # logit transform
  s0<-which(cmean==0)
  s1<-which(cmean==1)
  lc<-log(cmean/(1-cmean))
  lc[is.infinite(lc)]<-NA
  lc[s0]<-min(lc,na.rm=TRUE)
  lc[s1]<-max(lc,na.rm=TRUE)
  # resample
  lc<-.resample(.rast(lc,dtmc),dtmf)
  # Apply elevation adjustment to mean
  dc<-.is(.resample(dtmc,dtmf))
  ddif<-.is(dtmf)-dc
  ldif<- -9.938e-03-6.550e-04*ddif+4.966e-05*dc+1.358e-06*dc*ddif
  lcn<-.rast(.is(lc)+ldif,dtmf)
  lcn<-mask(lcn,dtmf)
  # Calculate difference in total between elevation adjusted and resampled / resampled
  mu<-.is(lcn-lc)
  mu<-.rta(mu,dim(csfc)[3])
  # Calculated individual logit transoformed radiation
  cf<-.resample(csfc,dtmf)
  cf<-.is(cf)
  s0<-which(cf==0)
  s1<-which(cf==1)
  cf<-log((cf)/(1-cf))
  cf[is.infinite(cf)]<-NA
  cf[s0]<-min(cf,na.rm=TRUE)
  cf[s1]<-max(cf,na.rm=TRUE)
  cf<-cf+mu
  # Back transform
  cf<-1/(1+exp(-cf))
  cf<-.rast(cf,dtmf)
  cf<-mask(cf,dtmf)
  return(cf)
}
#' @title Calculates wind altitude coefficient in specified direction
#' @noRd
.windz<-function(dtm1,dtm2,dtmr,wdir) {
  reso1<-res(dtm1)[1]
  reso2<-res(dtm2)[2]
  # Calculate shifts
  # Generate sequence of numbers distances over which to sample
  d<-exp(c(0:150)/10)
  s<-which(d >= reso1 & d <= reso2)
  d2<-d[s]
  xshift<-round((sin(wdir*pi/180)*d2)/reso1)
  yshift<-round((cos(wdir*pi/180)*d2)/reso1)
  # Calculate unique shifts
  df1<-data.frame(x=xshift,y=yshift)
  df1<-unique(df1)
  xshift<-df1$x
  yshift<-df1$y
  # Calculate buffer
  bdist<-max(abs(yshift),abs(xshift))+1
  # Create matrix with buffer
  m2<-array(NA,dim=c(dim(dtm1)[1]+bdist*2,dim(dtm1)[2]+bdist*2))
  m2[(bdist+1):(bdist+dim(dtm1)[1]),(bdist+1):(bdist+dim(dtm1)[2])]<-.is(dtm1)
  wc<-array(1,dim=dim(dtm1)[1:2])
  m<-.is(dtm1)
  m3<-.is(dtmr) # replacement raster
  for (i in 1:length(xshift)) {
    # elevation difference
    ed<-m-m2[(bdist+1-yshift[i]):(bdist-yshift[i]+dim(m)[1]),(bdist+1+xshift[i]):(bdist+xshift[i]+dim(m)[2])]
    mu<-suppressWarnings(log(67.8*(ed+2)-5.42)/4.8699)
    mu[ed<0]<-1
    mu[mu<1]<-1
    s<-which(is.na(mu))
    mu[s]<-m3[s]
    wc<-wc+mu
  }
  wc<-wc/length(xshift)
  wc<-.rast(wc,dtm1)
  return(wc)
}
#' @title Calculates wind shelter coefficient in specified direction
#' @noRd
.windcoef <- function(dtm, direction, hgt = 1) {
  reso<-res(dtm)[1]
  m<-.is(dtm)
  m[is.na(m)]<-0
  m<-m/reso
  hgt<-hgt/reso
  azi<-direction*(pi/180)
  horizon <- array(0, dim(m))
  m3 <- array(0, dim(m) + 200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  m3[101:(x + 100), 101:(y + 100)] <- m
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(m3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-m3[101:(x+100),101:(y+100)])/(step^2),na.rm=T)
    horizon[1:x,1:y]<-ifelse(horizon[1:x,1:y]<(hgt/step^2),0,horizon[1:x,1:y])
  }
  index<-1-atan(0.17*100*horizon)/1.65
  r<-.rast(index,dtm)
  return(r)
}

#' @title Applies height correction to wind speed measurements
#' @description cf microclima::windheight
#' wspeed numeric value(s) of measured wind speed (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{ m s^{-1}}}) at height `zi` (m).
#' @param zi a numeric value idicating the height (m) above the ground of `wspeed` input
#' @param zo a numeric value indicating the height (m) above ground level of output speeds
#' Assumes a logarithmic height profile and imposes a minimum zo of 0.2.
#' Equivalent of using a roughness length (z0) of ~0.2 where v = vref ln(z/z0)/ln(zref/z0) (https://www.rensmart.com/Information/WindSheer)
#' @keywords internal
.windhgt<-function (wspeed, zi, zo) {
  if (zo < 0.2 & zo > (5.42/67.8)){
    warning("Wind-height profile function performs poorly below 20 cm so output height converted to 20 cm")
    zo <- 0.2
  }
  return(wspeed * log(67.8 * zo - 5.42)/log(67.8 * zi - 5.42))
}


#' @title Calculates horizon angle
#' @noRd
.horizon <- function(dtm, azimuth) {
  reso<-res(dtm)[1]
  dtm<-.is(dtm)
  dtm[is.na(dtm)]<-0
  dtm<-dtm/reso
  azi<-azimuth*pi/180
  horizon<-array(0,dim(dtm))
  dtm3<-array(0,dim(dtm)+200)
  x<-dim(dtm)[1]
  y<-dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)]<-dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2),na.rm=T)
  }
  horizon
}

# ============================================================================ #
# ~~~~~~~~~ Basin delineation worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Wrapper for C++ function to delineate basins
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @noRd
.basindelinCpp<-function(dtm) {
  dm<-.is(dtm)
  dm[is.na(dm)]<-9999
  dm2<-array(9999,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # (2) create blank basin file
  bsn<-dm2*NA
  dun<-array(0,dim=dim(bsn))
  bsn<-basinCpp(dm2, bsn, dun)
  dd<-dim(bsn)
  bsn<-bsn[2:(dd[1]-1),2:(dd[2]-1)]
  if(class(bsn)[1]!='matrix') bsn<-matrix(bsn,ncol=ncol(dtm),nrow=nrow(dtm))
  r<-.rast(bsn,dtm)
  return(r)
}
#' @title Internal function for delineating basins with option for boundary > 0
#' @noRd
.basindelin<-function(dtm, boundary = 0) {
  # Delineate basins
  dm<-dim(dtm)
  me<-mean(as.vector(dtm),na.rm=TRUE)
  if (is.na(me) == FALSE) {
    bsn<-.basindelinCpp(dtm)
    # Merge basins if boundary > 0
    if (boundary > 0) {
      mx<-max(as.vector(bsn),na.rm=T)
      tst<-1
      while (tst == 1) {
        u<-unique(as.vector(bsn))
        u<-u[is.na(u) == F]
        if (length(u) > 1) {
          bsn<-.basinmerge(dtm,bsn,boundary)
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
          if (length(u) > 1) {
            bsn<-.basinmerge(dtm,bsn,boundary)
            u<-unique(as.vector(bsn))
            u<-u[is.na(u) == F]
          }
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
        }
        if (length(u) == 1) tst<-0
        mx2<-max(as.vector(bsn),na.rm=T)
        if (mx2 ==  mx) {
          tst<-0
        } else mx<-mx2
      } # end while
    } # end if boundary
  } else bsn<-dtm # end if boundary
  return(bsn)
}
#' @title function to identify which basin edge cells are less or equal to boundary
#' @noRd
.edge<-function(v) {
  o<-0
  if (is.na(v[1]) == FALSE) {
    if (max(v,na.rm=TRUE) > v[1]) o<-1
  }
  o
}
#' @title function to assign which surrounding cells should be merged
#' @noRd
.edgec<-function(v) {
  o<-v*0
  if (is.na(v[1]) == FALSE) {
    s<-which(v>v[1])
    o[s]<-1
  }
  o
}
#' @title function to grab neighbouring cells and reassign basin number
#' @noRd
.asign3<-function(bm2,bea,rw,cl) {
  b3<-bm2[rw:(rw+2),cl:(cl+2)]
  v<-bea[rw,cl,]
  if (is.na(v[2])==FALSE & v[2] > 0) b3[2,1]<-b3[2,2]
  if (is.na(v[3])==FALSE & v[3] > 0)  b3[2,3]<-b3[2,2]
  if (is.na(v[4])==FALSE & v[4] > 0)  b3[1,2]<-b3[2,2]
  if (is.na(v[5])==FALSE & v[5] > 0)  b3[1,1]<-b3[2,2]
  if (is.na(v[6])==FALSE & v[6] > 0)  b3[1,3]<-b3[2,2]
  if (is.na(v[7])==FALSE & v[7] > 0) b3[3,2]<-b3[2,2]
  if (is.na(v[8])==FALSE & v[8] > 0)  b3[3,1]<-b3[2,2]
  if (is.na(v[9])==FALSE & v[9] > 0)  b3[3,3]<-b3[2,2]
  b3
}
#' @title Merge basins based on specified boundary
#' @noRd
.basinmerge<-function(dtm,bsn,boundary=0.25) {
  # Put buffer around basin and dtn
  bm<-.is(bsn)
  bm2<-array(NA,dim=c(dim(bm)[1]+2,dim(bm)[2]+2))
  bm2[2:(dim(bm)[1]+1),2:(dim(bm)[2]+1)]<-bm
  dm<-.is(dtm)
  dm2<-array(NA,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # Create 3D array of  basin numbers  with adjoining cells
  bma<-array(NA,dim=c(dim(bm),9))
  bma[,,1]<-bm # rw, cl
  bma[,,2]<-bm2[2:(dim(bm)[1]+1),1:dim(bm)[2]] # rw, cl-1
  bma[,,3]<-bm2[2:(dim(bm)[1]+1),3:(dim(bm)[2]+2)] # rw, cl+1
  bma[,,4]<-bm2[1:dim(bm)[1],2:(dim(bm)[2]+1)] # rw-1, cl
  bma[,,5]<-bm2[1:dim(bm)[1],1:dim(bm)[2]] # rw-1, cl-1
  bma[,,6]<-bm2[1:dim(bm)[1],3:(dim(bm)[2]+2)] # rw-1, cl+1
  bma[,,7]<-bm2[3:(dim(bm)[1]+2),2:(dim(bm)[2]+1)] # rw+1, cl
  bma[,,8]<-bm2[3:(dim(bm)[1]+2),1:dim(bm)[2]] # rw+1, cl-1
  bma[,,9]<-bm2[3:(dim(bm)[1]+2),3:(dim(bm)[2]+2)] # rw+1, cl+1
  # Create 3D array of elevation differences with adjoining cells
  dma<-array(NA,dim=c(dim(dm),9))
  dma[,,1]<-dm # rw, cl
  dma[,,2]<-dm2[2:(dim(dm)[1]+1),1:dim(dm)[2]]-dm # rw, cl-1
  dma[,,3]<-dm2[2:(dim(dm)[1]+1),3:(dim(dm)[2]+2)]-dm  # rw, cl+1
  dma[,,4]<-dm2[1:dim(dm)[1],2:(dim(dm)[2]+1)]-dm  # rw-1, cl
  dma[,,5]<-dm2[1:dim(dm)[1],1:dim(dm)[2]]-dm  # rw-1, cl-1
  dma[,,6]<-dm2[1:dim(dm)[1],3:(dim(dm)[2]+2)]-dm  # rw-1, cl+1
  dma[,,7]<-dm2[3:(dim(dm)[1]+2),2:(dim(dm)[2]+1)]-dm  # rw+1, cl
  dma[,,8]<-dm2[3:(dim(dm)[1]+2),1:dim(dm)[2]]-dm  # rw+1, cl-1
  dma[,,9]<-dm2[3:(dim(dm)[1]+2),3:(dim(dm)[2]+2)]-dm  # rw+1, cl+1
  dma2<-dma*0
  dma2[abs(dma)<boundary]<-1
  bma<-bma*dma2
  bma[,,1]<-bm
  # identify edge and basin merge cells
  be<-apply(bma,c(1,2),.edge)
  bea<-aperm(apply(bma,c(1,2),.edgec),c(2,3,1))
  s<-which(be>0,arr.ind=TRUE)
  for (i in 1:dim(s)[1]) {
    rw<-as.numeric(s[i,1])
    cl<-as.numeric(s[i,2])
    b3<-.asign3(bm2,bea,rw,cl)
    bm2[rw:(rw+2),cl:(cl+2)]<-b3
  }
  # reassign basin number
  u<-unique(as.vector(bm2))
  u<-u[is.na(u)==FALSE]
  u<-u[order(u)]
  bm3<-bm2
  for (i in 1:length(u)) {
    s<-which(bm2==u[i])
    bm3[s]<-i
  }
  dd<-dim(bm3)
  bsn<-bm3[2:(dd[1]-1),2:(dd[2]-1)]
  # case where single column/row and bsn is a vector
  if(class(bsn)[1]!='matrix') bsn<-matrix(bsn,ncol=ncol(dtm),nrow=nrow(dtm))
  r<-.rast(bsn,dtm)
}

#' @title Mosaic tiled basins merging common joins
#' @noRd
.basinmosaic<-function(b1,b2) {
  e1<-ext(b1)
  e2<-ext(b2)
  reso<-res(b1)
  # *********** Do this if the tiles are vertically adjoined  *************** #
  if (abs(e1$ymax-e2$ymax) > reso[1]) {
    if (e2$ymax > e1$ymax) {  # b2 above b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b1 above b2
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on top row of b1
      v1<-m1[1,] # top row of b1
      n<-dim(m2)[1] # mumber of rows
      v2<-m2[n,] # bottom row of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$ymax > e1$ymax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  } else {# end do this if the tiles are vertically adjoined
    # *********** Do this if the tiles are horizontally adjoined  ************** #
    if (e2$xmax > e1$xmax) {  # b2 right of b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b2 left of b1
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on right hand column of b1
      n<-dim(m1)[2]
      v1<-m1[,n] # right-hand column of b1
      v2<-m2[,1] # left-hand column of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$xmax > e1$xmax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  }
  # ********************************** Mosaic ******************************* #
  bout<-mosaic(b1n,b2n)
  return(bout)
}
#' @title Do an entire column of tiled basins
#' @noRd
.docolumn<-function(dtm,tilesize,boundary,x) {
  e<-ext(dtm)
  reso<-res(dtm)
  ymxs<-as.numeric(ceiling((e$ymax-e$ymin)/reso[2]/tilesize))-1
  xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
  xmx<-xmn+reso[1]*tilesize
  ymn<-as.numeric(e$ymin)+reso[2]*tilesize*0
  ymx<-ymn+reso[2]*tilesize
  if (xmx > e$xmax) xmx<-e$xmax
  if (ymx > e$ymax) ymx<-e$ymax
  ec<-ext(xmn,xmx,ymn,ymx)
  dc<-crop(dtm,ec)
  bma<-basindelin(dc,boundary)
  # delineate basins for columns
  for (y in 1:ymxs) {
    xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
    xmx<-xmn+reso[1]*tilesize
    ymn<-as.numeric(e$ymin)+reso[2]*tilesize*y
    ymx<-ymn+reso[2]*tilesize
    if (xmx > e$xmax) xmx<-e$xmax
    if (ymx > e$ymax) ymx<-e$ymax
    ec<-ext(xmn,xmx,ymn,ymx)
    dc<-crop(dtm,ec)
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-basindelin(dc,boundary)+ta
    bma<-.basinmosaic(bma,bo)
  } # end y
  return(bma)
}
#' @title Function used for delineating basins with big dtms
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @noRd
.basindelin_big<-function(dtm, boundary = 0, tilesize = 100, plotprogress = FALSE) {
  # chop into tiles
  e<-ext(dtm)
  reso<-res(dtm)
  xmxs<-as.numeric(ceiling((e$xmax-e$xmin)/reso[1]/tilesize))-1
  bma<-.docolumn(dtm,tilesize,boundary,0)
  for (x in 1:xmxs) {
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-.docolumn(dtm,tilesize,boundary,x)+ta
    ed<-Sys.time()
    bma<-.basinmosaic(bma,bo)
    if (plotprogress) plot(bma,main=x)
  }
  m<-.is(bma)
  # renumber basins
  u<-unique(as.vector(m))
  u<-u[is.na(u) == FALSE]
  u<-u[order(u)]
  m<-array(renumberbasin(m,u),dim=dim(m))
  bout<-.rast(m,dtm)
  return(bout)
}
#' @title Calculate flow direction
# '@noRd
.flowdir <- function(md) {
  fd<-md*0
  md2<-array(NA,dim=c(dim(md)[1]+2,dim(md)[2]+2))
  md2[2:(dim(md)[1]+1),2:(dim(md)[2]+1)]<-md
  v<-c(1:length(md))
  v<-v[is.na(md) == F]
  x<-arrayInd(v,dim(md))[,1]
  y<-arrayInd(v,dim(md))[,2]
  for (i in 1:length(x)) {
    md9<-md2[x[i]:(x[i]+2),y[i]:(y[i]+2)]
    fd[x[i],y[i]]<-round(mean(which(md9==min(md9,na.rm=TRUE))),0)
  }
  fd
}
# ============================================================================ #
# ~~~~~~~~~~~~~~ Spatial downscale worker functions here ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #

# ============================================================================ #
# ~~~~~~~~~~~~ Temperature downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #

# ~~~~~~~~~~~~ Temperature altitude downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#' @title Downscales temperature for elevation effects
#' @keywords internal
.tempelev <- function(tc, lrf, lrc, dtmf, dtmc = NA) {
  if (class(dtmc)[1] == "logical")  dtmc<-.resample(dtmf,tc)
  if (class(tc)[1] == "array") tc<-.rast(tc,dtmc)

  # Convert NA to elevation of 0 in dtmc
  dtmc<-ifel(is.na(dtmc),0,dtmc)

  # Lapse rate multipliers x elev
  #lrc<-as.array(lrc)*.rta(dtmc,n)
  #lrf<-as.array(lrf)*.rta(dtmf,n)
  lrc<-lrc*dtmc
  lrf<-lrf*dtmf

  # Sea-level temperature from coarse resolution lapse rate
  #stc<-.rast(.is(tc)+lrc,dtmc)
  stc<-tc+lrc
  if (crs(dtmc) != crs(dtmf)) stc<-project(stc,crs(dtmf))
  stc<-.resample(stc,dtmf)

  # Actual temperature from resampled sea temps and fine resolution lapse rate
  #tcf<-suppressWarnings(stc-.rast(lrf,dtmf))
  tcf<-suppressWarnings(stc-lrf)
  return(tcf)
}

# ~~~~~~~~~~~~ Temperature Cold Air Drainage downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#' Calculate cold air drainage lapse rate multiplier
#'
#' @param dtmf - fine scale dtm spatraster
#' @param basins - fine scale map of basins - if NA will calculate from the dtmf
#' @param refhgt - Height of temperature measruements above ground
#'
#' @return spatRaster of lapse rate multiplier mu
#' @details calculates the lapse rate multiplier by the difference
#'  in elevation between each point and the highest in the
#' basin.
#' @keywords internal
#'
#' @examples
.cad_multiplier<-function(dtmf, basins = NA, refhgt = 2){
  # Calculate elevation difference between basin height point and pixel
  if (class(basins) == "logical") basins<-basindelin(dtmf,refhgt)
  b<-.is(basins)
  d<-.is(dtmf)
  u<-unique(as.vector(b))
  u<-u[is.na(u)==FALSE]
  bmx<-b*0
  for (i in 1:length(u)) {
    s<-which(b==u[i])
    mx<-max(d[s],na.rm=TRUE)
    bmx[s]<-mx
  }
  edif<-bmx-.is(dtmf)
  edif<-.rast(edif,dtmf)

  # Calculate lapse-rate multiplication factor (basin elev difference)
  mu<-edif*.cadpotential(dtmf,basins,refhgt)
  return(mu)
}
#' @title Calculate cold-air drainage potential (spatial)
#' @noRd
.cadpotential <- function(dtm, basins = NA, refhgt = 2) {
  if (class(basins) == "logical") basins<-basindelin(dtm,refhgt)
  fa<-flowacc(dtm, basins) - 1
  # Calculate basin size
  fre<-freq(basins)
  b<-.is(basins)
  bsize<-b*0
  for (i in 1:length(fre$value))  {
    s<-which(b==fre$value[i])
    bsize[s]<-fre$count[i]
  }
  cadfr<-.is(fa)/bsize
  cadfr<-.rast(cadfr,dtm)
  cadfr[cadfr>1]<-1
  return(cadfr)
}
#' Calculate whether conditions allow for cold air drainage
#'
#' @param tc - temperature in deg celsius
#' @param u2 - windspeed in m/s
#' @param swad - SW downward radiation
#' @param lwrad - LW downward radiation
#' @param dtmc - coarse resolution dtm matching climate variable inputs - not required if other parameters all spatRasters.
#' @param dtmf - fine scale dtm for outputs
#' @param refhgt - of temp above ground (m)
#'
#' @return a binary spatRaster matching dtmf where 1= suitable, 0 = not suitable for cold air drainage
#' @details temperature, wind  and radiation data used to determine
#' whether cold air drainage conditions are likely to occur and returns a binary
#' vector of the same length as `em` and `wind` indicating whether conditions
#' occur (1) or not (0). For daily measurements it may be appropriate to apply only to daily min
#' and  to set swrad to 0 assuming this minimum occurs during the night.
#' @keywords internal
#'
#' @examples
.cad_conditions<-function(tc,u2,swrad,lwrad,dtmc,dtmf,refhgt = 2){
  # determine whether cold-air drainage conditions exist
  d<-0.65*0.12
  zm<-0.1*0.12
  uf<-(0.4*u2)/log((refhgt-d)/zm)
  H<-(swrad+lwrad-(5.67*10^-8*0.97*(tc+273.15)^4))*0.5
  st<- -(0.4*9.81*(refhgt-d)*H)/(1241*(tc+273.15)*uf^3)
  if(!inherits(st,"SpatRaster")) st<-.rast(st,dtmc)
  if (crs(st) != crs(dtmf)) st<-project(st,crs(dtmf))
  st<-.resample(st,dtmf)
  st<-ifel(st>1,1,st)
  st<-ifel(st<1,0,st)
  return(st)
}

#' Calculate cold air drainage correction
#'
#' @param lrf spatraster of lapse rates
#' @param mu spatraster of lapse rate multiplier output by`.cad_multiplier`.
#' @param st spatraster of suitability for cold air drainage as output by `.cad_conditions`
#'
#' @return spatraster of cold air drainage correction in deg C
#' @keywords internal
#'
#' @examples
.apply_cad<-function(lrf,mu,st){
  cad<-lrf*-mu
  ce<-cad*st
  return(ce)
}

#' @title DEFUNCT Downscale temperature with cold air drainage effects
#' @export
#' @keywords internal
.tempcad<-function(climdata, tempvar, lrf, dtmf, basins = NA, refhgt = 2) {
  # Calculate elevation difference between basin height point and pixel
  if (class(basins) == "logical") basins<-basindelin(dtmf,refhgt)
  b<-.is(basins)
  d<-.is(dtmf)
  u<-unique(as.vector(b))
  u<-u[is.na(u)==FALSE]
  bmx<-b*0
  for (i in 1:length(u)) {
    s<-which(b==u[i])
    mx<-max(d[s],na.rm=TRUE)
    bmx[s]<-mx
  }
  edif<-bmx-.is(dtmf)
  edif<-.rast(edif,dtmf)

  # Calculate lapse-rate multiplication factor (basin elev difference)
  mu<-edif*.cadpotential(dtmf,basins,refhgt)

  # extract climate variables
  relhum<-climdata$relhum
  pk<-climdata$pres
  tc<-climdata[[tempvar]]
  dtmc<-rast(climdata$dtm)
  # Calculate lapse rate
  ea<-.satvap(tc)*(relhum/100)
  lr<-lapserate(tc, ea, pk)
  lr<-.rast(lr,dtmc)
  if (crs(lr) != crs(dtmf)) lr<-project(lr,crs(dtmf))
  lrf<-.resample(lr,dtmf)
  n<-dim(lrf)[3]
  cad<-.is(lrf)*-.rta(mu,n)

  # determine whether cold-air drainage conditions exist
  d<-0.65*0.12
  zm<-0.1*0.12
  # Extract additional climate variables
  u2<-climdata$windspeed
  swrad<-climdata$swrad
  lwrad<-climdata$lwrad
  uf<-(0.4*u2)/log((refhgt-d)/zm)
  H<-(swrad+lwrad-(5.67*10^-8*0.97*(tc+273.15)^4))*0.5
  st<- -(0.4*9.81*(refhgt-d)*H)/(1241*(tc+273.15)*uf^3)
  st<-.rast(st,dtmc)
  if (crs(st) != crs(dtmf)) st<-project(st,crs(dtmf))
  st<-.resample(st,dtmf)
  st<-.is(st)
  st[st>1]<-1
  st[st<1]<-0
  ce<-.rast(cad*st,dtmf)
  return(ce)
}

# ~~~~~~~~~~~~ Temperature Coastal Effect downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#' @title  applying correction factor to coastal effects to account for data resolution
#' @noRd
.correctcoastal<-function(r) {
  reso<-res(r)[1]
  am<- -6.74389+1.08141*log(reso)
  ap<- -3.1163+0.4818*log(reso)
  bm<- -0.6167+0.3197*log(reso)
  bp<- 1.5842-0.7922*log(reso)+0.1104*log(reso)^2
  m<-.is(r)
  l<-log(m/(1-m))
  sm<-which(m> 0 & m <= 0.5)
  sp<-which(m< 1 & m > 0.5)
  l[sm]<-am+bm*l[sm]
  l[sp]<-ap+bp*l[sp]
  pm<-1/(1+exp(-l))
  s0<-which(m==0)
  s1<-which(m==1)
  pm[s0]<-0
  pm[s1]<-1
  rp<-.rast(pm,r)
  return(rp)
}
#' @title Coastal temperature effects
#'
#' @param tc - downscaled temperature at same res as dtmf
#' @param sstf - downscaled sea surface temperature to dtmf resolution and extent - no NA and timeseries must match tc
#' @param u2 - downscaled windspeed at temprature height
#' @param wdir - wind direction (coarse resolution) - same value for all of dtmf extent will be used
#' @param dtmf - fine dtm spatraster
#' @param dtmm - fine-scale dtm covering wider area than dtmf (but same resolution!!)
#' @param dtmc - coarse dtm spatraster matching resolution and extent of wdir
#' @param ndir - number of directions to calculate for wind exposure
#' @param smooth - number of cells to use for smoothing matrix
#' @param correct - corrects each timestep  so that output mean temp of area matches input area mean temp
#' @return Spatraster of temperature that includes coastal effect
#' @export
#' @keywords internal
.tempcoastal<-function(tc, sstf, u2, wdir, dtmf, dtmm, dtmc,ndir=32,smooth=5,correct=FALSE) {
  # Resample dtmm to dtmf resolution
  if(any(res(dtmm)!=res(dtmf))) dtmm<-.resample(dtmm,dtmf,msk=TRUE)

  # Calculate coastal exposure for each wind direction - NOT very realistic for complex coasts if ndir=8
  landsea<-ifel(is.na(dtmm),NA,1)
  lsr<-array(NA,dim=c(dim(dtmf)[1:2],ndir))
  for (i in 0:(ndir-1)) {
    r<-coastalexposure(landsea, e=ext(dtmf), i%%ndir*(360/ndir)) # seems to return a raster 1 row and col too big!
    r<-.correctcoastal(r)
    lsr[,,i+1]<-.is(r)
  }
  # smooth
  lsr2<-lsr
  for (i in 0:(ndir-1)) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%ndir+1]+0.5*lsr[,,i%%ndir+1]+0.25*lsr[,,(i+1)%%ndir+1]
  lsr2<-.is( focal(.rast(lsr2,dtmf),w=smooth,fun="mean",na.policy="omit",na.rm=TRUE) )

  lsm<-apply(lsr,c(1,2),mean)
  tst<-min(lsm,na.rm=T)
  if (tst < 1) { # only apply coastal effects if there are coastal area
    # slot in wind speeds
    if(!inherits(wdir,"SpatRaster")) wdr<-.rast(wdir,dtmc) else wdr<-wdir
    ew<-ext(wdr)
    xy<-data.frame(x=(ew$xmin+ew$xmax)/2,y=(ew$ymin+ew$ymax)/2)
    wdir<-as.numeric(extract(wdr,xy))[-1]
    if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,median,na.rm=TRUE)
    # Calculate array land-sea ratios for every hour
    i<-round(wdir/(360/ndir))%%ndir
    lsr<-lsr2[,,i+1]
    # Calculate sstf weighting upwind
    # derive power scaling coefficient from wind speed
    p2<-0.10420*sqrt(.is(u2))-0.47852
    # calculate logit lsr and lsm
    llsr<-log(lsr/(1-lsr))
    llsm<-.rta(log(lsm/(1-lsm)),dim(lsr)[3])
    # Calculate mins and maxes
    s<-which(is.na(lsr) == FALSE & lsr > 0 & lsr < 1)
    llsr[lsr==0]<-log(min(lsr[s])/(1-min(lsr[s])))
    llsr[lsr==1]<-log(max(lsr[s])/(1-max(lsr[s])))
    s<-which(is.na(lsm) == FALSE & lsm > 0 & lsm < 1)
    llsm[lsm==0]<-log(min(lsm[s])/(1-min(lsm[s])))
    llsm[lsm==1]<-log(max(lsm[s])/(1-max(lsm[s])))
    # predict logit swgt
    lswgt<- -0.1095761+p2*(llsr+3.401197)-0.1553487*llsm
    swgt<-.rast(1/(1+exp(-lswgt)),tc)
    tcp<-swgt*sstf+(1-swgt)*tc
    # Correct area mean temps to equal area mean input temps???? - NOT suitable for tiles
    # Increases temperature range significantly!!
    if(correct){
      af<-res(dtmc)[1]/res(dtmf)[1]
      tcc<-resample(aggregate(tcp,af,na.rm=TRUE),tcp)
      tcp<-tc+(tcp-tcc)
    }
  } else tcp<-tc
  return(tcp)
}

# ** Following is a bit of a code dump. We won't need it all:
# NB:
#  ** (1) For several of these functions we'll need to add the appropriate imports
#     (e.g. raster terra) before they work (see example for .cropnc).
#  ** (2) it would be good to get rid of raster as a dependency as I suspect it
#         will be deprecated owing to rgdal dependency, but for now brick seems
#         to be faster than terra equivalent
#  ** (3) I've preceded all by a dot. But I wrote them before this was done, so they
#          they might still call functions without the dot.
# ============================================================================= #


# ============================================================================ #
# ~~~~~~~~~~~~~ Data handing worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #

#' @title Quick crop of one variable nc file
#' @param fi - filename of nc file
#' @param e - raster extent object of the extent to crop
#' @return a SpatRaster object
#' @import raster
#' @import terra
#' @noRd
.cropnc<-function(fi,e) {
  r<-suppressWarnings(brick(fi))
  r<-crop(r,e)
  r<-rast(r)
  r
}

#' @title extracts array data from nc file via spatRaster if cropping required
#' @param filein - file name of nc file
#' @param varid - name of variable in nc file - required
#' @param aoi - spatraster or extent to which nc data will be cropped
#' @details
#' Data values of new beta CDS ncdf format is not reliably read by terra,
#' therefore gets spatial info via rast but values via ncdf4::ncvar_get
#' @import terra
#' @import ncdf4
#' @noRd
.nctoarray <- function(filein, varid, aoi = NA) {
   nc <- nc_open(filein)
   if (class(aoi)[1] %in% c('SpatRaster','SpatExtent')){
    r<-rast(filein,varid)
    r<-setValues(r,aperm(ncvar_get(nc, varid = varid), c(2,1,3))) # new format era5 ncdf vals cannot be read by reliably terra
    a <-as.array(crop(r,aoi))
  } else a <- aperm(ncvar_get(nc, varid = varid), c(2,1,3))
  nc_close(nc)
  return(a)
}

#' @title extracts array data from nc file TO DELETE!!!???
#' @param filein - file name of nc file
#' @param varid - name of variable in nc file
#' @noRd
.nctoarray_old <- function(filein, varid = NA, aoi = NA) {

  a <- aperm(ncvar_get(nc, varid = varid), c(2,1,3))
  #a <- apply(a, c(2,3), rev)
  nc_close(nc)
  a
}
#' @title Converts matrix to array by replicating all entries for quick multiplication
#' @param m - a matrix
#' @param n - dim(3) of returned array
#' @return returns a 3D array
#' @noRd
.mta <- function(m,n) {
  a<-array(rep(m,n),dim=c(dim(m),n))
  a
}
#' @title Apply equivelent for arrays with NAs. Outperforms apply by an order of magnitude
#' @param a - expects a 3D array (typically with NAs for sea) but will convert 1D array or vector assuming a timeseries
#' @param fun - a function to apply
#' @returns a 3D array of e.g. daily values form hourly
#' @noRd
.applynotna<-function(a,fun) {
  if(inherits(a,'numeric')|length(dim(a))==1) a<-array(a,dim=c(1,1,length(a)))
  m<-matrix(a,ncol=dim(a)[1]*dim(a)[2],byrow=T)
  sel<-which(is.na(m[1,])==F)
  r<-apply(m[,sel,drop=FALSE],2,fun)
  n<-dim(r)[1]
  ao<-array(NA,dim=c(dim(a)[1:2],n))
  sel<-which(is.na(a[,,1:n])==F)
  ao[sel]<-aperm(r,c(2,1))
  ao
}
#' @title Applies e.g. a landsea mask to an array of data
#' @param a - a 3D array
#' @param mask - a raster or spatrast object the extent of which matches a
#' @returns a 3D array with masked areas set to NA
#' @noRd
.applymask <- function(a, mask) {
  m<-.is(mask)
  m[is.na(m)==F]<-1
  m<-.mta(m,dim(a)[3])
  a<-a*m
  a
}
#' @title Converts an array, vector or SpatRaster of hourly data to daily
#' @param a - a 3D array, SpatRaster or vector
#' @param fun - a function , typically mean, min, max or sum
#' @returns a 3D array, SpatRaster or vector of daily data matching class
#' @noRd
.hourtoday<-function(a,fun=mean) {
    if(inherits(a,"SpatRaster")){
      tem<-a[[1]]
      a<-.is(a)
      toSpatRaster<-TRUE
    } else toSpatRaster<-FALSE
    if(inherits(a,'numeric')|length(dim(a))==1){
      a<-array(a,dim=c(1,1,length(a)))
      toVector<-TRUE
    } else toVector<-FALSE

    .htd<-function(x) {
      y<-matrix(x,ncol=24,byrow=T)
      apply(y,1,fun,na.rm=T)
    }
    d<-.applynotna(a,.htd)
    if (toVector) d<-as.vector(d)
    if (toSpatRaster) d<-.rast(d,tem)
    return(d)
  }
#' @title Converts an array daily data to hourly using one of two methods.
#' @param name description a - an array or vector of daily data - assumes timeseries if only one dimension
#' @param Spline - optional logical indicating which method to use (see details)
#' @return Returns an array of hourly data such that dim(ah)[3] == 24 * dim(h)[3]
#' @details:
#'  If Sprine = TRUE data are spline interpolated using zoo::na.approx.
#'  If Spline = FALSE each hour is given the same value as the daily data
#' @noRd
.daytohour<-function(a, Spline = TRUE, toVector=FALSE) {
  if(inherits(a,'numeric')|length(dim(a))==1){
    a<-array(a,dim=c(1,1,length(a)))
    toVector<-TRUE
  }
  if (Spline) {
    sel<-c(1:dim(a)[3])*24-12
    ah<-array(NA,dim=c(dim(a)[1:2],dim(a)[3]*24))
    for (i in 1:dim(a)[1]) {
      for (j in 1:dim(a)[2]) {
        tst<-mean(a[i,j,],na.rm=T)
        if (is.na(tst) == F) {
          ah[i,j,sel]<-a[i,j,]
          ah[i,j,1]<-a[i,j,1]
          ah[i,j,dim(ah)[3]]<-a[i,j,dim(a)[3]]
          ah[i,j,]<-zoo::na.approx(ah[i,j,])
        }
      }
    }
  } else {
    n<-dim(a)[1]*dim(a)[2]
    o1<-rep(c(1:n),24*dim(a)[3])
    o2<-rep(c(1:dim(a)[3]),each=24*n)-1
    o2<-o2*max(o1)
    o<-o1+o2
    ah<-rep(a,24)
    ah<-ah[o]
    ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
    ah
  }
  if (toVector) ah<-as.vector(ah)
  return(ah)
}
#' @title Applies coastal correction to e.g. era5 diurnal temperature ranges
#' @param dtr - diurnal temperature range (deg C)
#' @param landsea - fraction of land relative to sea
#' @returns diurnal temperature range with coastal correction (deg C)
#' @noRd
.coastalcorrect<-function(dtr, landsea) {
  lp<-.is(landsea)
  m<-(1-lp)*1.285+1
  m<-.mta(m,dim(dtr)[3])
  dtr<-dtr*m
  return(dtr)
}
#' @title writes out an array as a terra SpatRaster object
#' @param a - a 3D array
#' @param  description template - a SpatRaster object with the extent and coordinate reference system matching a
#' @param path - directory in which to save data
#' @param varn - name of variable to be saved
#' @noRd
.saverast<-function(a,template,path,varn) {
  r<-rast(a)
  ext(r)<-ext(template)
  crs(r)<-crs(template)
  fo<-file.path(path,paste0(varn,".tif"))
  writeRaster(r,filename=fo,overwrite=TRUE)
}
#' @title Get latitude and longitude of centre of r
#' @returns dataframe of 1 row and 2 columns (lat,long)
#' @noRd
.latlongfromrast<-function (r) {
  e <- ext(r)
  xy <- data.frame(x = (e$xmin + e$xmax)/2, y = (e$ymin + e$ymax)/2)
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[2], long = sf::st_coordinates(ll)[1])
  return(ll)
}
# ============================================================================ #
# ~~~~~~~~~ Climate processing worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Calculate dewpoint
#' @noRd
.dewpoint <- function(ea, tc) {
  e0<-611.2/1000
  L<-(2.501*10^6)-(2340*tc)
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tdew<-1/it-273.15
  e0<-610.78/1000
  L<-2.834*10^6
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tfrost<-1/it-273.15
  sel<-which(Tdew<0)
  Tdew[sel]<-Tfrost[sel]
  Tdew
}
#' @title Calculate clear sky radiation
#' @noRd
.clearskyrad <- function(tme, lat, long, tc = 15, rh = 80, pk = 101.3) {
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
#' @title Calculates average daily clear sky radiation over all pixels of a SpatRaster object and all days in tme
#' @param tme - a POSIXlt object of dates
#' @param r - a terra::SpatRaster object
#' @return a 3D array of expected daily clear sky radiation values (W/m**2)
#' @noRd
#' # To work for r when not already lat/lon projection
.clearskyraddaily <- function(tme, r) {
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
# ============================================================================ #
# ~~~~~~~~~~~~~ Temporal downscale worker functions here ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Writes a ncd4 file to disk
#' @param r - a terra::SpatRaster object of data to be written out
#' @param year = the year for which data are to be written out (used for time stamp)
#' @param varn - climate variable name (e.g. tasmax)
#' @param varl - climate variable long name (Daily maximum temperature at 2m)
#' @param unit - units of varn (e.g. deg C)
#' @param asinteger - optional lgoical indicating whether values should be stored as integers (dee details)
#' @param mult - multiplier to use if data are stored as integers (dee details)
#' @details:
#' Setting asintiger to TRUE roughly halves the size of the file written out but
#' reduces precision. E.g If mult = 100, a temperature value of 9.87694 would be
#' written out as 988. The function automatically corrects the units written to
#' the nc file - e.g. if mult is 100 and unit is 'deg C', 'deg C * 100' is written
#' @return out as the unit in the nc file
#'@noRd
.writenc<-function(r,fo,year,varn,varl,unit,timestep=c('daily','hourly'),asinteger= TRUE,mult=100) {
  if (asinteger) {
    rte<-r
    a<-as.array(r)
    a<-array(as.integer(a*mult),dim=dim(a))
    r<-rast(a)
    ext(r)<-ext(rte)
    crs(r)<-crs(rte)
    if (mult != 1) {
      varl<-paste(varl,"*",mult)
      unit<-paste(unit,"*",mult)
    }
  }
  st<-as.POSIXlt(0,origin=paste0(year,"-01-01 00:00"),tz="UTC")
  fn<-as.POSIXlt(0,origin=paste0(year,"-12-31 23:00"),tz="UTC")
  if(timestep=='hourly') tme<-as.POSIXlt(seq(st,fn,by=3600)) else tme<-as.POSIXlt(seq(st,fn,by=3600*24))
  terra::time(r)<-tme
  writeCDF(r,fo,overwrite=TRUE,compression=9,varname=varn,
           longname=varl,unit=unit)
}

#' @title terra version of microclima get_dem
#' @import elevatr
#' @export
#' @noRd
get_dem<-function(r, zeroasna = TRUE) {
  if (!curl::has_internet()) {
    message("Please connect to the internet and try again.")
    return(NULL)
  }
  reso<-min(res(r))
  if (reso < 30) {
    warning("Higher resolution data only available for some locations. DEM
             may be derived by interpolation")
  }
  ll<-.latlongfromrast(r)
  z<-ceiling(log((cos(ll$lat*pi/180)*2*pi*6378137)/(256*reso),2))
  prj<-sf::st_crs(r)
  dtm<-elevatr::get_aws_terrain(r,z,prj)
  dtm<-resample(dtm,r)
  if (zeroasna) dtm[dtm<=0]<-NA
  return(dtm)
}
