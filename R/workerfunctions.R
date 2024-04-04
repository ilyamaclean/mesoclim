# ******** NB - suggest putting any functions we don't want to make
# available as stand-alone documented functions in here. Precede by dot as R
# studio won't then expect oxygen2 when documenting package
# ** These are functions that are definately used, checked and working
#' Check if input is a SpatRaster or PackedSpatRaster and convert to matrix or array
#' if it is
#' @import terra
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] == "SpatRaster") {
    if (dim(r)[3] == 1) {
      m<-as.matrix(r,wide=TRUE)
    } else m<-as.array(r)
  } else {
    m<-r
  }
  return(m)
}
#' Convert matrix or rast to array
.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}
#' Convert vector to array
.vta <- function(v,r) {
  m<-.is(r)
  va<-rep(v,each=dim(m)[1]*dim(m)[2])
  a<-array(va,dim=c(dim(m),length(v)))
  a
}
#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' expand daily array to hourly array
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
#' Calculate moving average
.mav <- function(x, n = 5) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}
#' Produces a matrix of latitudes form a terra::SpatRaster object
#' Inputs:
#' r - a terra::SpatRaster object
#' Returns a matrix of latidues
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Produces a matrix of longitudes form a terra::SpatRaster object
#' Inputs:
#' r - a terra::SpatRaster object
#' Returns a matrix of longitudes
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' @import terra
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
#' version of terra resample that equates to NA.RM = TRUE
#' r1 resampled to same geometry as r2
#' msk=TRUE if output to be masked out where r2 cells = NA
#' method for resample and project can be set
#' @import terra
.resample <- function(r1,r2, msk=FALSE, method='bilinear'){
  if (crs(r1) != crs(r2)) r1<-project(r1, crs(r2), method)
  af<-res(r2)[1] /res(r1)[1]
  if (round(af,10) > 1) {			         # If resolution different aggregate
    ro<-aggregate(r1, af, na.rm=TRUE)
    if (ext(ro) != ext(r2)){           # if extents different then also resample
      ro<-resample(ro, r2, method)
    }
  } else ro<-resample(r1, r2, method)  # if res same then just resample

  if(msk) ro<-mask(ro, r2)

  #ensure time and names of layers in new raster match r1
  names(ro)<-names(r1)
  terra::time(ro)<-terra::time(r1)
  return(ro)
}
# ============================================================================ #
# ~~~~~~~~~ Climate processing worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' Calculate saturated vapour pressure
.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
}
#' Calculates the astronomical Julian day
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
#' Calculates solar time
.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-localtime+(4*(long-merid)+eot)/60-dst
  st
}
#' Calculates the solar altitude
.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  sa<-(180*atan(sh/sqrt(1-sh^2)))/pi
  sa
}
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


#' Simulate cloud or rain patchiness
#' @import gstat
#' @import terra
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
  } else days <- ni-1
  yy1 <- predict(g1,newdata=xy,nsim=days+1,debug.level=0)
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
#' Appply elevation adjustment to cloud fractional cover
#' @import terra
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
#' Calculates wind altitude coefficient in specified direction
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
#' Calculates wind shelter coefficient in specified direction
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
# Calculates horizon angle
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
.skyview<-function(dtm,steps=36) {
  r<-dtm
  dtm[is.na(dtm)]<-0
  ha <- array(0, dim(dtm)[1:2])
  for (s in 1:steps) {
    ha<-ha+atan(.horizon(dtm,s*360/steps))
  }
  ha<-ha/steps
  ha<-tan(ha)
  svf<-0.5*cos(2*ha)+0.5
  svf<-.rast(svf,dtm)
  svf<-mask(svf,r)
  return(svf)
}
# ============================================================================ #
# ~~~~~~~~~ Basin delineation worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' Wrapper for C++ function
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
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
  r<-.rast(bsn,dtm)
  return(r)
}
#' Internal function for delineating basins with option for boundary > 0
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
# function to identify which basin edge cells are less or equal to boundary
.edge<-function(v) {
  o<-0
  if (is.na(v[1]) == FALSE) {
    if (max(v,na.rm=TRUE) > v[1]) o<-1
  }
  o
}
# function to assign which surrounding cells should be merged
.edgec<-function(v) {
  o<-v*0
  if (is.na(v[1]) == FALSE) {
    s<-which(v>v[1])
    o[s]<-1
  }
  o
}
# function to grab neighbouring cells and reassign basin number
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
#' Merge basins based on specified boundary
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
  r<-.rast(bsn,dtm)
}

#' Mosaic tiled basins merging common joins
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
#' Do an entire column of tiled basins
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
#' Function used for delineating basins with big dtms
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
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
  bout<-.rast(m,bout)
  return(bma)
}
#' Calculate flow direction
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
#' Calculate cold-air drainage potential (spatial)
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
#' function for applying correction factor to coastal effects to account for
#' data resolution
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
  af<-1000/reso
  if (af > 1) {
    rpa<-.resample(aggregate(rp,af,na.rm=TRUE),rp)
    rpa<-mask(rpa,rp)
  } else rpa<-rp
  return(rpa)
}

# Component temperature effect spatial downscale functions
# ============================================================================ #
# ~~~~~~~~~~~~ Temperature downscale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' Downscales temperature for elevation effects
.tempelev <- function(tc, dtmf, dtmc = NA, rh = NA, pk = NA) {
  if (class(dtmc)[1] == "logical")  dtmc<-.resample(dtmf,tc)
  if (class(tc)[1] == "array") tc<-.rast(tc,dtmc)

  # Convert NA to elevation of 0 in dtm's
  dtmc<-ifel(is.na(dtmc),0,dtmc)

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
      lrf<-.resample(lrcp,dtmf)
    } else lrf<-.resample(lrc,dtmf)
    lrc<-as.array(lrc)*.rta(dtmc,n)
    lrf<-as.array(lrf)*.rta(dtmf,n)
  }
  # Sea-level temperature
  stc<-.rast(.is(tc)+lrc,dtmc)
  if (crs(dtmc) != crs(dtmf)) stc<-project(stc,crs(dtmf))
  stc<-.resample(stc,dtmf)
  # Actual temperature
  tcf<-suppressWarnings(stc-.rast(lrf,dtmf))
  return(tcf)
}
#' Downscale temperature with cold air drainage effects
.tempcad<-function(climdata, dtmf, basins = NA, refhgt = 2) {
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
  tc<-climdata$climarray$temp
  dtmc<-rast(climdata$dtmc)
  # Calculate lapse rate
  ea<-.satvap(tc)*(relhum/100)
  lr<-lapserate(tc, ea, pk)
  lr<-.rast(lr,dtmc)
  if (crs(lr) != crs(dtmf)) lr<-project(lr,crs(dtmf))
  lr<-.resample(lr,dtmf)
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
  st<-.resample(st,dtmf)
  st<-.is(st)
  st[st>1]<-1
  st[st<1]<-0
  ce<-.rast(cad*st,dtmf)
  return(ce)
}
.tempcoastal<-function(tc, SST, u2, wdir, dtmf, dtmm, dtmc) {
  # produce land sea mask
  if (crs(dtmm) != crs(dtmf)) dtmm<-project(dtmm,crs(dtmf))
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  if (crs(SST) != crs(dtmf)) SST<-project(SST,crs(dtmf))
  landsea<-.resample(dtmm,dtmf)
  landsea<-mask(landsea,dtmf)
  lsr<-array(NA,dim=c(dim(dtmf)[1:2],8))
  for (i in 0:7) {
    r<-coastalexposure(landsea, ext(dtmf), i%%8*45)
    r<-.correctcoastal(r)
    lsr[,,i+1]<-.is(r)
  }
  # smooth
  lsr2<-lsr
  for (i in 0:7) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%8+1]+0.5*lsr[,,i%%8+1]+0.25*lsr[,,(i+1)%%8+1]
  lsm<-apply(lsr,c(1,2),mean)
  tst<-min(lsm,na.rm=T)
  if (tst < 1) { # only apply coastal effects if there are coastal area
    # slot in wind speeds
    wdr<-.rast(wdir,dtmc)
    ll<-.latlongfromrast(wdr)
    xy<-data.frame(x=ll$long,y=ll$lat)
    wdir<-as.numeric(xx<-extract(wdr,xy))[-1]
    if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,median,na.rm=TRUE)
    # Calculate array land-sea ratios for every hour
    i<-round(wdir/45)%%8
    lsr<-lsr2[,,i+1]
    # Calculate SST weighting upwind
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
    tcp<-swgt*SST+(1-swgt)*tc
    # calculate aggregation factor
    af<-res(dtmc)[1]/res(dtmf)[1]
    tcc<-resample(aggregate(tcp,af,na.rm=TRUE),tcp)
    tcp<-tc+(tcp-tcc)
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

# Quick crop of one variable nc file
# fi - filename of nc file
# e - raster extent object of the extent to crop to
# returns a SpatRaster object
#' @import raster
#' @import terra
.cropnc<-function(fi,e) {
  r<-suppressWarnings(brick(fi))
  r<-crop(r,e)
  r<-rast(r)
  r
}

# extracts array data from nc file
# filein - file name of nc file
# varid - name of variable in nc file
.nctoarray <- function(filein, varid = NA) {
  nc <- nc_open(filein)
  a <- aperm(ncvar_get(nc, varid = varid), c(2,1,3))
  #a <- apply(a, c(2,3), rev)
  nc_close(nc)
  a
}
# Converts matrix to array by replicating all entries for quick multiplication
# m - a matrix
# n - dim(3) of returned array
# returns a 3D array
.mta <- function(m,n) {
  a<-array(rep(m,n),dim=c(dim(m),n))
  a
}
# Apply equivelent for arrays with NAs. Outperforms apply by an order of magnitude
# a - a 3D array (typically with NAs for sea)
# fun - a function to apply
# returns a 3D array of e.g. daily values form hourly
.applynotna<-function(a,fun) {
  m<-matrix(a,ncol=dim(a)[1]*dim(a)[2],byrow=T)
  sel<-which(is.na(m[1,])==F)
  r<-apply(m[,sel],2,fun)
  n<-dim(r)[1]
  ao<-array(NA,dim=c(dim(a)[1:2],n))
  sel<-which(is.na(a[,,1:n])==F)
  ao[sel]<-aperm(r,c(2,1))
  ao
}
# Applies e.g. a landsea mask to an array of data
# a - a 3D array
# mask - a raster object the extent of which matches a
# returns a 3D array with masked areas set to NA
.applymask <- function(a, mask) {
  m<-getValues(mask,format="matrix")
  m[is.na(m)==F]<-1
  m<-.mta(m,dim(a)[3])
  a<-a*m
  a
}
# Converts an array of hourly data to daily
# a - a 3D array
# fun - a function , typically mean, min, max or sum
# returns a 3D array of daily data - e.g. daily mean, max or min temperature or total rainfall
.hourtoday<-function(a,fun=mean) {
  .htd<-function(x) {
    y<-matrix(x,ncol=24,byrow=T)
    apply(y,1,fun,na.rm=T)
  }
  d<-.applynotna(a,.htd)
  d
}
# Converts an array daily data to hourly using one of two methods.
# Inputs:
# a - an array of daily data
# Spline - optional logical indicating which method to use (see details)
# Returns an array of hourly data such that dim(ah)[3] == 24 * dim(h)[3]
# Details:
#  If Sprine = TRUE data are spline interpolated using zoo::na.approx.
#  If Spline = FALSE each hour is given the same value as the daily data
.daytohour<-function(a, Spline = TRUE) {
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
          ah[i,j,]<-na.approx(ah[i,j,])
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
  return(ah)
}
# Applies coastal correction to e.g. era5 diurnal temperature ranges
# dtr - diurnal temperature range (deg C)
# landsea - fraction of land relative to sea
# returns diurnal temperature range with coastal correction (deg C)
.coastalcorrect<-function(dtr, landsea) {
  lp<-getValues(landsea,format="matrix")
  m<-(1-lp)*1.285+1
  m<-.mta(m,dim(dtr)[3])
  dtr<-dtr*m
  return(dtr)
}
# writes out an array as a terra SpatRaster object
# a - a 3D array
# template - a SpatRaster object with the extent and coordinate reference system matching a
# path - directory in which to save data
# varn - name of variable to be saved
.saverast<-function(a,template,path,varn) {
  r<-rast(a)
  ext(r)<-ext(template)
  crs(r)<-crs(template)
  fo<-paste0(path,varn,".tif")
  writeRaster(r,filename=fo,overwrite=TRUE)
}
#' Get latitude and longitude of centre of r
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
#' Calculate dewpoint
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
#' Calculate clear sky radiation
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
# Calculates average daily clear sky radiation over all pixels of a SpatRaster object and all days in tme
# Inputs:
# tme - a POSIXlt object of dates
# r - a terra::SpatRaster object
# Returns a 3D array of expected daily clear sky radiation values (W/m**2)
.clearskyraddaily <- function(tme, r) {
  dmean<-function(x) {
    x<-matrix(x, ncol = 1440, byrow=T)
    apply(x, 1, mean, na.rm=T)
  }
  e <- ext(r)
  lats <- seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1])
  jd <- rep(julday(tme$year + 1900, tme$mon + 1, tme$mday),  each = 1440)
  lt <- rep(c(0:1439)/60,length(tme))
  # Create matrices of times and latitudes
  n1 <- length(lats)
  n2 <- length(jd)
  lats <- matrix(rep(lats, n2), ncol = n2)
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
# ~~~~~~~~~~~~~ Bias-correct worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# Converts era5 hourly data to daily to enable bias correction to be applied to ukcp data
# filein - filename (including path )of era5 nc file with all variables
# pathout - directory in which to save data
# landsea - a landsea raster object of land fractions that must match the extent of the era5 data
# The function saves individual terra SpatRaster objects to disk (one for each variable)
.era5todaily<-function(filein,pathout,landsea) {
  landsea[landsea==0]<-NA
  te<-rast(landsea)
  suppressWarnings(dir.create(pathout))
  # dtr #
  a<-nctoarray(filein,"t2m")
  a<-applymask(a,landsea)
  tmx<-hourtoday(a,max)
  tmn<-hourtoday(a,min)
  dtr<-tmx-tmn
  # apply coastal correction
  dtr<-coastalcorrect(dtr,landsea)
  tmean<-(tmx+tmn)/2
  tmx<-(tmean+0.5*dtr)-273.15
  tmn<-(tmean-0.5*dtr)-273.15
  # Save tmx and tmn
  .saverast(tmx,te,pathout,"tasmax")
  .saverast(tmn,te,pathout,"tasmin")
  # psl
  a2<-nctoarray(filein,"sp")/1000
  a2<-applymask(a2,landsea)
  psl<-hourtoday(a2)
  .saverast(psl,te,pathout,"psl")
  # huss #
  a3<-nctoarray(filein,"d2m")
  a3<-applymask(a3,landsea)
  ea<-satvap(a3-273.15)
  es<-satvap(a-273.15)
  rh<-(ea/es)*100
  a3<-suppressWarnings(converthumidity(rh, intype = "relative", outtype = "specific",
                                       tc = a-273.15, pk = a2))
  huss<-hourtoday(a3)
  .saverast(huss,te,pathout,"huss")
  # u10
  a2<-nctoarray(filein,"u10")
  a2<-applymask(a2,landsea)
  uas<-hourtoday(a2)
  .saverast(uas,te,pathout,"uas")
  #v10
  a2<-nctoarray(filein,"v10")
  a2<-applymask(a2,landsea)
  vas<-hourtoday(a2)
  .saverast(vas,te,pathout,"vas")
  # rss
  a2<-nctoarray(filein,"ssrd")
  a2<-applymask(a2,landsea)
  rss<-hourtoday(a2)/3600
  .saverast(rss,te,pathout,"rss")
  # skyem
  lwdn<-nctoarray(filein,"msdwlwrf")
  lwme<-nctoarray(filein,"msnlwrf")
  lwup<-(-lwme+lwdn)
  skyem<-lwdn/lwup
  skyem<-applymask(skyem,landsea)
  skyem[skyem>1]<-1
  skyem<-hourtoday(skyem)
  .saverast(skyem,te,pathout,"skyem")
  # pr
  a<-nctoarray(filein,"tp")
  a[is.na(a)]<-0
  pr<-hourtoday(a,sum)*1000
  pr<-applymask(pr,landsea)
  .saverast(pr,te,pathout,"pr")
}
# Takes global ukcp data as inputs, ensures correct number of dates in each year,
# interpolating missing values, resamples data to match era5 and returns data for a whole true decade
# ukcpfile1 - filename of first file coinciding wiht decade - e.g. for decade 2010-2020,
# data with an extension 20091201-20191130.nc
# ukcpfile2 - filename of second file coinciding wiht decade - e.g. for decade 2010-2020,
# data with an extension 20191201-20291130.nc
# landsea - a landsea raster object of land fractions that must match the extent of the era5 data.
# must also match the resolution of era5 data as it is used for resampling ukcp data
# decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
# returns a SpatRaster object for an entire decade
.cropandsortUKCPone <- function(ukcpfile1,ukcpfile2,landsea,decade=1) {
  # Get and crop nc files
  ecrop<-extent(landsea)
  r1<-cropnc(ukcpfile1,ecrop)
  r2<-cropnc(ukcpfile2,ecrop)
  # Find out which entries belong to true decade
  st<-as.POSIXlt(0,origin=paste0(2000+(decade-1)*10+9,"-12-01 00:00"),tz="UTC")
  fn<-as.POSIXlt(0,origin=paste0(2000+decade*10+9,"-11-30 00:00"),tz="UTC")
  dt1<-round(as.numeric(seq(st,fn,length.out=3600))/(24*3600),0)
  dt1<-as.POSIXlt(dt1*24*3600,origin="1970-01-01 00:00",tz="UTC")
  st<-as.POSIXlt(0,origin=paste0(2000+decade*10+9,"-12-01 00:00"),tz="UTC")
  fn<-as.POSIXlt(0,origin=paste0(2000+(decade+1)*10+9,"-11-30 00:00"),tz="UTC")
  dt2<-round(as.numeric(seq(st,fn,length.out=3600))/(24*3600),0)
  dt2<-as.POSIXlt(dt2*24*3600,origin="1970-01-01 00:00",tz="UTC")
  styear<-2000+decade*10
  edyear<-styear+9
  sel1<-which(dt1$year+1900 >=styear & dt1$year+1900<=edyear)
  sel2<-which(dt2$year+1900 >=styear & dt2$year+1900<=edyear)
  a1<-as.array(r1)[,,sel1]
  a2<-as.array(r2)[,,sel2]
  a<-abind(a1,a2)
  # ** Create data.frame of 360 dates for real decade
  df1<-data.frame(v1="f360",dates=c(dt1[sel1],dt2[sel2]))
  # ** Create data.frame of 365/6 dates for real decade
  st<-as.POSIXlt(0,origin=paste0(2000+decade*10,"-01-01"),tz="UTC")
  fn<-as.POSIXlt(0,origin=paste0(2000+decade*10+9,"-12-31"),tz="UTC")
  df2<-data.frame(v2="f365",dates=seq(st,fn,by=24*3600))
  dfo<-merge(df1,df2,by="dates",all=TRUE)
  sna<-which(is.na(dfo$v1)==F)
  ao<-array(NA,dim=c(dim(a)[1:2],dim(dfo)[1]))
  ao[,,sna]<-a
  ao<-apply(ao,c(1,2),na.approx)
  ao<-aperm(ao,c(2,3,1))
  # resample and mask out coastal area
  ro<-rast(ao)
  ext(ro)<-ext(r1)
  rte<-rast(landsea)
  crs(ro)<-crs(rte)
  ro<-resample(ro,rte)
  a<-as.array(ro)
  landsea[landsea==0]<-NA
  a<-applymask(a,landsea)
  ro2<-rast(a)
  ext(ro2)<-ext(ro)
  crs(ro2)<-crs(rte)
  return(ro2)
}
# For one tile an one model run and one decade applies function cropandsortUKCPone to all variables
# pathtoUKCP - directory with raw UKCP data for that model run and tile
# pathout - directory in which to save data
# landsea - a landsea raster object of land fractions that must match the extent of the era5 data
# decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
# modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
# data, so file naming and foldr convention assumed to match that of data supplied via dropbox
# saves SpatRaster objects for each variable to disk in directory pathout
.cropandsortUKCPdecade <- function(pathtoUKCP,pathout,landsea,decade,modelrun) {
  mrtxt<-ifelse(modelrun<10,paste0("0",modelrun),paste0("",modelrun))
  suppressWarnings(dir.create(pathout))
  txt1<-paste0(2000+(decade-1)*10+9,"1201")
  txt2<-paste0(2000+decade*10+9,"1130")
  txt3<-paste0(2000+decade*10+9,"1201")
  txt4<-paste0(2000+(decade+1)*10+9,"1130")
  to<-paste0(2000+decade*10,"_",2000+decade*10+9)
  # tmx
  fi1<-paste0(pathtoUKCP,"tasmax/","tasmax_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"tasmax/","tasmax_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"tasmax_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # tmn
  fi1<-paste0(pathtoUKCP,"tasmin/","tasmin_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"tasmin/","tasmin_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"tasmin_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # tas
  fi1<-paste0(pathtoUKCP,"tas/","tas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"tas/","tas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  rot<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"tas_",mrtxt,"_",to,".tif")
  writeRaster(rot,filename=fo,overwrite=T)
  # psl
  fi1<-paste0(pathtoUKCP,"psl/","psl_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"psl/","psl_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  ro<-ro/10
  fo<-paste0(pathout,"psl_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # huss
  fi1<-paste0(pathtoUKCP,"huss/","huss_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"huss/","huss_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"huss_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # u10
  fi1<-paste0(pathtoUKCP,"uas/","uas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"uas/","uas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"uas_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # v10
  fi1<-paste0(pathtoUKCP,"vas/","vas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"vas/","vas_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"vas_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # rss
  fi1<-paste0(pathtoUKCP,"rss/","rss_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"rss/","rss_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  fo<-paste0(pathout,"rss_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
  # skyem
  lwup<-0.97*5.67*10^-8*(rot+273.15)^4
  fi1<-paste0(pathtoUKCP,"rls/","rls_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"rls/","rls_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  lwdn<-lwup+ro
  skyem<-lwdn/lwup
  skyem[skyem>1]<-1
  fo<-paste0(pathout,"skyem_",mrtxt,"_",to,".tif")
  writeRaster(skyem,filename=fo,overwrite=T)
  # pr
  fi1<-paste0(pathtoUKCP,"pr/","pr_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt1,"-",txt2,".nc")
  fi2<-paste0(pathtoUKCP,"pr/","pr_rcp85_land-gcm_global_60km_",mrtxt,"_day_",txt3,"-",txt4,".nc")
  ro<-cropandsortUKCPone(fi1,fi2,landsea,decade)
  ro[ro<0]<-0
  fo<-paste0(pathout,"pr_",mrtxt,"_",to,".tif")
  writeRaster(ro,filename=fo,overwrite=T)
}
# For one tile an one model run extracts the 2018 data and saves this to disk
# pathin - directory with data - should match pathour of cropandsortUKCPdecade function
# pathout - directory in which to save data - usually a subdirectory of pathin
# modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
# data, so file naming and foldr convention assumed to match that of data supplied via dropbox
# saves 2018 data to disk as terra SpatRaster objects
.crop2018UKCP <- function(pathin,pathout,modelrun) {
  saver<-function(varn,pathin,pathout,mrtxt,to,sel) {
    fi<-paste0(pathin,varn,"_",mrtxt,"_",to,".tif")
    fo<-paste0(pathout,varn,"_2018.tif")
    r<-rast(fi)
    a<-as.array(r)
    a<-a[,,sel]
    ro<-rast(a)
    ext(ro)<-ext(r)
    crs(ro)<-crs(r)
    writeRaster(ro,file=fo,overwrite=TRUE)
  }
  # find out wich entries are 2018
  st<-as.POSIXlt(0,origin="2010-01-01",tz="UTC")
  fi<-as.POSIXlt(0,origin="2019-12-31",tz="UTC")
  tme<-as.POSIXlt(seq(st,fi,by=3600*24))
  sel<-which(tme$year+1900==2018)
  # create text strings for filenames
  decade<-1
  mrtxt<-ifelse(modelrun<10,paste0("0",modelrun),paste0("",modelrun))
  suppressWarnings(dir.create(pathout))
  to<-paste0(2000+decade*10,"_",2000+decade*10+9)
  # tmx
  saver("tasmax",pathin,pathout,mrtxt,to,sel)
  saver("tasmin",pathin,pathout,mrtxt,to,sel)
  saver("tas",pathin,pathout,mrtxt,to,sel)
  saver("psl",pathin,pathout,mrtxt,to,sel)
  saver("huss",pathin,pathout,mrtxt,to,sel)
  saver("uas",pathin,pathout,mrtxt,to,sel)
  saver("vas",pathin,pathout,mrtxt,to,sel)
  saver("rss",pathin,pathout,mrtxt,to,sel)
  saver("skyem",pathin,pathout,mrtxt,to,sel)
  saver("pr",pathin,pathout,mrtxt,to,sel)
}
# ============================================================================ #
# ~~~~~~~~~~~~~ Temporal downscale worker functions here ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# Writes a ncd4 file to disk
# NB - we might want to document this function and move it out of workers??
# Inputs:
# r - a terra::SpatRaster object of data to be written out
# year = the year for which data are to be written out (used for time stamp)
# varn - climate variable name (e.g. tasmax)
# varl - climate variable long name (Daily maximum temperature at 2m)
# unit - units of varn (e.g. deg C)
# asinteger - optional lgoical indicating whether values should be stored as integers (dee details)
# mult - multiplier to use if data are stored as integers (dee details)
# Details:
# Setting asintiger to TRUE roughly halves the size of the file written out but
# reduces precision. E.g If mult = 100, a temperature value of 9.87694 would be
# written out as 988. The function automatically corrects the units written to
# the nc file - e.g. if mult is 100 and unit is 'deg C', 'deg C * 100' is written
# out as the unit in the nc file
.writenc<-function(r,fo,year,varn,varl,unit,asinteger= TRUE,mult=100) {
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
  tme<-as.POSIXlt(seq(st,fn,by=3600))
  terra::time(r)<-tme
  writeCDF(r,fo,overwrite=TRUE,compression=9,varname=varn,
           longname=varl,unit=unit)
}
#' terra version of microclima get_dem
#' require(elevatr)
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
