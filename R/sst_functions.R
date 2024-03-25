#####################################################################
# Alternatives to SSTinterpolate
# spatial_interpNA is a direct copy of SSTinterpolate code
# Temporal interpolation may be limited in terms of length of tmeout as uses terra::approximate
#####################################################################
#' Spatially interpolate missing values in Spatraster
#'
#' @param r 2 or 3D Spatraster
#'
#' @return Spatraster of same dim as r but without NA
#' @export
#' @import terra
#' @examples
#' r <- rast(ncols=5, nrows=5, vals=rep(c(1,2,NA,4,5),5))
#' rout<-spatial_interpNA(r)
#' plot(c(r,rout))
spatial_interpNA<-function(r){
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
      me<-apply(.is(r),3,mean,na.rm=T)
      m2<-.vta(me,r[[1]])
      m1[s]<-m2[s]
    }
    rout<-.rast(m1,r)
  }
  return(rout)
}
#' Temporal interpolation of missing times in Spatraster
#'  Works for any time step not just daily or hourly eg monthly SST
#' @param rin = 3D Spatraster
#' @param tmeout = Date or POSIX vector of datetimes required
#' @param tmein = Date or POSIX vector of datetimes of rin - if NA assumes time(rin) provides valules
#'
#' @return Spatraster stack of nlyrs=length(tmeout)
#' @export
#' @import terra
#' @import lubridate
#' @examples
#' r <- rast(ncols=5, nrows=5, vals=rep(1,25))
#' rin<-c(r,r*2,r*3,r*4,r*5)
#' tmein<-as.Date(seq(17532,17652,30))
#' terra::time(rin)<-tmein
#' tmeout<-as.Date(seq(17532,17651,1))
#' rout<-time_interp(rin,tmeout, NA)
#' plot(rout[[c(1,15,31,95,119,120)]])

time_interp<-function(rin, tmeout, tmein=NA){
  # Check time inputs - converts all date/datetime values to POSIXlt
  if (class(tmein)[1]=='logical') tmein<-as.POSIXlt(terra::time(rin))
  if(class(tmein)[1]!='POSIXlt') tmein<-as.POSIXlt(tmein)
  if(class(tmeout)[1]!='POSIXlt') tmeout<-as.POSIXlt(tmeout)
  if (length(tmein) != dim(rin)[3]) stop("length of tmein must match dim 3 of rin")
  if (tmeout[length(tmeout)] < tmein[1] || tmeout[1] > tmein[length(tmein)]) stop ("tmeout and tmein do NOT overlap!!")

  # Find time step (in sec)
  if(length(unique(diff(tmeout)))>1) stop("tmeout does NOT have regular time steps!!") else tstep<-unique(as.numeric(diff(tmeout),units='secs'))

  # Extend desired output timeseries to next available input data before/after output timeseries
  if(any(tmein<tmeout[1])) tstart<- tmein[which(tmein<tmeout[1])][1] else tstart<-as.POSIXlt(tmeout[1])
  if(any(tmein>tmeout))  tend<-tmein[which(tmein>tmeout[length(tmeout)])][1] else tend<-as.POSIXlt(tmeout[length(tmeout)])
  ext_tmeout<-as.POSIXlt(seq(tstart,tend,tstep))

  # Create output raster
  rout<-rast(rin[[1]], nl=length(ext_tmeout),vals=NA)
  terra::time(rout)<-ext_tmeout

  # Add rin layers to rout
  selout<-which(ext_tmeout %in% tmein)
  selin<-which(tmein %in% ext_tmeout)
  if(length(selin)>0) rout[[selout]]<-rin[[selin]]

  # Interpolate linearly other times
  rout<-approximate(rout,rule=2)

  # Return only values in tmeout
  sel<-which(as.POSIXlt(terra::time(rout)) %in% tmeout)
  return(rout[[sel]])
}
