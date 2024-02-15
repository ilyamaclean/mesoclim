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
# Produces a matrix of latitudes form a terra::SpatRaster object
# Inputs:
# r - a terra::SpatRaster object
# Returns a matrix of latidues
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
# Produces a matrix of longitudes form a terra::SpatRaster object
# Inputs:
# r - a terra::SpatRaster object
# Returns a matrix of longitudes
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
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
#' Simulate cloud patchiness
#' @import gstat
.simpatch<-function(dtmf,n,mn=0.5,mx=2) {
  reso<-res(dtmf)[1]
  if (reso <= 10) dtmc<-aggregate(dtmf,100)
  if (reso <= 100) dtmc<-aggregate(dtmf,10)
  reso2<-res(dtmc)[1]
  sill<-1000/reso2
  rge<-sill
  xy <- expand.grid(1:dim(dtmc)[2], 1:dim(dtmc)[1])
  names(xy) <- c('x','y')
  g1 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0,
              model=vgm(psill = sill, range = rge, nugget = 3, model='Sph'), nmax = 40)
  yy1 <- predict(g1, newdata=xy, nsim=n)
  r<-rast(yy1)
  # adjust to required range
  rge<-max(as.vector(r))-min(as.vector(r))
  nrge<-log(mx)-log(mn)
  mu<-nrge/rge
  r<-r*mu
  r<-exp(r-mean(r))
  ext(r)<-ext(dtmc)
  crs(r)<-crs(dtmc)
  r<-resample(r,dtmf)
  return(r)
}

# basin code
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
# ============================================================================ #
# ~~~~~~~~~ Basin delineation worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
.basindelin<-function(dtm) {
  # (1) stick buffer around dtm
  dm<-as.matrix(dtm,wide=TRUE)
  dm2<-array(NA,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # (2) create blank basin file
  bsn<-dm2*NA
  dun<-array(1,dim=dim(bsn))
  # Iterate to delineate all basins
  tsta<-1
  bn<-1 # basin number
  while (tsta == 1) {
    # First iteration
    mu<-bsn*NA
    mu[is.na(bsn)]<-1
    s<-suppressWarnings(which(dm2*mu*dun==min(dm2*mu*dun,na.rm=T),arr.ind=T))
    if (length(s) > 0) {
      rw<-as.numeric(s[1,1])
      cl<-as.numeric(s[1,2])
      bsn[rw,cl]<-bn
      dun[rw,cl]<-NA
      # (4) select 3 x 3 matrix around grid cell
      m3<-dm2[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]
      b3<-bsn[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]
      # (5) identify  which grid cells in 3 x 3 undone and higher and assign to basin i
      sb<-which(m3>m3[2,2])
      b3[sb]<-bn
      bsn[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]<-b3
    } else tsta<-0
    # Subsequent iterations
    if (tsta > 0) {
      tst<-1
      while (tst == 1) {
        mu<-bsn*NA
        mu[bsn==bn]<-1
        s<-suppressWarnings(which(dm2*mu*dun==min(dm2*mu*dun,na.rm=T),arr.ind=T))
        if (length(s) > 0) {
          rw<-as.numeric(s[1,1])
          cl<-as.numeric(s[1,2])
          bsn[rw,cl]<-bn
          dun[rw,cl]<-NA
          # (4) select 3 x 3 matrix around grid cell
          d3<-dun[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]
          m3<-dm2[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]*d3
          b3<-bsn[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]
          m3[2,2]<-dm2[rw,cl]
          # (5) identify  which grid cells in 3 x 3 undone and higher and assign to basin bn
          sb<-which(m3>m3[2,2])
          b3[sb]<-bn
          bsn[(c(rw-1):(rw+1)),(c(cl-1):(cl+1))]<-b3
        } else tst<-0
      }
      bn<-bn+1  # go on to next basin
    }
  }
  dd<-dim(bsn)
  bsn<-bsn[2:(dd[1]-1),2:(dd[2]-1)]
  r<-.rast(bsn,dtm)
  return(r)
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
  # ********************** Mosaic and renumber ******************************* #
  bout<-mosaic(b1n,b2n)
  m<-.is(bout)
  # renumber basins
  u<-unique(as.vector(m))
  u<-u[is.na(u) == FALSE]
  u<-u[order(u)]
  for (i in 1:length(u)) m[m==u[i]]<-i
  bout<-.rast(m,bout)
  return(bout)
}
#' Do an entire column of tiled basins
.docolumn<-function(dtm,tilesize,boundary,x) {
  e<-ext(dtm)
  reso<-res(dtm)
  ymxs<-as.numeric(ceiling((e$ymax-e$ymin)/reso[2]/100))-1
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
    xmn<-as.numeric(e$xmin)+reso[1]*100*x
    xmx<-xmn+reso[1]*100
    ymn<-as.numeric(e$ymin)+reso[2]*100*y
    ymx<-ymn+reso[2]*100
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
#' @import raster, terra
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
# ============================================================================ #
# ~~~~~~~~~~~~~~ Spatial downscale worker functions here ~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
