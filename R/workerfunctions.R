# ******** NB - suggest putting any functions we don't want to make
# available as stand-alone documented functions in here. Precede by dot as R
# studio won't then expect oxygen2 when documenting package
# ============================================================================ #
# ~~~~~~~~~~~~~ Data handing worker functions here ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
# ** Following is a bit of a code dump. We may not need it all:
# NB:
#  ** (1) For several of these functions we'll need to add the appropriate imports
#     (e.g. raster terra) before they work (see example for .cropnc).
#  ** (2) it would be good to get rid of raster as a dependency as I suspect it
#         will be deprecated owing to rgdal dependency, but for now brick seems
#         to be faster than terra equivalent
#  ** (3) I've preceded all by a dot. But I wrote them before this was done, so they
#          they might still call functions without the dot.
# ============================================================================= #
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
# Crops a large ERA5 nc file into smaller e.g. subtiles in a compuationally efficient manner
# filein - filename of nc file to be chopped
# ecrop - raster extent object of area to which nc file is to be cropped
# fileout - filename (and path) of chopped nc file
# This function writes the chopped nc file to disk
#' @import terra, raster
.cropera5nc<-function(filein, ecrop, fileout) {
  vars<-c("t2m","d2m","sp","u10","v10","tp","tcc","msnlwrf","msdwlwrf","fdir","ssrd")
  unit<-c("K","K","Pa","m s**-1","m s**-1","m","0 - 1","W m**-2","W m**-2","J m**-2","J m**-2")
  lnames<-c("2 metre temperature","2 metre dewpoint temperature","Surface pressure",
            "10 metre U wind component","10 metre V wind component","Total precipitation",
            "Total cloud cover","Mean surface net long-wave radiation flux",
            "Mean surface downward long-wave radiation flux","Total sky direct solar radiation at surface",
            "Surface solar radiation downwards")

  # (1)
  r1<-brick(filein,varname=vars[1])
  r1<-crop(r1,ecrop)
  r1<-rast(r1)
  # (2)
  r2<-brick(filein,varname=vars[2])
  r2<-crop(r2,ecrop)
  r2<-rast(r2)
  # (3)
  r3<-brick(filein,varname=vars[3])
  r3<-crop(r3,ecrop)
  r3<-rast(r3)
  # (4)
  r4<-brick(filein,varname=vars[4])
  r4<-crop(r4,ecrop)
  r4<-rast(r4)
  # (5)
  r5<-brick(filein,varname=vars[5])
  r5<-crop(r5,ecrop)
  r5<-rast(r5)
  # (6)
  r6<-brick(filein,varname=vars[6])
  r6<-crop(r6,ecrop)
  r6<-rast(r6)
  # (7)
  r7<-brick(filein,varname=vars[7])
  r7<-crop(r7,ecrop)
  r7<-rast(r7)
  # (8)
  r8<-brick(filein,varname=vars[8])
  r8<-crop(r8,ecrop)
  r8<-rast(r8)
  # (9)
  r9<-brick(filein,varname=vars[9])
  r9<-crop(r9,ecrop)
  r9<-rast(r9)
  # (10)
  r10<-brick(filein,varname=vars[10])
  r10<-crop(r10,ecrop)
  r10<-rast(r10)
  # (11)
  r11<-brick(filein,varname=vars[11])
  r11<-crop(r11,ecrop)
  r11<-rast(r11)
  # Sort out times
  nc<-nc_open(filein)
  tme<-ncvar_get(nc,var="time")
  nc_close(nc)
  tme<-as.POSIXlt(tme*3600,origin="1900-01-01 00:00", tz = "UTC")
  terra::time(r1)<-tme
  terra::time(r2)<-tme
  terra::time(r3)<-tme
  terra::time(r4)<-tme
  terra::time(r5)<-tme
  terra::time(r6)<-tme
  terra::time(r7)<-tme
  terra::time(r8)<-tme
  terra::time(r9)<-tme
  terra::time(r10)<-tme
  terra::time(r11)<-tme
  # stick together
  a<-sds(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11)
  names(a) <- vars
  terra::longnames(a)<-lnames
  terra::units(a)<-unit
  writeCDF(a, filename = fileout, compression = 9, overwrite=TRUE)
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
