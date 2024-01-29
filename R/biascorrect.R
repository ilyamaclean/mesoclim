# Applies bias correction to one variable
# rera - a SpatRaster object of era5 data for 2018
# rukcpf - a SpatRaster object of ukcp data for 2018
# rukcpa   a SpatRaster object of ukcp data to which correction is to be applied - typically a decade's worth of data
# - returns a bias corrected SpatRaster object
# NB suggest including function call to rainadjust to make below properly generic (and
# suggest calling it precipadjust).
biascorrectukcpone <- function(rera, rukcpf, rukcpa, rangelims = NA) {
  a1<-as.array(rera)
  a2<-as.array(rukcpf)
  a3<-as.array(rukcpa)
  msk1<-apply(a1,c(1,2),mean,na.rm=T)
  msk2<-apply(a2,c(1,2),mean,na.rm=T)
  msk<-msk1*msk2
  ao<-array(NA,dim=dim(a3))
  for (i in 1:dim(a1)[1]) {
    for (j in 1:dim(a1)[2]) {
      if (is.na(msk[i,j]) == F) {
        v1<-a1[i,j,]
        v2<-a2[i,j,]
        v3<-a3[i,j,]
        v1 <- v1[order(v1)]
        v2 <- v2[order(v2)]
        m1 <- gam(v1~s(v2))
        xx <- predict.gam(m1, newdata = data.frame(v2 = v3))
        if (is.na(rangelims) == F) {
          me <- mean(v1)
          mx <- me + rangelims * (max(v1) - me)
          mn <- me - rangelims * (me - min(v1))
          xx[xx>mx]<-mx
          xx[xx<mn]<-mn
        }
        ao[i,j,]<-xx
      }
    }
  }
  ro<-rast(ao)
  ext(ro)<-ext(rukcpa)
  crs(ro)<-crs(rukcpa)
  return(ro)
}
# Adjust a vector of daily precipitation to ensure bias correction
# v1 - era5 2018 data
# v2 - ukcp 2018 data
# v3 - ucp data to which correction is applied
# returns a corrected vector of rainfall data of length mathcing v3
rainadjust<-function(v1, v2, v3) {
  raindays<-function(v) {
    u<-v>0
    length(u[u==TRUE])
  }
  ma <- function(x){
    y<-filter(x, rep(1 / 10, 10), sides = 2)
    y[1:5]<-y[6]
    n<-length(y)
    y[(n-5):n]<-y[n-6]
    y
  }
  v1[v1<0.01]<-0
  v2[v2<0.01]<-0
  v3[v3<0.01]<-0
  rd1<-raindays(v1)
  rd2<-raindays(v2)
  rd3<-raindays(v3)
  if (rd1 > rd2) {
    # find out how many extra raindays need asigning
    rdextra<-round(rd3*(rd1/rd2-1),0)
    # find out which dry days are most like to have rain
    rav<-ma(v3)
    rav[v3==0]<-0
    o<-order(rav,decreasing=TRUE)[1:rdextra]
    v3[o]<-rav[o]
  } else {
    rdextra<-round(rd3*(rd2/rd1-1),0)
    # get rid of driest rain days
    sel<-which(v3>0)
    o<-order(v3[sel])[1:rdextra]
    v3[sel[o]]<-0
  }
  rt1<-sum(v1)
  rt2<-sum(v2)
  v3<-v3*(sum(v1)/sum(v2))
  return(v3)
}
# Applies the rainadjust to SpatRaster arrays of data
# rera - a SpatRaster object of era5 2018 data
# rukcpf - a SpatRaster object of ukcp 2018 data
# rukcpa - a SpatRaster object of ukcp data to which the correction is to be applied
# returns a bias corrected SpatRaster object of rainfall
precipcorrect<-function(rera, rukcpf, rukcpa) {
  a1<-as.array(rera)
  a2<-as.array(rukcpf)
  a3<-as.array(rukcpa)
  msk<-apply(a1,c(1,2),mean,na.rm=T)
  ao<-array(NA,dim=dim(a3))
  for (i in 1:dim(a1)[1]) {
    for (j in 1:dim(a1)[2]) {
      if (is.na(msk[i,j]) == F) {
        ao[i,j,]<-rainadjust(a1[i,j,],a2[i,j,],a3[i,j,])
      }
    }
  }
  ro<-rast(ao)
  ext(ro)<-ext(rukcpa)
  crs(ro)<-crs(rukcpa)
  return(ro)
}
# Saves a SpatRaster object as an nc file
# r - a SpatRaster object
# baseyear - used for adding time stamp to nc file (asumes daily)
# shortname - Name of the dataset
# longname - Long name of the dataset
# unit - character. Unit of the data
# fileout - name of file to be saved (includign path)
# saves a compressed netCDF4 file to disk
savenc<-function(r,baseyear,shortname,longname,unit,fileout) {
  terra::time(r)<-as.Date(paste0(baseyear,"-01-01")) + 0:(nlyr(r)-1)
  writeCDF(r,fileout,overwrite=TRUE,compression=9,varname=shortname,
           longname=longname,unit=unit)
}
# Applies biascorrectukcpone (precipcorrect) to all ukcp variables for one tile and model run
# and saves data to disk as nc files per decade
# pathtoera - directory with 2018 daily era5 data (as returned by era5todaily)
# pathtoukcp18 - directory with 2018 ukcp data (as returned by crop2018UKCP)
# pathtoukcpdecade - directory with decadal ukcp data (as returned by cropandsortUKCPdecade)
# pathout - directory in which to save corrected data
# decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
# modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
# data, so file naming and foldr convention assumed to match that of data supplied via dropbox
# saves bias corrected daily ukcp data as compressed ncdf4 files
biascorrectukcpall<-function(pathtoera,pathtoukcp18,pathtoukcpdecade,pathout,decade,modelrun) {
  # file names
  mtxt<-ifelse(modelrun<10,paste0("0",modelrun),paste0("",modelrun))
  to<-paste0(2000+decade*10,"_",2000+decade*10+9)
  byr<-2000+decade*10
  # Temperature
  # ** era5
  fi<-paste0(pathtoera,"tasmax.tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoera,"tasmin.tif")
  tmn<-rast(fi)
  edtr<-tmx-tmn
  etme<-(tmx+tmn)/2
  # ** ukcp fit
  fi<-paste0(pathtoukcp18,"tasmax_2018.tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoukcp18,"tasmin_2018.tif")
  tmn<-rast(fi)
  udtr<-tmx-tmn
  utme<-(tmx+tmn)/2
  # ** ukcp apply
  fi<-paste0(pathtoukcpdecade,"tasmax_",mtxt,"_",to,".tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"tasmin_",mtxt,"_",to,".tif")
  tmn<-rast(fi)
  uadtr<-tmx-tmn
  uatme<-(tmx+tmn)/2
  dtrc<-biascorrectukcpone(edtr, udtr, uadtr, rangelims = 1.1)
  tmec<-biascorrectukcpone(etme, utme, uatme, rangelims = 1.1)
  tmxc<-tmec+0.5*dtrc
  tmnc<-tmec-0.5*dtrc
  fo1<-paste0(pathout,"tasmax_",mtxt,"_",to,".nc")
  fo2<-paste0(pathout,"tasmin_",mtxt,"_",to,".nc")
  savenc(tmxc,byr,"tasmax","Maximum air temperature at 2 m","deg C",fo1)
  savenc(tmnc,byr,"tasmin","Minimum air temperature at 2 m","deg C",fo2)
  # Pressure
  fi<-paste0(pathtoera,"psl.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"psl_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"psl_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  psl<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  fo<-paste0(pathout,"psl_",mtxt,"_",to,".nc")
  savenc(psl,byr,"psl","Surface level pressure","kPa",fo)
  # Humidity
  fi<-paste0(pathtoera,"huss.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"huss_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"huss_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  huss<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  h<-as.array(huss)
  relh<-suppressWarnings(converthumidity(h, intype = "specific", outtype = "relative",
                                         tc = as.array(tmec), pk = as.array(psl)))
  relh[relh>100]<-100
  relh[relh<20]<-20
  relh<-rast(relh)
  ext(relh)<-ext(huss)
  crs(relh)<-crs(huss)
  fo<-paste0(pathout,"relhum_",mtxt,"_",to,".nc")
  savenc(relh,byr,"relhum","Relative humidity at 2 m","Percentage",fo)
  # Wind speed and direction
  # ** u wind
  ufi<-paste0(pathtoera,"uas.tif")
  ur1<-rast(ufi)
  ufi<-paste0(pathtoukcp18,"uas_2018.tif")
  ur2<-rast(ufi)
  ufi<-paste0(pathtoukcpdecade,"uas_",mtxt,"_",to,".tif")
  ur3<-rast(ufi)
  # ** v wind
  vfi<-paste0(pathtoera,"vas.tif")
  vr1<-rast(vfi)
  vfi<-paste0(pathtoukcp18,"vas_2018.tif")
  vr2<-rast(vfi)
  vfi<-paste0(pathtoukcpdecade,"vas_",mtxt,"_",to,".tif")
  vr3<-rast(vfi)
  # ** Calculate wind speed
  ws1<-sqrt(ur1^2+vr1^2)
  ws2<-sqrt(ur2^2+vr2^2)
  ws3<-sqrt(ur3^2+vr3^2)
  wsc<-biascorrectukcpone(ws1, ws2, ws3, rangelims = 1.1)
  # ** Adjust to 2 m
  wsc<-wsc*(4.87/log(67.8*10-5.42))
  wsc[wsc<0.001]<-0.001
  wdir<-(atan2(ur3,vr3)*(180/pi)+180)%%360
  # ** Save files
  fo<-paste0(pathout,"winds_",mtxt,"_",to,".nc")
  savenc(wsc,byr,"winds","Wind speed at 2 m","m/s",fo)
  fo<-paste0(pathout,"windd_",mtxt,"_",to,".nc")
  savenc(wdir,byr,"windd","Wind direction","degrees",fo)
  # ** Sky emmisivity
  fi<-paste0(pathtoera,"skyem.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"skyem_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"skyem_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  skyem<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  skyem[skyem>1]<-1
  skyem[skyem<0]<-0
  fo<-paste0(pathout,"skyem_",mtxt,"_",to,".nc")
  savenc(skyem,byr,"skyem","Sky emissivity","0-1",fo)
  # ** Shortwave radiation
  fi<-paste0(pathtoera,"rss.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"rss_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"rss_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  rss<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  rss[rss<0]<-0
  fo<-paste0(pathout,"rss_",mtxt,"_",to,".nc")
  savenc(rss,byr,"rss","Total downward shortwave radiation","W/m**2",fo)
  # ** Precipitation
  fi<-paste0(pathtoera,"pr.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"pr_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"pr_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  pr<-precipcorrect(r1,r2,r3)
  fo<-paste0(pathout,"pr_",mtxt,"_",to,".nc")
  savenc(pr,byr,"pr","Precipitation rate","mm/day",fo)
}
