#' @title Applies bias correction to one climate variable
#' @description The function `biascorrect` applies bias correction to spatial
#' climate dataset
#' @param hist_obs a stacked SpatRaster of climate observations.
#' @param hist_mod a stacked SpatRaster of modelled climate data for
#' for the same period as `hist_obs`. Reprojected and resampled to match `hist_obs`
#' if necessary.
#' @param fut_mod optionally, SpatRaster of modelled climate data for
#' for e.g. a future period to which corrections are applied. Reprojected and
#' resampled to match `hist_obs` if necessary. If not supplied and `mod_out = FALSE`
#' corrections are applied to `hist_mod`
#' @param mod_out optional logical indicating whether to return an object
#' of class `biascorrectmodels`, which allows the same form of bias correction
#' to be applied to datasets over multiple time periods using [biascorrect_apply()].
#' @param rangelims fractional amount above or below the range of values in
#' `hist_obs` by which values in the returned dataset can vary. I.e. if set to
#' 1.1, values in any given grid cell of the returned dataset cannot exceed
#' 1.1 times the range of values in the corresponding grid cell in `hist_obs`.
#' If not supplied (the default), returned values are unbounded.
#' @return if `mod_out = FALSE`, a SpatRaster of bias corrected data. If
#' `mod_out = TRUE`  an object of class `biascorrectmodels`,
#' @import terra
#' @import mgcv
#' @export
#' @keywords biascorrect
#' @rdname biascorrect
#' @seealso [precipcorrect()] for applying corrections to precipitation data and
#' [biascorrect_apply()] for applying corrections to multiple datasets.
biascorrect <- function(hist_obs, hist_mod, fut_mod = NA, mod_out = FALSE, rangelims = NA) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (mod_out == FALSE & inherits(fut_mod, "logical")) fut_mod<-hist_mod
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (mod_out == FALSE) if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  if (mod_out == FALSE) fut_mod<-resample(fut_mod,hist_obs)
  # Convert to arrays
  a1<-.is(hist_obs)
  a2<-.is(hist_mod)
  if (mod_out == FALSE) a3<-.is(fut_mod)
  # Mask out any cells that are missing
  msk1<-apply(a1,c(1,2),mean,na.rm=T)
  msk2<-apply(a2,c(1,2),mean,na.rm=T)
  if (mod_out == FALSE) {
    msk3<-apply(a3,c(1,2),mean,na.rm=T)
    msk<-msk1*msk2*msk3
  } else msk<-msk1*msk2
  # Check whether dataset has more than 1000 entries per time-series and
  # subset if it does
  n<-dim(a1)[3]
  if (n > 1000) s<-sample(0:n,998,replace = FALSE)  # 998 so min and max can be tagged on
  if (mod_out) {
    mods<-list()
    m<-a1[,,1]
  }  else {
    ao<-array(NA,dim=dim(a3))
  }
  # Create array for storing data
  counter<-1
  for (i in 1:dim(a1)[1]) {
    for (j in 1:dim(a1)[2]) {
      if (is.na(msk[i,j]) == F) {
        v1<-a1[i,j,]
        v2<-a2[i,j,]
        v1 <- v1[order(v1)]
        v2 <- v2[order(v2)]
        if (n > 1000) {
          v1<-c(v1[1],v1[s],v1[length(v1)]) # tags on min and max value
          v2<-c(v2[1],v2[s],v2[length(v2)]) # tags on min and max value
        }
        # Apply gam
        m1 <- gam(v1~s(v2))
        if (mod_out) {
          m[i,j]<-counter
          mods[[counter]]<-m1
        } else {
          v3<-a3[i,j,]
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
  }
  if (mod_out) {
    rindex<-.rast(m,hist_obs[[1]])
    rindex<-wrap(rindex)
    ao<-list(models=mods,rindex=rindex)
    class(ao) <- "biascorrectmodels"
  } else {
    ao<-.rast(ao,hist_obs)
  }
  return(ao)
}
#' @title Applies bias correction model to climate variable
#' @description The function `biascorrect_apply` applies a set of bias corrections
#' derived using `bias_correct` to a spatial climate dataset
#' @param fut_mod a SpatRaster of modelled climate data for for e.g. a future
#' period to which corrections are applied.
#' @param biasmods an object of class `biascorrectmodels` derived using [biascorrect()].
#' @return a SpatRaster matching `fut_mod` with climate data bias corrected.
#' @details The `biascorrect_apply` function allows the same bias correction to be applied
#' to multiple datasets, spanning e.g. several future time periods or climate scenarious.
#' The spatial extent of `fut_mod` must at least partially overlap with that of
#' `hist_obs` when running `biascorrect` as corrected values are only returned for the overlapping area
#' @import terra
#' @import mgcv
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#' @keywords biascorrect
#' @rdname biascorrect_apply
#' @seealso [biascorrect()] for deriving `biasmods` and `precip_correct` for
#' applying corrections to precipitation datasets.
biascorrect_apply<-function(fut_mod, biasmods, rangelims = NA) {
  if (!inherits(fut_mod, "SpatRaster")) stop("fut_mod must be a SpatRaster")
  if (!inherits(biasmods, "biascorrectmodels")) stop("biasmods must be an object of class biascorrectmodels")
  r<-rast(biasmods$rindex)
  mods<-biasmods$models
  if (crs(fut_mod) != crs(r)) fut_mod<-project(fut_mod,r)
  fut_mod<-resample(fut_mod,r)
  # Populate m wiht i and j
  ri<-.rast(populatematrix(.is(r)*0, "i"),r)
  rj<-.rast(populatematrix(.is(r)*0, "j"),r)
  # do xy match
  xy1<-data.frame(xyFromCell(r, 1:ncell(r)))
  index<-as.vector(extract(r,xy1)[,2])
  i<-as.vector(extract(ri,xy1)[,2])
  j<-as.vector(extract(rj,xy1)[,2])
  xy1<-cbind(xy1,index,i,j)
  xy2<-data.frame(xyFromCell(fut_mod, 1:ncell(fut_mod)))
  xy1<-merge(xy1,xy2,by=c("x","y"),all=TRUE)
  ao<-.is(fut_mod)
  for (i in 1:dim(ao)[1]) {
    for (j in 1:dim(ao)[[2]]) {
      s<-which(xy1$i == i & xy1$j == j)
      index<-xy1[s,3]
      if (is.na(index) == FALSE) {
        v3<-ao[i,j,]
        ao[i,j,]<-predict.gam(mods[[index]], newdata = data.frame(v2 = v3))
      } else ao[i,j,]<-NA
    }
  }
  ao<-.rast(ao,fut_mod[[1]])
  ao<-mask(ao,r)
  return(ao)
}
#' @title Applies bias correction to precipitation data
#' @description The function `precipcorrect` applies bias correction to spatial
#' climate dataset
#' @param hist_obs a stacked SpatRast of precipitation data.
#' @param hist_mod a stacked SpatRast of modelled precipitation data for
#' for the same period and spatial area as `hist_obs`. Reprojected and resampled
#' to match `hist_obs` if necessary.
#' @param fut_mod a stacked SpatRast of modelled precipitation data for
#' for e.g. a future period to which corrections are applied. eprojected and resampled
#' to match `hist_obs` if necessary.
#' @param mod_out optional logical indicating whether to return model correction
#' coefficients, which can be applied to datasets over multiple time periods
#' using [precipcorrect_apply()].
#' @param rangelim fractional amount above the maximum values of precipitation
#' in `hist_obs` in any given cell by which values in the returned dataset can
#' attain. I.e. if set to 2, values in any given grid cell of the returned dataset
#' cannot exceed twice the maximum amount in the corresponding grid cell in
#' `hist_obs`.If not supplied (the default), returned values are unbounded.
#' @return If `mod_out = FALSE` (the default), a SpatRaster of bias corrected
#' precipitation. If `mod_out = TRUE` a list containing two wrapped SpatRasters
#' representing the amount by which to adjust precipitation totals and precipitation
#' day fractions.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#' @keywords biascorrect
#' @rdname precipcorrect
#' @seealso [biascorrect()] for applying corrections to other climate variables and
#' [precipcorrect_apply()] for applying correction coefficients to multiple datasets.
precipcorrect <- function(hist_obs, hist_mod, fut_mod, mod_out = FALSE, rangelim = NA) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (mod_out == FALSE & inherits(fut_mod, "logical")) fut_mod<-hist_mod
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (mod_out == FALSE) if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  if (mod_out == FALSE) fut_mod<-resample(fut_mod,hist_obs)
  # mask data
  hist_mod<-mask(hist_mod,hist_obs)
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_obs
  rcount[rcount > 0] <-1
  rtot1<-apply(.is(hist_obs),c(1,2),sum)
  tfrac1<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_mod
  rcount[rcount > 0] <-1
  rtot2<-apply(.is(hist_mod),c(1,2),sum)
  tfrac2<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate ratios
  mu_tot<-rtot1/rtot2
  mu_frac<-tfrac1/tfrac2
  if (mod_out) {
    r1<-.rast(mu_tot,hist_mod[[1]])
    r2<-.rast(mu_frac,hist_mod[[1]])
    r1<-wrap(r1)
    r2<-wrap(r2)
    out<-list(mu_tot=r1,mu_frac=r2)
  } else {
    # Apply correction to modelled future data
    rcount<-fut_mod
    rcount[rcount > 0] <-1
    rtot<-apply(.is(fut_mod),c(1,2),sum)*mu_tot
    tfrac<-(apply(.is(rcount),c(1,2),sum)/dim(rcount)[3])*mu_frac
    # Calculate regional precipitation
    rrain<-apply(.is(fut_mod),3,sum,na.rm=TRUE)
    rr2<-as.numeric(.mav(rrain,10))
    s<-which(rrain==0)
    rrain[s]<-rr2[s]
    s<-which(rrain==0)
    rrain[s]<-0.1
    rrain<-(rrain/max(rrain))*10
    # Adjust rainfall
    fut_mod<-mask(fut_mod,hist_obs)
    a<-as.array(fut_mod)
    mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
    rtot<-.rast(rtot,hist_obs)
    tfrac<-.rast(tfrac,hist_obs)
    rtot<-as.vector(t(rtot))
    rfrac<-as.vector(t(tfrac))
    mm<-rainadjustm(mm,rrain,rfrac,rtot)
    a2<-array(mm,dim=dim(a))
    if (is.na(rangelim) == F) {
      rmax<-apply(.is(hist_obs),c(1,2),max)
      rmax<-.rta(rmax*rangelim,dim(a2)[3])
      s<-which(a2>rmax)
      a2[s]<-rmax[s]
    }
    # convert to SpatRast
    out<-.rast(a2,fut_mod)
  }
  return(out)
}
#' @title Applies bias correction model to precipitation data
#' @description The function `precipcorrect_apply` applies a set of bias corrections
#' derived using `precip_correct` to a spatial climate dataset
#' @param fut_mod a stacked SpatRast of modelled precipitation data for
#' for e.g. a future period to which corrections are applied. reprojected and
#' resampled to match data held in `mod_out` if necessary.
#' @param a list containing two wrapped SpatRasters derived using [precipcorrect()]
#' representing the amount by which to adjust rainfall totals and rainfall days,
#' which is used as a basis for adjusting individual rainfall events. `
#' @return a SpatRaster matching `fut_mod` with climate data bias corrected.
#' @details The `precipcorrect_apply` function allows the same bias correction to
#' be applied to multiple datasets, spanning e.g. several future time periods or
#' climate scenarious. The spatial extent of `fut_mod` must at least partially
#' overlap with that of `hist_obs` when running `precipcorrect` as corrected values are only returned for the overlapping area
#' @import terra
#' @import mgcv
#' @importFrom Rcpp sourceCpp
#' @useDynLib mesoclim, .registration = TRUE
#' @export
#' @keywords biascorrect
#' @rdname precipcorrect_apply
#' @seealso [precipcorrect()] for deriving `biasmods` and `biascorrect_apply` for
#' applying corrections to other climate variables
precipcorrect_apply<-function(fut_mod, biasmods) {
  if (!inherits(fut_mod, "SpatRaster")) stop("fut_mod must be a SpatRaster")
  mu_tot<-rast(biasmods$mu_tot)
  mu_frac<-rast(biasmods$mu_frac)
  # reproject and crop if necessary
  if (crs(fut_mod) != crs(mu_tot)) fut_mod<-project(fut_mod,mu_tot)
  fut_mod<-resample(fut_mod,mu_tot)
  # Calculate and adjust rainfall total and rain day frac
  rcount<-fut_mod
  rcount[rcount > 0] <-1
  rtot<-apply(.is(fut_mod),c(1,2),sum)*.is(mu_tot)
  tfrac<-(apply(.is(rcount),c(1,2),sum)/dim(rcount)[3])*.is(mu_frac)
  # Calculate regional precipitation
  rrain<-apply(.is(fut_mod),3,sum,na.rm=TRUE)
  rr2<-as.numeric(.mav(rrain,10))
  s<-which(rrain==0)
  rrain[s]<-rr2[s]
  s<-which(rrain==0)
  rrain[s]<-0.1
  rrain<-(rrain/max(rrain))*10
  # Adjust rainfall
  fut_mod<-mask(fut_mod,hist_obs)
  a<-as.array(fut_mod)
  mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
  rtot<-.rast(rtot,hist_obs)
  tfrac<-.rast(tfrac,hist_obs)
  rtot<-as.vector(t(rtot))
  rfrac<-as.vector(t(tfrac))
  mm<-rainadjustm(mm,rrain,rfrac,rtot)
  a2<-array(mm,dim=dim(a))
  # convert to SpatRast
  out<-.rast(a2,fut_mod)
  return(out)
}


#' @title Applies bias correction to UK era5 temperature data
#' @description The function `correct_era5temps` applies automatic bias correction to era5
#' temperature data data to correct for unaturally low diurnal temperature
#' fluctuations in coastal grid cells.
#' @param era5hourly a stacked SpatRast of hourly ERA5 temperature data for any part
#' of the UK.
#' @param era5correctmodels a list of model correction coefficients for neach UK
#' grid cell. Available: https://universityofexeteruk-my.sharepoint.com/:f:/g/personal/i_m_d_maclean_exeter_ac_uk/EjxJdJq0MClFqUU3LATlcnEBo1cGFiUAxqLQALNNxvdZaw?e=wLR2Rf
#' Zenodo link to be added dreckly.
#' @return a SpatRast of bias corrected temperature data
#' @details `era5correctmodels` was derived by applying [biascorrect()] to a 2018
#' dataset of era5 temperature data, calibrating against Met office data
#' @import terra
#' @import mgcv
#' @export
#' @keywords biascorrect
#' @rdname correct_era5temps
correct_era5temps<-function(era5hourly,era5correctmodels) {
  # Check whether in Kelvin or degrees
  me<-mean(.is(era5hourly[[1]]),na.rm=TRUE)
  if (me > 150) era5hourly<-era5hourly-273.15
  # Calculate tmx, tmn and dtr
  tmx_era5<-hourtodayCpp(.is(era5hourly), "max")
  tmn_era5<-hourtodayCpp(.is(era5hourly), "min")
  dtr_era5<-tmx_era5-tmn_era5
  # Match cells in era5hourly to model list index
  # get x & ys of model list index
  rindex<-rast(era5correctmodels$indexr)
  xy1<-data.frame(xyFromCell(rindex, 1:ncell(rindex)))
  z<-as.vector(extract(rindex,xy1)[,2])
  xyz1 <- cbind(xy1,z)
  s<-which(is.na(z)==F)
  xyz1<-xyz1[s,]
  xyz1$x<-round(xyz1$x*4,0)/4
  xyz1$y<-round(xyz1$y*4,0)/4
  # get x & ys of model era5hourly grid cells
  xy2<-data.frame(xyFromCell(era5hourly[[1]], 1:ncell(era5hourly[[1]])))
  v<-c(1:(dim(xy2)[1]))
  xyz2 <- cbind(xy2,v)
  xyz2$x<-round(xyz2$x*4,0)/4
  xyz2$y<-round(xyz2$y*4,0)/4
  vindex<-rast(xyz2)
  mindex<-.is(vindex)
  # merge datasets
  xy<-merge(xyz1,xyz2,by=c("x","y"),all=TRUE)
  adtr<-array(NA,dim=dim(dtr_era5))
  atmn<-adtr
  for (i in 1:dim(adtr)[1]) {
    for (j in 1:dim(adtr)[2]) {
      v1<-dtr_era5[i,j,]
      v2<-tmn_era5[i,j,]
      # get model
      v<-mindex[i,j]
      s<-which(xy$v==v)
      index<-xy$z[s]
      if (is.na(index) == FALSE) {
        mod1<-era5correctmodels$dtrmod[[index]]
        mod2<-era5correctmodels$tmnmod[[index]]
        adtr[i,j,] <- predict.gam(mod1, newdata = data.frame(v2 = v1))
        atmn[i,j,] <- predict.gam(mod2, newdata = data.frame(v2 = v2))
      }
    }
  }
  # Apply dtr and min correction to original datasets
  # ~~ Calculate hourly dtr fractions
  dtr_era5h<-.ehr(dtr_era5)
  tmn_era5h<-.ehr(tmn_era5)
  tc<-.is(era5hourly)
  tfrac<-(tc-tmn_era5h)/dtr_era5h
  # Apply correct data
  adtrh<-.ehr(adtr)
  atmnh<-.ehr(atmn)
  tcn<-(tfrac*adtrh)+atmnh
  tcn<-.rast(tcn,era5hourly)
  return(tcn)
}

# ================================================================= #
# ~~~~~~~~~~~~~ Code dump from here - probably not all needed ~~~~~~~~~ #
# ================================================================= #
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
