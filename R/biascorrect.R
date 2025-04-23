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
biascorrect <- function(hist_obs, hist_mod, fut_mod = NA, mod_out = FALSE, rangelims = NA, samplenum=1000) {
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
  if (n > samplenum) s<-sample(0:n,samplenum-2,replace = FALSE)  # two less so min and max can be tagged on
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
        if (n > samplenum) {
          v1<-c(v1[1],v1[s],v1[length(v1)]) # tags on min and max value
          v2<-c(v2[1],v2[s],v2[length(v2)]) # tags on min and max value
        }
        # Apply gam
        m1 <- gam(v1~s(v2))
        if (mod_out) {
          m[i,j]<-counter
          mods[[counter]]<-m1
          counter<-counter+1
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
precipcorrect <- function(hist_obs, hist_mod, fut_mod=NA, mod_out = FALSE, rangelim = NA) {
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
  # Calculate modelled rainfall total and rain day frac
  rcount<-hist_mod
  rcount[rcount > 0] <-1
  rtot2<-apply(.is(hist_mod),c(1,2),sum)
  tfrac2<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate ratios obs/modelled
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
#' @param biasmods list containing two wrapped SpatRasters derived using [precipcorrect()]
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
  fut_mod<-mask(fut_mod,mu_tot)
  a<-as.array(fut_mod)
  mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
  rtot<-.rast(rtot,mu_tot)
  tfrac<-.rast(tfrac,mu_tot)
  rtot<-as.vector(t(rtot))
  rfrac<-as.vector(t(tfrac))
  mm<-rainadjustm(mm,rrain,rfrac,rtot)
  a2<-array(mm,dim=dim(a))
  # convert to SpatRast
  out<-.rast(a2,fut_mod)
  return(out)
}

#' @title Converts era5 hourly data to daily to enable bias correction to be applied to ukcp data
#' @param filein - filename (including path ) of era5 nc file with required variables
#' @param pathout - directory in which to save data
#' @param landsea - a landsea raster object of land fractions that must match the extent (or less) and crs of the era5 data
#'  era5 data will be cropped to lsm extent
#' @param elev - a spatRaster of surface elevations to correct surface pressure to sea level pressure
#' @param output - character of value tifs, singlenc or none - determines type of files (if any) written
#' @param vars - variables for which daily values calculated
#' @return list of daily spatRasters of each variable requested. Also write .tiff or .nc files if requested
#' The function saves individual terra SpatRaster objects to disk (one for each variable)
#' @keywords internal
#' @examples
#'  \dontrun{
#' vars<-c('psl','hurs','u10','v10','swrad','difrad','lwrad')
#' daily<-.era5todaily(filein=file.path(dir_era5,f),tempdir(),landsea=lsm,output='none',vars)
#' }
era5todaily<-function(filein,pathout,landsea, elev, output=c('tifs','singlenc','none'),
                       vars=c('tasmax','tasmin','sp','psl','hurs','huss','u10','v10','swrad','difrad','lwrad','skyem','pr')) {
  landsea[landsea==0]<-NA # used to mask and crop
  suppressWarnings(dir.create(pathout))

  # Find out what vars in era5 data and get rast template
  nc<-nc_open(filein)
  era5vars<-names(nc$var)[which(names(nc$var) %in% c("number","expver")==FALSE)]
  validvar<-era5vars[length(era5vars)]
  r<-rast(filein,subds=validvar)
  te<-crop(r[[1]],landsea) # use as template
  if("expver" %in% names(nc$var)) tme<-as.POSIXlt(nc$var$expver$dim[[1]]$vals,tz='GMT',origin="1970-01-01") else tme<-time(r)
  nc_close(nc)

  daily_tme<-as.POSIXlt(seq(trunc(tme[1],"days"),trunc(tme[length(tme)],"days"),by=3600*24))
  daily_list<-list()
  vnames<-c()
  lngnms<-c()
  uts<-c()

  # temperature - also required for other spec humidity
  if('tasmax' %in% vars || 'tasmin' %in% vars){
    if('t2m'%in% era5vars ==FALSE) stop('Temperature required in era5 data file!!!')
    a<-.nctoarray(filein,"t2m",landsea)
    a<-.applymask(a,landsea)
    tmx<-.hourtoday(a,max)
    tmn<-.hourtoday(a,min)
    dtr<-tmx-tmn
    # apply coastal correction
    dtr<-.coastalcorrect(dtr,landsea)
    tmean<-(tmx+tmn)/2
    tmx<-(tmean+0.5*dtr)-273.15
    tmn<-(tmean-0.5*dtr)-273.15
    # Save tmx and tmn
    vnames<-c(vnames,'tasmax','tasmin')
    lngnms<-c(lngnms,'Maximum temperature','Minimum temperature')
    uts<-c(uts,'degC','degC')
    if(output=='tifs') .saverast(tmx,te,pathout,"tasmax") else daily_list[["tasmax"]]<-.rast(tmx,te)
    if(output=='tifs') .saverast(tmn,te,pathout,"tasmin") else daily_list[["tasmin"]]<-.rast(tmn,te)
  }

  # Output surface pressure
  if('sp' %in% vars){
    if('sp'%in% era5vars ==FALSE) stop('Surface pressure required in era5 data file!!!')
    a2<-.nctoarray(filein,"sp",landsea)/1000
    a2<-.applymask(a2,landsea)
    sp<-.hourtoday(a2)
    vnames<-c(vnames,'sp')
    lngnms<-c(lngnms,'Surface pressure')
    uts<-c(uts,'kPa')
    if(output=='tifs') .saverast(sp,te,pathout,"sp") else daily_list[["sp"]]<-.rast(sp,te)
  }

  # Output sea level pressure
  if('psl' %in% vars){
    if('msl' %in% era5vars){ # sea level pressure in inputs variables
      a3<-.nctoarray(filein,"msl",landsea)/1000
      a3<-.applymask(a3,landsea)
      psl<-.hourtoday(a3)
    }
    if('msl' %in% era5vars ==FALSE){ #  no sea level pressure in inputs calculate from surface pressure
      if('sp'%in% era5vars ==FALSE) stop('Surface pressure required in era5 data file!!!')
      a2<-.nctoarray(filein,"sp",landsea)/1000
      a2<-.applymask(a2,landsea)

      dtm<-as.matrix(elev,wide=TRUE)
      psl<-.is(a2)/(((293-0.0065*.rta(dtm,dim(a2)[3]))/293)^5.26)
      psl<-.hourtoday(psl)
    }
    vnames<-c(vnames,'psl')
    lngnms<-c(lngnms,'Sea surface pressure')
    uts<-c(uts,'kPa')
    if(output=='tifs') .saverast(psl,te,pathout,"psl") else daily_list[["psl"]]<-.rast(psl,te)
  }

  # Humidity - specific or relative from dewpoint temp and temperature
  if('huss' %in% vars || 'hurs' %in% vars){
    if('d2m'%in% era5vars ==FALSE) stop('Dewpoint temperature required in era5 data file!!!')
    if(!exists("a")){
      if('t2m'%in% era5vars ==FALSE) stop('Temperature required in era5 data file!!!')
      a<-.nctoarray(filein,"t2m",landsea)
      a<-.applymask(a,landsea)
    }

    a3<-.nctoarray(filein,"d2m",landsea)
    a3<-.applymask(a3,landsea)
    ea<-.satvap(a3-273.15)
    es<-.satvap(a-273.15)
    rh<-(ea/es)*100

    # Relative
    if('hurs' %in% vars){
      hurs<-.hourtoday(rh)
      vnames<-c(vnames,'hurs')
      lngnms<-c(lngnms,'Relative humidity')
      uts<-c(uts,'%')
      if(output=='tifs') .saverast(hurs,te,pathout,"hurs") else daily_list[["hurs"]]<-.rast(hurs,te)
    }

    # Specific humidity requires temp and surface pressure to convert from relative
    if('huss' %in% vars){
      if(!exists("a2")){
        if('sp'%in% era5vars ==FALSE) stop('Surface pressure required in era5 data file!!!')
        a2<-.nctoarray(filein,"sp",landsea)/1000
        a2<-.applymask(a2,landsea)
      }
      a3<-suppressWarnings(converthumidity(rh, intype = "relative", outtype = "specific",
                                           tc = a-273.15, pk = a2))
      huss<-.hourtoday(a3)
      vnames<-c(vnames,'huss')
      lngnms<-c(lngnms,'Specific humidity')
      uts<-c(uts,'1')
      if(output=='tifs') .saverast(huss,te,pathout,"huss") else daily_list[["huss"]]<-.rast(huss,te)
    }
  }

  # Wind u10 NOT converted to 2m height
  if('u10' %in% vars){
    if('u10'%in% era5vars ==FALSE) stop('Wind u10 required in era5 data file!!!')
    a2<-.nctoarray(filein,"u10",landsea)
    a2<-.applymask(a2,landsea)
    uas<-.hourtoday(a2)
    vnames<-c(vnames,'u10')
    lngnms<-c(lngnms,'10 metre U wind component')
    uts<-c(uts,'m/s')
    if(output=='tifs') .saverast(uas,te,pathout,"uas") else daily_list[["u10"]]<-.rast(uas,te)
  }
  # Wind v10 NOT converted to 2m height
  if('v10' %in% vars){
    if('v10'%in% era5vars ==FALSE) stop('Wind v10 required in era5 data file!!!')
    a2<-.nctoarray(filein,"v10",landsea)
    a2<-.applymask(a2,landsea)
    vas<-.hourtoday(a2)
    vnames<-c(vnames,'vas')
    lngnms<-c(lngnms,'10 metre V wind component')
    uts<-c(uts,'m/s')
    if(output=='tifs') .saverast(vas,te,pathout,"vas") else daily_list[["v10"]]<-.rast(vas,te)
  }

  # swdown - total downward shortwave from total or mean radiation
  if('swrad' %in% vars){
    if(any(c('ssrd' ,'msdwswrf') %in% era5files[1:8])) stop('Total or mean flux SW downward required in era5 data file!!!')
    if('ssrd' %in% era5vars) a2<-.nctoarray(filein,"ssrd",landsea) else .nctoarray(filein,"msdwswrf",landsea)
    a2<-.applymask(a2,landsea)
    if('ssrd' %in% era5vars) rss<-.hourtoday(a2)/3600 else rss<-.hourtoday(a2) # convert to Wm2 if required
    vnames<-c(vnames,'swrad')
    lngnms<-c(lngnms,'Total downward shortwave radiation')
    uts<-c(uts,'watt/m^2')
    if(output=='tifs') .saverast(rss,te,pathout,"swrad") else daily_list[["swrad"]]<-.rast(rss,te)
  }
  # Downward diffuse radiation from ssrd and fdir in Jm-2 or from mean fluxes
  if('difrad' %in% vars){
    if((all(c('ssrd','fdir') %in% era5vars) || all(c('msdwswrf','msdrswrf') %in% era5vars))==FALSE) stop('Total and direct SW radiation (as flux or total) required in era5 data file!!!')
    if(all(c('ssrd','fdir') %in% era5vars)) a1<-.nctoarray(filein,"ssrd",landsea) else a1<-.nctoarray(filein,"msdwswrf",landsea)
    if(all(c('ssrd','fdir') %in% era5vars)) a2<-.nctoarray(filein,"fdir",landsea) else a1<-.nctoarray(filein,"msdrswrf",landsea)
    a3<-a1-a2
    a3<-.applymask(a3,landsea)
    if(all(c('ssrd','fdir') %in% era5vars)) difsw<-.hourtoday(a3)/3600 else difsw<-.hourtoday(a3)
    vnames<-c(vnames,'difrad')
    lngnms<-c(lngnms,'Downward diffuse radiation')
    uts<-c(uts,'watt/m^2')
    if(output=='tifs') .saverast(difsw,te,pathout,"difrad") else daily_list[["difrad"]]<-.rast(difsw,te)
  }

  # Downward LW radiation from msdwlwrf in Wm-2
  if('lwrad' %in% vars){
    if('msdwlwrf' %in% era5vars ==FALSE) stop('Mean downward LW flux required in era5 data file!!!')
    a2<-.nctoarray(filein,"msdwlwrf",landsea)
    a2<-.applymask(a2,landsea)
    lwrad<-.hourtoday(a2)
    vnames<-c(vnames,'lwrad')
    lngnms<-c(lngnms,'Total downward longwave radiation')
    uts<-c(uts,'watt/m^2')
    if(output=='tifs') .saverast(lwrad,te,pathout,"lwrad") else daily_list[["lwrad"]]<-.rast(lwrad,te)
  }

  # skyem
  if('skyem' %in% vars){
    if(all(c('msdwlwrf','msnlwrf') %in% era5vars)==FALSE) stop('Mean downward and net LW fluxes required in era5 data file!!!')
    lwdn<-.nctoarray(filein,"msdwlwrf",landsea)
    lwme<-.nctoarray(filein,"msnlwrf",landsea)
    lwup<-(-lwme+lwdn)
    skyem<-lwdn/lwup
    skyem<-.applymask(skyem,landsea)
    skyem[skyem>1]<-1
    skyem<-.hourtoday(skyem)
    vnames<-c(vnames,'skyem')
    lngnms<-c(lngnms,'Sky emmisivity')
    uts<-c(uts,'watt/m^2')
    if(output=='tifs') .saverast(skyem,te,pathout,"skyem") else daily_list[["skyem"]]<-.rast(skyem,te)
  }

  # pr
  if('pr' %in% vars){
    if('tp' %in% era5vars ==FALSE) stop('Total precipitation required in era5 data file!!!')
    a<-.nctoarray(filein,"tp",landsea)
    a[is.na(a)]<-0
    pr<-.hourtoday(a,sum)*1000
    pr<-.applymask(pr,landsea)
    vnames<-c(vnames,'pr')
    lngnms<-c(lngnms,'Precipitation')
    uts<-c(uts,'mm/day')
    if(output=='tifs') .saverast(pr,te,pathout,"pr") else daily_list[["pr"]]<-.rast(pr,te)
  }

  # Add metadata to list
  for(n in 1:length(daily_list)){
    time(daily_list[[n]])<-daily_tme
    units(daily_list[[n]])<-uts[n]
    longnames(daily_list[[n]])<-lngnms[n]
  }

  # Write nc of all variables if requested
  if(output=='singlenc'){
    dat<-terra::sds(daily_list)
    names(dat)<-vnames
    longnames(dat)<-lngnms
    units(dat)<-uts
    filename<-paste0(tools::file_path_sans_ext(basename(filein)),"_daily.nc")
    writeCDF(dat,file.path(pathout,filename),overwrite=TRUE,prec='double',compression=9)
  }
  return(daily_list)
}


# ============================================================================ #
# ~~~~~~~~~~~~~ Bias-correct worker functions here - CHECK REQUIRESD ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #

#' @title Takes global ukcp data as inputs, ensures correct number of dates in each year,
#' interpolating missing values, resamples data to match era5 and returns data for a whole true decade
#' @param ukcpfile1 - filename of first file coinciding wiht decade - e.g. for decade 2010-2020,
#' @param data with an extension 20091201-20191130.nc
#' @param ukcpfile2 - filename of second file coinciding wiht decade - e.g. for decade 2010-2020,
#' data with an extension 20191201-20291130.nc
#' @param landsea - a landsea raster object of land fractions that must match the extent of the era5 data.
#' must also match the resolution of era5 data as it is used for resampling ukcp data
#' @param decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
#' returns a SpatRaster object for an entire decade
#' @noRd
.cropandsortUKCPone <- function(ukcpfile1,ukcpfile2,landsea,decade=1) {
  # Get and crop nc files
  ecrop<-extent(landsea)
  r1<-.cropnc(ukcpfile1,ecrop)
  r2<-.cropnc(ukcpfile2,ecrop)
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
  a<-abind::abind(a1,a2)
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
  ao<-apply(ao,c(1,2),zoo::na.approx)
  ao<-aperm(ao,c(2,3,1))
  # resample and mask out coastal area
  ro<-rast(ao)
  ext(ro)<-ext(r1)
  rte<-rast(landsea)
  crs(ro)<-crs(rte)
  ro<-resample(ro,rte)
  a<-as.array(ro)
  landsea[landsea==0]<-NA
  a<-.applymask(a,landsea)
  ro2<-rast(a)
  ext(ro2)<-ext(ro)
  crs(ro2)<-crs(rte)
  return(ro2)
}
#' @title For one tile an one model run and one decade applies function cropandsortUKCPone to all variables
#' @param pathtoUKCP - directory with raw UKCP data for that model run and tile
#' @param pathout - directory in which to save data
#' @param landsea - a landsea raster object of land fractions that must match the extent of the era5 data
#' @param decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
#' @param modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
#' data, so file naming and foldr convention assumed to match that of data supplied via dropbox
#' @return saves SpatRaster objects for each variable to disk in directory pathout
#' @noRd
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
#' @title For one tile an one model run extracts the 2018 data and saves this to disk
#' @param pathin - directory with data - should match pathour of cropandsortUKCPdecade function
#' @param pathout - directory in which to save data - usually a subdirectory of pathin
#' @param modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
#' data, so file naming and foldr convention assumed to match that of data supplied via dropbox
#' @return saves 2018 data to disk as terra SpatRaster objects
#' @noRd
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
    writeRaster(ro,filename=fo,overwrite=TRUE)
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

# ================================================================= #
# ~~~~~~~~~~~~~ Code dump from here - probably not all needed ~~~~~~~~~ #
# ================================================================= #
#' @title Saves a SpatRaster object as an nc file
#' @param r - a SpatRaster object
#' @param baseyear - used for adding time stamp to nc file (asumes daily)
#' @param shortname - Name of the dataset
#' @param longname - Long name of the dataset
#' @param unit - character. Unit of the data
#' @param fileout - name of file to be saved (includign path)
#' saves a compressed netCDF4 file to disk
#' @noRd
savenc<-function(r,baseyear,shortname,longname,unit,fileout) {
  terra::time(r)<-as.Date(paste0(baseyear,"-01-01")) + 0:(nlyr(r)-1)
  writeCDF(r,fileout,overwrite=TRUE,compression=9,varname=shortname,
           longname=longname,unit=unit)
}



