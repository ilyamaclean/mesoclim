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
  if (!inherits(fut_mod, "logical") & crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
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
    ao<-crop(ao,fut_mod)
  }
  return(ao)
}
#' @title Applies bias correction model to climate variable
#' @description The function `biascorrect_apply` applies a set of bias corrections
#' derived using `bias_correct` to a spatial climate dataset
#' @param fut_mod a SpatRaster of modelled climate data for for e.g. a future
#' period to which corrections are applied.
#' @param biasmods an object of class `biascorrectmodels` derived using [biascorrect()].
#' @return a SpatRaster matching the overlap of `fut_mod` and `biasmods$rindex`
#' @details The `biascorrect_apply` function allows the same bias correction to be applied
#' to multiple datasets, spanning e.g. several future time periods or climate scenarious.
#' The spatial extent of `fut_mod` must at least partially overlap with that of area used in
#' `biasmods` when running `biascorrect` as corrected values are only returned for the overlapping area
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
  if (!inherits(fut_mod, "SpatRaster")) stop("In biascorrect_apply the parameter fut_mod must be a SpatRaster")
  if (!inherits(biasmods, "biascorrectmodels")) stop("biasmods must be an object of class biascorrectmodels")
  if (inherits(biasmods$rindex,"PackedSpatRaster")) r<-unwrap(biasmods$rindex) else r<-biasmods$rindex
  mods<-biasmods$models
  if (crs(fut_mod) != crs(r)) fut_mod<-project(fut_mod,crs(r))
  if(!relate(ext(r), ext(fut_mod),"intersects")) warning("Data and bias correction model in biascorrect_apply do not overlap!!!")
  if(!relate(ext(fut_mod),ext(r),"coveredby")) warning("Data in biascorrect_apply exceeds model area - outputs will be cropped!!!")
  te<-fut_mod
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
  ao<-crop(ao,te)
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
#' @param prec_thold - numeric value (mm per day) below which daily precipitation is taken to be zero.
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
precipcorrect <- function(hist_obs, hist_mod, fut_mod=NA, mod_out = FALSE, rangelim = NA, prec_thold=0.01) {
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
  rcount[rcount >= prec_thold] <-1
  rtot1<-apply(.is(hist_obs),c(1,2),sum)
  tfrac1<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate modelled rainfall total and rain day frac
  rcount<-hist_mod
  rcount[rcount >= prec_thold] <-1
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
    out<-crop(out,fut_mod)
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
#' which is used as a basis for adjusting individual rainfall events.
#' @param prec_thold numeric value (mm per day) below which precipitation is taken to be zero
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
precipcorrect_apply<-function(fut_mod, biasmods, prec_thold=0.01) {
  if (!inherits(fut_mod, "SpatRaster")) stop("fut_mod must be a SpatRaster")
  if (inherits(biasmods$mu_tot,"PackedSpatRaster")) mu_tot<-unwrap(biasmods$mu_tot) else mu_tot<-biasmods$mu_tot
  if (inherits(biasmods$mu_frac,"PackedSpatRaster")) mu_frac<-unwrap(biasmods$mu_frac) else mu_frac<-biasmods$mu_frac
  # reproject and crop if necessary
  if (crs(fut_mod) != crs(mu_tot)) fut_mod<-project(fut_mod,crs(mu_tot))
  if(!relate(ext(mu_tot), ext(fut_mod),"intersects")) warning("Data and bias correction model in precipcorrect_apply do not overlap!!!")
  if(!relate(ext(fut_mod),ext(mu_tot),"coveredby")) warning("Data in precipcorrect_apply exceeds model area - outputs will be cropped!!!")
  te<-fut_mod[[1]]
  fut_mod<-resample(fut_mod,mu_tot)
  # Calculate and adjust rainfall total and rain day frac
  rcount<-fut_mod
  rcount[rcount >= prec_thold] <-1
  rtot<-apply(.is(fut_mod),c(1,2),sum,na.rm=TRUE)*.is(mu_tot)
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
  # convert to SpatRast and crop to input fut_mod
  out<-.rast(a2,fut_mod)
  out<-crop(out,te)
  return(out)
}

#' Apply bias correction models to all climate variables
#'
#' @param climdata - a list of climate and associated data as output by `ukcp18toclimarray` function
#' @param model_list - a list of bias correction models as output by `biascorrect` and `precipcorrect` functions
#' names of list elements must match variable names in `climdata`.
#' @param rangelim fractional amount above the maximum values of precipitation
#' in model observed data in any given cell by which values in the returned dataset can
#' attain. I.e. if set to 2, values in any given grid cell of the returned dataset
#' cannot exceed twice the maximum amount in the corresponding grid cell in
#' observed data of model.If not supplied (the default), returned values are unbounded.
#' @param prec_thold - numeric value (mm per day) below which daily precipitation is taken to be zero.
#'
#' @return list of bias corected spatRasters and other variables matching `climdata` input
#' @details
#' SpatRaster output may be cropped if climdata input extent exceeds the extent over which bias correction
#' models were developed.
#'
#' @export
#'
#' @examples
biascorrect_climdata<-function(climdata, model_list, prec_thold=0.01, rangelims = NA, fillna=FALSE){
  # Check inputs - names match and convert to spatrasters if required
  vars<-c("relhum","pres","lwrad","swrad","tmax","tmin","windspeed","prec" )
  input_vars <-names(climdata)
  model_vars<-names(model_list)
  if(!all(vars %in% input_vars)) stop("Cannot find all input variable names in climdata parameter!!!")
  if(!all(vars %in% model_vars)) stop("Cannot find all input variable names in model_list!!!")
  input_class<-lapply(lapply(climdata,class),`[[`, 1)
  if(any(input_class=="PackedSpatRaster")) climdata[which(input_class=="PackedSpatRaster")]<-lapply(climdata[which(input_class=="PackedSpatRaster")],unwrap)
  if(any(input_class=="array")) climdata[which(input_class=="array")]<-lapply(climdata[which(input_class=="array")],.rast,tem=climdata$dtm)

  # Bias correct variables
  for (v in vars){
    if(v!="prec") climdata[[v]]<-biascorrect_apply(fut_mod=climdata[[v]], biasmods=model_list[[v]], rangelims)
    if(v=="prec") climdata[[v]]<-precipcorrect_apply(fut_mod=climdata[[v]], biasmods=model_list[[v]], prec_thold)
    if(fillna) climdata[[v]]<-.spatinterp(climdata[[v]])
  }
  ## Correct other variables to match output of bias correction (may crop)
  if(ext(climdata$dtm)!=ext(climdata[[v]])){
    r<-ext(climdata[[v]])
    climdata$winddir<-crop(climdata$winddir,r)
    if("cloud" %in% input_vars) climdata$cloud<-crop(climdata$cloud,r)
    climdata$dtm<-crop(climdata$dtm,r)
  }
  return(climdata)
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



