
#' @title Get HadUK data from ceda archive on jasmine
#'
#' @param startdate - POSIXlt value indicating start date for data required
#' @param enddate - POSIXlt value indicating end date for data required
#' @param dtmc - spatraster to which data will be resampled/cropped
#' @param filepath - dir where HADUK source files located
#' @param varn  variable required, one or more of: 'rainfall','tasmax','tasmin' available at daily time steps
#'
#' @return spatRaster
#' @export
#' @keywords jasmin
#' @seealso [download_hadukdaily()]
#' \dontrun{
#'
#' }
addtrees_hadukdata<-function(startdate, enddate, dtmc, filepath="/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda/1km/rainfall/day/latest",
                             var=c('rainfall','tasmax','tasmin')) {
  # Checks
  varn<-match.arg(var)
  if (any(!var %in% c("rainfall","tasmax","tasmin")) ) stop("Chosen variables are not available as daily data!!" )
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")

  # Derive months of data required
  dateseq<-seq(as.Date(startdate) + 1, as.Date(enddate) + 1, by = "1 months") - 1
  yrs<-year(dateseq)
  mnths<-month(dateseq)

  # Get date text used in file names
  mtxt<-ifelse(mnths<10,paste0("0",mnths),paste0("",mnths))
  daysofmonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdays<-daysofmonth[mnths]
  mdays<-ifelse((yrs%%4==0 & mnths==2),29,mdays)
  files<-paste0(var,"_hadukgrid_uk_1km_day_",yrs,mtxt,"01-",yrs,mtxt,mdays,".nc")

  # Get aoi in projection of HadUK data OS projection and extend by one cell
  r<-rast(file.path(filepath,files[1]))
  if(!compareCRS(dtmc,r)) dtmc<-terra::project(dtmc,r)
  dtmc_ext<-crop(extend(dtmc,1),dtmc)

  # Load all monthly files to spatRaster then resample
  r_list<-list()
  n<-1
  for (f in files){
    r<-rast(file.path(filepath,f))
    #terra::crs(r)<-'epsg:27700'
    e<-terra::crop(terra::extend(terra::crop(r[[1]],dtmc,snap='out'),1),r[[1]])
    r<-terra::crop(r,e)
    r<-.resample(r,dtmc,msk=TRUE)
    r_list[[n]]<-r
    n<-n+1
  }
  rout<-rast(r_list)
  return(rout)
}


#' @title Get ukcp18 RCM 12km DTM from ceda archive
#'
#' @param aoi - used to cut dtm
#' @param basepath - "" when using on jasmin
#'
#' @return spatrast object cropped to extended aoi
#' @export
#' @import terra
#' @import units
#' @keywords jasmin
#' @examples
#' \dontrun{
#' dtmc<-get_ukcp_dtm(aoi, basepath=ceda_basepath)
#' }
get_ukcp_dtm<-function(aoi, basepath=""){
  # Load dtm
  fname<-"orog_land-rcm_uk_12km_osgb.nc"
  path<-file.path(basepath,"badc/ukcp18/data/land-rcm/ancil/orog")
  dtmc<-terra::rast(file.path(path,fname))
  #terra::crs(dtmc)<-'epsg:27700'
  # Convert aoi to dtmc projection and if a vector convert to bounding box vect
  aoiproj<-terra::project(aoi,terra::crs(dtmc))
  if(class(aoi)[1]=='SpatVector'){
    aoiproj<-terra::vect(terra::ext(aoiproj))
    terra::crs(aoiproj)<-terra::crs(dtmc)
  }
  # Crop dtmc to AOI then add a cell to all surrounds (if within original dtmc)
  dtmc_ext<-terra::crop(terra::extend(terra::crop(dtmc,aoiproj,snap='out'),1),dtmc)
  dtmc_out<-terra::crop(dtmc,dtmc_ext)
  units(dtmc)<-'m'
  names(dtmc)<-'Elevation'
  return(dtmc_out)
}


#' @title Calculate medium resolution dtm from aoi, coarse (ukcp) dtmc and uk wide fine scale dtmuk
#'
#' @param aoi - area of interest spatvector, raster or sf
#' @param dtmc - coarse dtm spatraster
#' @param dtmuk - uk wide fine resolution dtm spatraster
#' @param basepath - "" when using on jasmin
#'
#' @return spatraster object
#' @export
#' @import terra
#' @keywords preprocess
#' @examples
#' \dontrun{
#' dtmm<-get_dtmm(aoi,dtmc,dtmuk,basepath=ceda_basepath)
#' }
get_dtmm<-function(aoi,dtmc,dtmuk,basepath=""){
  dtmm_res<-round(exp( ( log(terra::res(dtmc)[1]) + log(terra::res(dtmf)[1]) ) / 2 ))
  dtmm_f<-terra::mask(terra::crop(terra::crop(dtmuk,dtmc),dtmuk),coast_v)
  dtmm<-terra::mask(terra::aggregate(dtmm_f,dtmm_res / res(dtmf),  na.rm=TRUE),coast_v, touches=FALSE)
  return(dtmm)
}

#' @title Source UKCP files from ceda archive and preprocess
#'
#' @param aoi - area of interest as spatvector, raster or sf
#' @param startdate POSIXlt class defining start date of required timeseries
#' @param enddate POSIXlt class defining end date of required timeseries
#' @param collection text string defining UKCP18 collection, either 'land-gcm' or 'land-rcm'
#' @param domain text string defining UKCP18 domain, either 'uk' or 'eur'(land-rcm collection only) or 'global'
#' @param member string defining the climate model member to be used for the timeseries. Available members vary between UKCP18 collections.
#' @param basepath - "" for jasmin use
#' @param wsalbedo - white sky albedo value
#' @param bsalbedo - black sky albedo value
#' @param ukcp_vars UKCP18 variable names to be extracted DO NOT CHANGE
#' @param ukcp_units units of the UKCP18 variables extracted DO NOT CHANGE
#' @param output_units units required for output DO NOT CHANGE CHECK CORRECT
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasts.
#' @param sampleplot if TRUE plots examples of interpolated dates when converting from 360 to 366 day years
#'
#' @return a list of the following:
#' \describe{
#'    \item{dtm}{Digital elevation of downscaled area in metres (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2)}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed}{At `windheight_m` above ground` (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @import terra
#' @import units
#' @keywords jasmin
#' @seealso [ukcp18toclimarray()]
#' @examples
#' \dontrun{
#' climdata<-addtrees_climdata(aoi,ftr_sdate,ftr_edate,collection='land-rcm',domain='uk',member='01',basepath=ceda_basepath)
#' }
addtrees_climdata <- function(aoi,  startdate, enddate,
                              collection=c('land-gcm','land-rcm'),
                              domain=c('uk','eur','global'),
                              member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                                       '16','17','18','19','20','21','22','23','24','25','26','27','28'),
                              basepath="",
                              wsalbedo=0.19, bsalbedo=0.22,
                              ukcp_vars=c('clt','hurs','pr','psl','rls','rss',
                                          'tasmax','tasmin','uas','vas'),
                              ukcp_units=c('%','%','mm/day','hPa','watt/m^2','watt/m^2',
                                           'degC','degC','m/s','m/s'),
                              output_units=c('%','%','mm/day','kPa','watt/m^2','watt/m^2',
                                             'degC','degC','m/s','m/s'),
                              toArrays=TRUE, sampleplot=FALSE){
  # Parameter check
  collection<-match.arg(collection)
  domain<-match.arg(domain)
  ukcp_vars<-match.arg(ukcp_vars,several.ok=TRUE)
  ukcp_units<-match.arg(ukcp_units,several.ok=TRUE)
  output_units<-match.arg(output_units,several.ok=TRUE)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  if(!class(aoi)[1] %in% c("SpatRaster","SpatVector","sf") ) stop("aoi parameter is NOT a SpatRaster, SpatVector or sf object!!!")
  if(class(aoi)[1] =="sf") aoi<-terra::vect(aoi)


  # Check member in chosen collection
  member<-match.arg(member)
  if(collection =='land-rcm' &
     !member %in% c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15')) stop(paste("Model member",member,"NOT available in land-rcm collection - ignoring!!"))

  # Add rcp and collection resolution
  rcp<-'rcp85'
  if(collection=='land-gcm') collres<-'60km' else collres<-'12km'

  # Identify which decades are required
  decades<-.find_ukcp_decade(collection,startdate,enddate)

  # Create coarse-resolution dtm to use as template for cropping etc !!!
  dtmc<-get_ukcp_dtm(aoi,basepath)

  # Jasmin basepath
  basepath<-file.path(basepath,"badc","ukcp18","data",collection,domain,collres,rcp)

  # Load spatrasters from ukcp.nc file, crop to dtmc, convert to normal calendar and add to output list
  clim_list<-list()
  for (n in 1:length(ukcp_vars)){
    v<-ukcp_vars[n]
    filepath<-file.path(basepath,member,v,'day','latest')
    ukcp_u<-ukcp_units[n]
    out_u<-output_units[n]
    var_r<-terra::rast()
    for (d in decades){
      filename<-paste0(v,'_rcp85_',collection,'_',domain,'_',collres,'_',member,'_day_', d,".nc")
      ncfile<-file.path(filepath, filename)
      message(paste("Loading",ncfile))

      # Load and crop data in native crs then project to dtmc crs and recrop
      r <- terra::rast(ncfile, subds=v)
      terra::crs(r)<-terra::crs(dtmc)
      r<-terra::crop(r,dtmc)

      # If requested then convert to real calendar dates and fill missing dates
      ukcp_dates<-.get_ukcp18_dates(ncfile)
      real_dates<-.correct_ukcp_dates(ukcp_dates)
      r<-.fill_calendar_data(r, real_dates, testplot=sampleplot)

      # Correct units
      terra::units(r)<-ukcp_u

      # PROBLEM with u v wind vectors in global data - based on DIFFERENT grid to other variables!!!
      if(round(terra::ext(r),2)!=round(terra::ext(dtmc),2)){
        message("Extents of input data DIFFERENT - correcting!!")
        r<-terra::resample(r,dtmc)
      }
      # Join if multiple decades of data
      terra::add(var_r)<-r
    }
    # Check & convert units, check all rasters geoms are same
    if(ukcp_u!=out_u) var_r<-.change_rast_units(var_r, out_u)
    if(!terra::compareGeom(dtmc,var_r)) warning(paste(v,"Spatrast NOT comparable to DTM!!") )
    clim_list[[v]]<-var_r
  }

  # Calculate derived variables: wind
  clim_list$windspeed<-terra::rast(sqrt(as.array(clim_list$uas)^2+as.array(clim_list$vas)^2)*log(67.8*2-5.42)/log(67.8*10-5.42),extent=terra::ext(dtmc),crs=terra::crs(dtmc))
  clim_list$winddir<-terra::rast(as.array((terra::atan2(clim_list$uas,clim_list$vas)*180/pi+180)%%360),crs=terra::crs(dtmc)) #  (deg from N)
  units(clim_list$windspeed)<-'m/s'
  units(clim_list$winddir)<-'deg'
  terra::time(clim_list$windspeed)<-terra::time(clim_list$uas)
  terra::time(clim_list$winddir)<-terra::time(clim_list$uas)
  names(clim_list$windspeed)<-terra::time(clim_list$windspeed)
  names(clim_list$winddir)<-terra::time(clim_list$winddir)

  # Calculate derived variables: longwave downward
  tmean<-(clim_list$tasmax+clim_list$tasmin)/2
  lwup<-terra::app(tmean, fun=.lwup)
  clim_list$lwdown<-clim_list$rls+lwup

  # Calculate derived variables: shortwave downward from white & black sky albedo as rast timeseries or fixed land and sea values
  clim_list$swdown<-.swdown(clim_list$rss, clim_list$clt, dtmc, wsalbedo, bsalbedo)

  # Select and rename climate output rasts MIGHT NEED TO CHANGE THESE TO MATCH THOSE USED BY MESOCLIM FUNCTIONS
  clim_list<-clim_list[c("clt","hurs","pr","psl","lwdown","swdown","tasmax","tasmin", "windspeed","winddir")]
  names(clim_list)<-c('cloud','relhum','prec','pres','lwrad','swrad','tmax','tmin','windspeed','winddir')

  # Restrict to dates requested
  filter_times<-function(x,startdate,enddate) x[[which(time(x) >= startdate & time(x) <= enddate)]]
  for(v in names(clim_list)) clim_list[[v]]<-filter_times(clim_list[[v]],startdate,enddate)

  ### Create output list
  output_list<-list()
  output_list$dtm<-dtmc
  output_list$tme<-as.POSIXlt(time(clim_list[[1]]),tz="UTC")
  output_list$windheight_m<-10 # ukcp windspeed at 10 metres height
  output_list$tempheight_m<-1.5 # ukcp air temp at 1.5 metres height

  # Convert climate data to arrays if required
  if(toArrays) clim_list<-lapply(clim_list,as.array)

  output_list<-c(output_list,clim_list)
  return(output_list)
}


#' @title Create UKCP sea surface temperature rast stack from ceda archive
#'
#' @param startdate start date as POSIXlt
#' @param enddate end date as POSIXlt
#' @param aoi SpatRaster, SpatVector or sf object defining area of interest to crop data. NA no cropping occurs
#' @param members model members to be included
#' @param v = SST sea surface temperature
#' @param basepath - "" for jasmin use
#'
#' @return Spatraster timeseris of sea surface temperatures in original projection
#' @export
#' @import terra
#' @import lubridate
#' @keywords jasmin
#' @seealso [create_ukcpsst_data()]
#' @examples
#' \dontrun{
#' sstdata<-addtrees_sstdata(ftr_sdate,ftr_edate,aoi=climdata$dtm,member='01',basepath=ceda_basepath)
#' }
addtrees_sstdata<-function(
    startdate,
    enddate,
    aoi=NA,
    member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15') ,
    basepath="",
    v='SST' ){
  # Check parameters
  member<-match.arg(member)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  all_land<-FALSE

  # Derive months of data required - will output month before and after start and end dates for interpolation
  start<-startdate %m-% months(1)
  end<-enddate %m+% months(1)
  yrs<-unique(c(year(start):year(end)))

  # Get member ID used in file names from mesoclim lookup table
  # Ref: https://www.metoffice.gov.uk/binaries/content/assets/metofficegovuk/pdf/research/ukcp/ukcp18-guidance-data-availability-access-and-formats.pdf
  modelid<-ukcp18lookup$PP_ID[which(ukcp18lookup$Member_ID == member)]
  if("" %in% modelid) stop(paste("Model NOT available for sea surface temperature!!"))

  # Get filepaths and names of sst data required
  filepath<-file.path(basepath,'badc/deposited2023/marine-nwsclim/NWSPPE',modelid,'annual')
  ncfiles<-do.call(paste0, c(expand.grid('NWSClim_NWSPPE_',modelid,'_',yrs,'_gridT.nc') ))
  ncfiles<-file.path(filepath,ncfiles)
  not_present<-which(!file.exists(ncfiles))
  if (length(not_present)>0) stop(paste("Input .nc files required are NOT present: ",ncfiles[not_present]," ") )

  # Create generous cropping area and check not all land
  if(!class(aoi)=='logical'){
    if(!class(aoi)[1] %in% c("SpatRaster","SpatVector","sf")) stop("Parameter aoi NOT of suitable spatial class ")
    if(class(aoi)[1]=="sf") aoi<-vect(aoi)
    r<-terra::rast(ncfiles[1], subds=v)[[1]]
    # Crop to area of dtmc plus one cell on all sides
    aoi_e<-ext(project(aoi,terra::crs(r)))
    aoi_r<-crop(extend(crop(r,aoi_e,snap='out'),1),r)
    aoi_r<-crop(r,aoi_r)
    # Check if any sea cells in crop area
    if(all(is.na(values(r)))){
      warning('No sea surface temperature data in dtmc area or adjacent cells!!! Returning NA...')
      all_land<-TRUE
      var_r<-NA
    }
  }
  # Get cropped spatrast stack of sea surface temperatures
  if(!all_land){
    var_r<-terra::rast()
    for(f in ncfiles){
      r<- rast(f, subds = v, drivers="NETCDF")
      if(!class(aoi)[1]=='logical') r<-crop(r,aoi_e)
      units(r)<-'degC'
      # Join if multiple decades of data
      terra::add(var_r)<-r
    }
    # Select relevant months
    tme<-terra::time(var_r)
    var_r<-var_r[[which(terra::time(var_r) %within%  interval(start-month(1), end+month(1)) ) ]]
    names(var_r)<-terra::time(var_r)
  }
  return(var_r)
}
