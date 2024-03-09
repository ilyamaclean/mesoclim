#' @title Downloads 1km historic climata data for the UK
#' @description The function `download_hadukdaily` downloads HadUK-Grid Gridded
#' Climate Observation data on a 1km grid over the UK from CEDA.
#' @param dir_out directory to which to save data
#' @param cedausr string of ceda username
#' @param cedapwd string of ceda username
#' @param year numeric value indicating the year for which data are required
#' @param month numeric value (1-12) indicating the month for which data are required
#' @param varn. Variable required. One of `rainfall` (precipitation in mm), `tasmax`
#' (maximum daily temperature, deg C) or `tasmin` (minimum daily temperature, deg C).
#' @details To obtain a username and password, first register at https://services.ceda.ac.uk/.
#' Your username will be the same as that given for your main CEDA account.
#' Your FTP password is separate from the password for your CEDA web account. If
#' you do not already have a password, or if you want to reset your FTP password
#' go to MyCEDA and click on the 'Configure FTP account' button.
#' @import curl
#' @export
download_hadukdaily<-function(dir_out, cedausr, cedapwd, year, month, varn) {
  # Checks
  if (varn != "rainfall" & varn != "tasmax" & varn != "tasmin") stop("daily data are not available for this variable")
  tme<-as.POSIXlt(Sys.time())
  yr<-tme$year+1900
  if (year > yr) stop("Function downloads observed data. Cannot be used for future years")
  if (year == yr) warning("Data may not yet be available for current year")
  if (varn == "rainfall") {
    if (year < 1891) stop("Rainfall data not available prior to 1891")
  } else {
    if (year < 1960) stop("Temperature data not available prior to 1960")
  }
  mtxt<-ifelse(month<10,paste0("0",month),paste0("",month))
  mdays<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  if (year%%4==0) mdays[2]<-29
  mday<-mdays[month]
  # Create url
  urlbase<-paste0("ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda/1km/")
  varbase<-paste0(urlbase,varn,"/day/latest/")
  fname<-paste0(varn,"_hadukgrid_uk_1km_day_",year,mtxt,"01-",year,mtxt,mday,".nc")
  destfile<-paste0(dir_out,fname)
  dload_url<-paste0(varbase,fname)
  h <- curl::new_handle()
  curl::handle_setopt(h, userpwd = paste0(cedausr,":",cedapwd))
  curl::curl_download(dload_url, destfile, handle = h)
}
#' ILya: Have changed this to an internal function, including calls to it.
.multi_download<-function(f,url,dir_out,h){
  dload_url<-paste(url, f, sep = "")
  destfile<-file.path(dir_out,basename(f))
  curl::curl_download(dload_url, destfile, handle = h)
}
#' Download UKCP18 climate data
#' @description
#' Using parameters the function will downloaded from ftp.ceda.ac.uk all available UCKP18 files containing relevant data to the user defined output directory
#'
#' @param dir_out Output directory to which files are downloaded
#' @param cedausr string of ceda user name
#' @param cedapwd string of ceda user pasword
#' @param startdate POSIXlt class defining start date of required timeseries
#' @param enddate POSIXlt class defining end date of required timeseries
#' @param collection text string defining UKCP18 collection, either 'land-gcm' or 'land-rcm'
#' @param domain text string defining UKCP18 domain, either 'uk' or 'eur'(land-rcm collection only) or 'global'
#' @param rcp text string iof RCP scenario to be downloaded
#' @param member vector of strings defining the climate model member to be downloaded. Available members vary between UKCP18 collections.
#' @param vars vector of strings corresponding to UKCP18 short variable names to download
#' @param cedaprot string of ceda protocol set to ftp DO NOT CHANGE
#' @param cedaserv string of the ceda server url DO NOT CHANGE
#' @return downloads files to `dir_out`
#' @details To obtain a username and password, first register at https://services.ceda.ac.uk/.
#' Your username will be the same as that given for your main CEDA account.
#' Your FTP password is separate from the password for your CEDA web account. If
#' you do not already have a password, or if you want to reset your FTP password
#' go to MyCEDA and click on the 'Configure FTP account' button.
#' @import curl
#' @export
#' #' @examples
#' dir_out <- tempdir()
#' download_ukcp18(dir_out,'ceda_username','cedapwd'ceda_password',as.POSIXlt('2018-05-01'),as.POSIXlt('2018-05-31'),'land-rcm','uk','rcp85',c('01'),c('tasmax','tasmin'))
download_ukcp18<-function(
    dir_out,
    cedausr,
    cedapwd,
    startdate,
    enddate,
    collection=c('land-gcm','land-rcm'),
    domain=c('uk','eur','global'),
    rcp=c('rcp85', 'rcp26'),
    member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                '16','17','18','19','20','21','22','23','24','25','26','27','28'),
    vars=c('clt','hurs','huss','pr','prsn','psl','rls','rss','snw','tasmax','tasmin','uas','vas'), # 'tas',
    cedaprot="ftp://",
    cedaserv="ftp.ceda.ac.uk"
){
  # Check parameters
  if(!dir.exists(dir_out)) stop(paste("Output directory does not exist:",dir_out))
  collection<-match.arg(collection)
  domain<-match.arg(domain)
  rcp<-match.arg(rcp)
  member<-match.arg(member,several.ok=TRUE)
  vars<-match.arg(vars,several.ok=TRUE)
  if(collection=='land-rcm' & !member %in% c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15')){
    warning('Invalid member for land-rcm - retricting to only those valid!!')
    member<-member[which(modelsruns %in% c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'))]
  }
  if(collection=='land-gcm' & !domain %in% c('uk','global')){
    warning('Invalid area for global data - downloading global data!!')
    domain<-'global'
  }
  if(collection=='land-rcm' & !domain %in% c('uk','eur')){
    warning('Invalid area for regional data - downloading Europe data!!')
    domain<-'eur'
  }

  # Identify which decades are required
  decades<-find_ukcp_decade(collection,startdate,enddate)

  # Get tiled data according to resolution
  if(collection=='land-rcm') res<-"12km"
  if(collection=='land-gcm') res<-"60km"

  url<-paste0(cedaprot,cedaserv)
  dir_ukcp18<-paste0("/badc/ukcp18/data/",collection,"/",domain,"/",res,"/",rcp)

  print('Parameters chosen:')
  print(collection)
  print(domain)
  print(rcp)
  print(member)
  print(vars)
  print(decades)
  print(paste('Downloading from',dir_ukcp18))

  # Download files - loop over model runs and vars - use lapply to download all decades
  n<-1
  for(run in member){
    for(v in vars){
      print(n)
      # Get pathways
      dir_dload<-file.path(dir_ukcp18,run,v,'day','latest')
      files<-paste0(v,'_',rcp,'_',collection,'_',domain,'_',res,'_',run,'_day_',decades,'.nc')
      dload_files<-file.path(dir_dload,files)
      print(dload_files)
      n<-n+1
      # Handles Multifile Download from server
      h <- curl::new_handle()
      curl::handle_setopt(h, userpwd = paste0(cedausr,":",cedapwd))
      lapply(dload_files, FUN=.multi_download, url=url,dir_out=dir_out,h=h)
    }
  }
  # TO ADD Get ancillary data if available?? ie original DEM for regional data
}

#' @title Gets dates of UKCP18 ncdf file
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by [mcera5::request_era5()] to the correct formal required for subsequent modelling.
#' @param ncfile filename of UKCP18 ncdf data file
#' @return a vector of date strings in the form 'yyy-dd-mm'
#' @import ncdf4, lubridate
#' @export
#' @details Creates a time series from the ncdf file time variable which is not correctly read by R terra package.
#' UKCP18 time values expressed as hours since 1/1/1970 12.00. Output used by function `correct_ukcp_dates()`
#' @examples
get_ukcp18_dates<-function(ncfile){
  netcdf_data <-ncdf4::nc_open(ncfile)
  time_hours <- ncdf4::ncvar_get(netcdf_data,"time")
  ncdf4::nc_close(netcdf_data)
  years<-floor(time_hours/(360*24))+1970
  months<-ceiling((time_hours-(years-1970)*(360*24)) / (30*24) )
  days<-ceiling( (time_hours-((years-1970)*(360*24)) - ((months-1)*30*24) )/ 24 )
  hours<-time_hours - ((years-1970)*(360*24)) - ((months-1)*30*24) - ((days-1)*24)
  ukcp_dates<-paste(years,sprintf("%02d",months),sprintf("%02d",days),sep="-")
  return(ukcp_dates)
}


#' @title Corrects time series of UKCP18 360 day year data to valid calendar dates
#' @description Converts 12 x 30 day monthly data provided by UKCP18 to data corresponding to actual calendar dates
#' by reassigning invalid dates (eg 29-30 Feb). Does NOT add missing dates (eg 31st May). Output used by `fill_calendar_data()`
#' @param ukcp_dates vector of POSIXlt UKCP18 dates  as returned by 'get_ukcp18_dates()'.
#' @return a vector of valid date strings in the form 'yyy-dd-mm'
#' @import lubridate
#' @export
#' @examples
correct_ukcp_dates<-function(ukcp_dates){
  years<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 1))
  months<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 2))
  days<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 3))

  # Shift March day by plus one
  sel<-which(months==3)
  days[sel] <- days[sel] + 1

  # Shift Feb day by minus 1
  sel <- which(months==2)
  days[sel] <-  days[sel] - 1

  # 1st of Feb as 31 Jan
  sel <- which(days == 0)
  days[sel] <- 31
  months[sel] <- 1

  # 29 of Feb as 1 Mar
  sel <- which(months == 2 & days == 29)
  months[sel] <- 3
  days[sel] <- 1

  # Construct new date variable
  real_dates<-as.POSIXlt(ISOdate(years, months, days))
  real_dates<-trunc(real_dates,"day")
  return(real_dates)
}

#' @title Converts 360-day yearly timeseries of UKCP18 raster data to actual calendar data
#' @description
#' Reassigns data and interpolates missing days of UKCP18 (360 day years) SpatRast data to create a new set of rasters corresponding to actual calendar dates.
#' @param ukcp_r SpatRaster timeseries stack loaded from UKCP18 .ncdf file
#' @param real_dates vector of valid date strings, as output by `correct_ukcp_dates()` corresponding to the time period of 'ukcp_r'
#' @param testplot logical if TRUE plots example data of missing then interpolated days
#' @param plotdays the example days of years to plot DO NOT CHANGE
#' @return a SpatRaster timeseries corresponding to actual calendar dates with invalid, no missing/empty days
#' @import terra
#' @export
#' @examples
fill_calendar_data<-function(ukcp_r, real_dates, testplot=FALSE, plotdays=c(89:91,242:244)){
  # Assign real dates to layer names and time values
  terra::time(ukcp_r)<-real_dates
  names(ukcp_r)<-real_dates
  # Insert and fill missing date layers of spatrast - could use zoo::na.approx or na.spline
  ukcp_r<-terra::fillTime(ukcp_r)
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'before filling'))
  ukcp_r<-terra::approximate(ukcp_r,method="linear")
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'after filling'))
  return(ukcp_r)
}


#' @title Convert units and values of a SpatRast
#' @param r SpatRaster with defined units
#' @param to_unit units of output SpatRaster
#' @return SpatRaster in new units
#' @import terra, units
#' @export
#' @examples
change_rast_units<-function(r,to_unit){
  in_unit<-terra::units(r)
  # Get and convert values
  v<-units::set_units(terra::values(r),terra::units(r)[1],mode="standard")
  v<-units::set_units(v,to_unit,mode="standard")
  # return new raster
  newr<-terra::setValues(r,v,keepunits=FALSE)
  terra::units(newr)<-to_unit
  return(newr)
}


#' @title Find UKCP18 decades containing timeseries of data
#' @description Returns UKCP18 decade text (used in file names) containing required data defined by start and end dates
#' @param collection UKCP18 collection, either global ('land-gcm') or regional ('land-rcm'), as file naming conventions vary between each collection
#' @param startdate POSIXlt class defining starting date
#' @param enddate POSIXlt class defining starting date
#' @return a vector of strings corresponding to the decade part of UKCP18 files containing the entire requested timeseries
#' @export
#' @examples
#' find_ukcp_decade('land-rcm',as.POSIXlt("2030/01/01"),as.POSIXlt("2039/12/31"))
#'
find_ukcp_decade<-function(collection=c('land-gcm','land-rcm'),startdate,enddate){
  collection<-match.arg(collection)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  rcm_decades<-c('19801201-19901130','19901201-20001130','20001201-20101130','20101201-20201130',
                 '20201201-20301130','20301201-20401130','20401201-20501130','20501201-20601130',
                 '20601201-20701130','20701201-20801130')
  gcm_decades<-c('19791201-19891130','19891201-19991130','19991201-20091130','20091201-20191130',
                 '20191201-20291130','20291201-20391130','20391201-20491130','20491201-20591130',
                 '20591201-20691130','20691201-20791130')
  if(collection=='land-gcm') ukcp_decades<-gcm_decades else ukcp_decades<-rcm_decades
  decade_start<-ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 1))
  decade_end<-ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 2))
  decades<-ukcp_decades[which(decade_end>startdate & decade_start<enddate)]
  return(decades)
}

#' @title convert UKCP18  ncdf4 files to format required for model
#' @description Converts UKCP18 global or regional land data in the form of netCDF4 file (as returned
#' by [ukcp18_downscale()] to the correct data form required for subsequent modelling.
#' @param dir_ukcp directory holding ALL ncdf files of UKCP18 data required to extract the requested variables and timeseries.
#' @param dtm SpatRaster object of elevations covering the geographical extent required (see details)
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasts.
#' @param sampleplot if TRUE plots examples of interpolated dates when converting from 360 to 366 day years
#' @param startdate POSIXlt class defining start date of required timeseries
#' @param enddate POSIXlt class defining end date of required timeseries
#' @param collection text string defining UKCP18 collection, either 'land-gcm' or 'land-rcm'
#' @param domain text string defining UKCP18 domain, either 'uk' or 'eur'(land-rcm collection only) or 'global'
#' @param member string defining the climate model member to be used for the timeseries. Available members vary between UKCP18 collections.
#' @param ukcp_vars UKCP18 variable names to be extracted DO NOT CHANGE
#' @param ukcp_units units of the UKCP18 variables extracted DO NOT CHANGE
#' @param output_units units required for output DO NOT CHANGE CHECK CORRECT
#'
#' @return a list of the following:
#' \describe{
#'   \item{tme}{POSIXlt object of times corresponding to climate observations}
#'   \item{climarray}{ either a list of arrays or a list of SpatRasts of daily weather variables
#'   \item{dtmc}{a coarse resolution digital elevation dataset matching the resolution and crs of the input
#'   climate data, and an extent matching `dtm`}
#' }
#' @import terra
#' @export
#' @details The function converts 360 day UKCP18 to actual calendar data by reassignment and linear interpolation of missing data. Returned climate data will be at the resolution, corrdinate reference system of the UKCP18 data requested, which can vary between domains and collections.
#' The extent will be provided by the `dtm` which, if necessary, will be reprojected and resampled to the CRS and resolution of UKCP18 data requested.
#' For regional UKCP18 data, it is recommended that 'dtm' is derived from the original orography data available for download.
#' @examples
#' dir_out <- tempdir()
#' download_ukcp18(dir_out,'ceda_username','cedapwd'ceda_password',as.POSIXlt('2018-05-01'),as.POSIXlt('2018-05-31'),'land-rcm','uk','rcp85',c('01','02','03'),c('clt','hurs','huss','pr','psl','rls','rss','tas','tasmax','tasmin','uas','vas'))
#' climdata<-ukcp18toclimarray()
ukcp18toclimarray <- function(dir_ukcp, dtm,
                              toArrays=TRUE, sampleplot=FALSE, startdate, enddate,
                              collection=c('land-gcm','land-rcm'),
                              domain=c('uk','eur','global'),
                              member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                                      '16','17','18','19','20','21','22','23','24','25','26','27','28'),
                              ukcp_vars=c('clt','hurs','huss','pr','psl','rls','rss',
                                         'tasmax','tasmin','uas','vas'),
                              ukcp_units=c('%','%','1','mm/day','hPa','watt/m^2','watt/m^2',
                                           'degC','degC','m/s','m/s'),
                              output_units=c('%','%','1','mm/day','kPa','watt/m^2','watt/m^2',
                                              'degC','degC','m/s','m/s') ){
  # Parameter check
  collection<-match.arg(collection)
  domain<-match.arg(domain)
  ukcp_vars<-match.arg(ukcp_vars,several.ok=TRUE)
  ukcp_units<-match.arg(ukcp_units,several.ok=TRUE)
  output_units<-match.arg(output_units,several.ok=TRUE)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  if(!class(dtm)[1] == "SpatRaster") stop("dtm parameter is NOT a SpatRaster!!!")

  # Check member in chosen collection
  member<-match.arg(member)
  if(collection =='land-rcm' &
     !member %in% c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15')) stop(paste("Model member",member,"NOT available in land-rcm collection - ignoring!!"))

  # Check dir_ukcp exists
  if(!dir.exists(dir_ukcp)) stop(paste("Directory does NOT exist:",dir_ukcp))

  # Identify which decades are required
  decades<-find_ukcp_decade(collection,startdate,enddate)

  # Get ukcp files names and check all exist in dir_ukcp before loading
  if(collection=='land-gcm') collres<-'60km' else collres<-'12km'
  ncfiles<-do.call(paste0, c(expand.grid(ukcp_vars,'_rcp85_',collection,'_',domain,'_',collres,'_',member,'_day_', decades,".nc")))
  ncfiles<-file.path(dir_ukcp,ncfiles)
  not_present<-which(!file.exists(ncfiles))
  if (length(not_present)>0) stop(paste("File NOT present in ukcp directory: ",ncfiles[not_present]," ") )

  # Get res, crs from sample nc raster as can vary between collections and domains - lat/lon or uk OS
  sample_ukcp18<-terra::rast(ncfiles[1], subds=strsplit(basename(ncfiles[1]),'_')[[1]][1])[[1]]
  ukcp_crs<-terra::crs(sample_ukcp18)
  ukcp_res<-terra::res(sample_ukcp18)
  ukcp_e<-terra::ext(sample_ukcp18)

  # Create output list
  output_list<-list()

  # Create coarse-resolution dtm to use as template for resampling - weighted !!!
  dtm <- subst(dtm,NA, 0)
  output_list$dtmc<-terra::trim(terra::resample(terra::project(dtm,ukcp_crs),sample_ukcp18,method='bilinear'))
  units(output_list$dtmc)<-'m'
  names(output_list$dtmc)<-'Elevation'

  # Load spatrasters from ukcp.nc file, crop to dtmc, convert to normal calendar and add to output list
  clim_list<-list()
  for (n in 1:length(ukcp_vars)){
    v<-ukcp_vars[n]
    print(v)
    ukcp_u<-ukcp_units[n]
    out_u<-output_units[n]
    var_r<-terra::rast()
    for (d in decades){
      ncfile<-file.path(dir_ukcp,paste0(v,'_rcp85_',collection,'_',domain,'_',collres,'_',member,'_day_', d,".nc") )
      print(paste("Loading",ncfile))

      # Load and crop data in native crs then project to dtmc crs and recrop
      r <- terra::rast(ncfile, subds=v)
      r<-terra::crop(r,output_list$dtmc)

      # If requested then convert to real calendar dates and fill missing dates
      ukcp_dates<-get_ukcp18_dates(ncfile)
      real_dates<-correct_ukcp_dates(ukcp_dates)
      r<-fill_calendar_data(r, real_dates, testplot=sampleplot)

      # Correct units
      terra::units(r)<-ukcp_u

      # PROBLEM with u v wind vectors in global data - based on DIFFERENT grid to other variables!!!
      if(terra::ext(r)!=terra::ext(output_list$dtmc)){
        print("Extents of input data DIFFERENT - correcting!!")
        r<-terra::resample(r,output_list$dtmc)
      }
      # Join if multiple decades of data
      terra::add(var_r)<-r
    }
    # Check & convert units, check all rasters geoms are same
    if(ukcp_u!=out_u) var_r<-change_rast_units(var_r, out_u)
    if(!terra::compareGeom(output_list$dtmc,var_r)) warning(paste(v,"Spatrast NOT comparable to output_list$dem!!") )
    clim_list[[v]]<-var_r
  }

  # Calculate derived variables
  clim_list$windspeed<-.rast(sqrt(as.array(clim_list$uas)^2+as.array(clim_list$vas)^2)*log(67.8*2-5.42)/log(67.8*10-5.42),output_list$dtmc) # Wind speed (m/s)
  clim_list$winddir<-.rast(as.array((terra::atan2(clim_list$uas,clim_list$vas)*180/pi+180)%%360),output_list$dtmc) # Wind direction (deg from N - from)

  # Set units, times and layer names
  units(clim_list$windspeed)<-'m/s'
  units(clim_list$winddir)<-'deg'
  terra::time(clim_list$windspeed)<-terra::time(clim_list$uas)
  terra::time(clim_list$winddir)<-terra::time(clim_list$uas)
  names(clim_list$windspeed)<-terra::time(clim_list$windspeed)
  names(clim_list$winddir)<-terra::time(clim_list$winddir)

  # Select and rename climate output rasts MIGHT NEED TO CHANGE THESE TO MATCH THOSE USED BY MESOCLIM FUNCTIONS
  # *** Ilya comment - will need to integrate some conversion functons and change names
  # mesclim names and variables are: temp, (or tmin / tmax), relhum, pres,swrad,     "difrad"    "lwrad"     "windspeed"
  # winddir,  prec.
  clim_list<-clim_list[c("clt","hurs","huss","pr","psl","rls","rss","tasmax","tasmin", "windspeed","winddir")]
  names(clim_list)<-c('cloud','relhum','spechum','prec','pres','lwnet','swnet','tmax','tmin','windspeed','winddir')

  # Add other data
  output_list$windheight_m<-10. # windspeed at 10 metres height
  output_list$tempheight_m<-1.5 # air temp at 1.5 metres height

  # Restrict to dates requested
  filter_times<-function(x,startdate,enddate) x[[which(time(x) >= startdate & time(x) <= enddate)]]
  for(v in names(clim_list))clim_list[[v]]<-filter_times(clim_list[[v]],startdate,enddate)

  # Generate POSIXlt object of times and add to output list
  output_list$tme<-as.POSIXlt(time(clim_list$tmax),tz="UTC")

  # Convert climate data to arrays if required
  if(toArrays) clim_list<-lapply(clim_list,as.array)

  output_list$climarray<-clim_list

  return(output_list)
}

