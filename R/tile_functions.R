
#' @title Spatially downscale all climate variables by tiles
#' @description Spatially downscales coarse-resolution climate data over large areas by:
#' (i) Downscaling in monthly steps
#' (ii) Dividing spatial area into a set of user defined tiles
#' @param climdata a `climdata` model object containing climate data of the same format as `era5climdata`
#' @param sst a SpatRast of sea-surface temperature data (deg C) that overlaps with climdata$tme
#' @param dtmf a high-resolution SpatRast of elevations
#' @param dtmm a medium-resolution SpatRast of elevations covering a larger area
#' than dtmf (only needed for coastal effects - see details under [temphrly_downscale()]).
#' @param basins optionally, a fine-resolution SpatRast of basins as returned by [basindelin()]
#' matching the coordinate reference system and extent of `dtmf`. Calculated if
#' not supplied.
#' @param wca optionally, an array of wind shelter coeffs in 8 directions
#' matching the x,y dimensions of `dtmf` as calculated by `calculate_windcoeffs`. Calculated if
#' not supplied.
#' @param cad optional logical indicating whether to calculate cold-air drainage effects
#' @param coastal optional logical indicating whether to calculate coastal effects
#' @param thgto height above ground of temperature output.
#' @param whgto height above ground of wind speed output.
#' @param include_tmean if TRUE and daily data will output mean daily temperature based on hourly downsccaling
#' @param rhmin minimum relative humidity (set to avoid relative humidity dropping too low
#' during down-scaling). Default 20 (percent).
#' @param pksealevel optional logical indicating whether input pressure data represent
#' sea-level pressure (default TRUE).
#' @param pathsim optional logical indicating whether to simulate cloud cover and
#' rainfall patchiness see details under [swdownscale()] and [precipdownscale()].
#' @param difr optional logical indicating whether to return diffuse radiation (
#' default TRUE).
#' @param terrainshade optional logical indicating whether to account for terrain shading
#' during shortwave radiation downscale.
#' @param precipmethod One of `Tps` or `Elev` indicating whether to account for
#' elevation effects using a Thin-plate spline model or en emprical adjustment
#' (see details under [precipdownscale()])
#' @param fast optional logical indicating whether to use fast Thin-plate spline
#' down-scaling (see details under [precipdownscale()]).
#' @param noraincut numeric value below which low precipitation amounts are set to
#' zero (see details under [precipdownscale()]).
#' @param toArrays if TRUE climate variables are output as 3D arrays if FALSE as Spatarasters
#' @param overlap - overlap of adjacent tiles (in metres or projection unit)
#' @param sz - size of tiles (in metres or projection unit)
#' @return a list of multi-layer SpatRast or 3D arrays of downscaled climate variables and metadata as follows:
#' \describe{
#'    \item{dtm}{Digital elevation of downscaled area in metres (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C) - is replaced with `tmax` and `tmin` if input data is daily and `tmean` if `include_tmean`=TRUE.}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2) - not returned unless `terrainshade = TRUE`}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed}{at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @keywords spatial
#' @examples
#'  dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#'  dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#'  climdata<-read_climdata(mesoclim::ukcpinput)
#'  sst<-terra::unwrap(mesoclim::ukcp18sst)
#' # Less efficient over such a small area but demonstrates use of multiple tiles
#'  mesodat<-spatialdownscale_tiles(climdata,sst,dtmf,dtmm,include_tmean=TRUE,noraincut=0.01, overlap=200, sz=2000)
#' for(n in 5:length(mesodat)) terra::plot(mesodat[[n]][[12]],main=paste(names(mesodat)[n],mesodat$tme[n]))
#' # A single tile can be used but in monthly steps by setting sz to greater than x and y dimensions of dtmf
#'  mesodat<-spatialdownscale_tiles(climdata,sst,dtmf,dtmm,include_tmean=TRUE,noraincut=0.01, overlap=0, sz=20000)
spatialdownscale_tiles<-function(climdata, sst, dtmf, dtmm = NA, basins = NA, wca=NA, skyview=NA, horizon=NA,
                                 cad = TRUE,coastal = TRUE, thgto =2, whgto=2,include_tmean=FALSE,
                                 rhmin = 20, pksealevel = TRUE, patchsim = FALSE, terrainshade = TRUE,
                                 precipmethod = "Elev",fast = TRUE, noraincut = 0, toArrays=FALSE, overlap=1000, sz=10000){

  # Convert variables - unpack any wrapped spatRasters and convert arrays to spatraster
  input_class<-lapply(lapply(climdata,class),"[",1)
  if(any(input_class=="PackedSpatRaster")) climdata[which(input_class=="PackedSpatRaster")]<-lapply(climdata[which(input_class=="PackedSpatRaster")],unwrap)
  if(any(input_class=="array")) climdata[which(input_class=="array")]<-lapply(climdata[which(input_class=="array")],.rast,tem=climdata$dtm)
  if(inherits(sst,"PackedSpatRaster")) sst<-unwrap(sst)
  dtmc<-climdata$dtm

  # Calculate topographical properties -for whole area if not supplied
  if(is.logical(wca)) wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)
  if(is.logical(basins) & cad==TRUE) basins<-basindelin(dtmf, boundary = 2)
  if(is.logical(skyview) | is.logical(horizon)){
    results<-calculate_terrain_shading(dtmf,steps=24,toArrays=FALSE)
    skyview<-results$skyview
    horizon<-results$horizon
  }

  # Calculate overlapping tile set (nb size can vary)
  tileset<-create_overlapping_tiles(dtmf,overlap,sz)
  tiles<-tileset[["tile_extents"]][which(tileset$tile_land=="y")]
 
  # Calculate year/months of data supplied
  start<-climdata$tme[1]
  end<-climdata$tme[length(climdata$tme)]
  yrs<-unique(year(climdata$tme))
  message(paste("Running for",length(tiles),"tiles from",format(start),"to",format(end)))

  allmonths<-list()
  t0<-now()
  for(yr in yrs){
    yrstart<-as.POSIXlt(paste0(yr,'/01/01'),tz=tz(climdata$tme))
    yrend<-as.POSIXlt(paste0(yr,'/12/31'),tz=tz(climdata$tme))
    climdata_y<-subset_climdata(climdata,yrstart,yrend)
    mnths<-unique(month(climdata_y$tme))
    for(m in mnths){
      message(paste0("Downscaling ",sprintf("%02d", m),"/",yr))
      # Get this yr/month data
      sdatetime<-as.POSIXlt(paste0(yr,'/',sprintf("%02d", m),'/01'),tz=tz(climdata$tme))
      dys<-lubridate::days_in_month(sdatetime)
      edatetime<-as.POSIXlt(paste0(yr,'/',sprintf("%02d", m),'/',sprintf("%02d", dys)),tz=tz(climdata$tme))
      climdata_m<-subset_climdata(climdata_y,sdatetime,edatetime)

      # Downscale by tile
      mesoclimate_tiles<-list()
      for(n in 1:length(tiles)){
        # for(n in 1:length(tiles[c(1,2,4,5,7,8)])){
        t<-tiles[[n]]
        print(paste("Downscaling tile",n))
        xmx<- t$xmax
        xmn<- t$xmin
        ymx<- t$ymax
        ymn<- t$ymin
        dtmf_tile<-crop(dtmf,t)
        if(cad) basins_tile<-crop(basins,t)
        wca_tile<-.is(crop(.rast(wca,dtmf),t))
        sky_tile<-crop(skyview,t)
        hor_tile<-crop(horizon,t)

        mesoclimate<-spatialdownscale(climdata=climdata_m, sst=sst, dtmf=dtmf_tile, dtmm=dtmm,
                                      basins = basins_tile, wca=wca_tile, skyview=sky_tile, horizon=hor_tile,
                                      cad=cad, coastal=coastal, thgto=thgto, whgto=whgto, include_tmean=include_tmean,
                                      rhmin=rhmin, pksealevel=pksealevel, patchsim=patchsim,
                                      terrainshade=terrainshade, precipmethod=precipmethod, fast=fast, noraincut=noraincut)

        mesoclimate_tiles[[length(mesoclimate_tiles)+1]]<-mesoclimate
      } # downscale tiles

      # Recreate single month data of merged tiles
      mesomonth<-list()
      mesomonth$dtm<-dtmf
      mesomonth$tme<-mesoclimate_tiles[[1]]$tme
      mesomonth$windheight_m<-mesoclimate_tiles[[1]]$windheight_m
      mesomonth$tempheight_m<-mesoclimate_tiles[[1]]$tempheight_m
      append_vars<-names(mesoclimate_tiles[[1]])[c(5:length(mesoclimate_tiles[[1]]))]
      if(length(mesoclimate_tiles)==1) for(v in append_vars) mesomonth[[v]]<-mesoclimate_tiles[[1]][[v]]
      if(length(mesoclimate_tiles)>1) for(v in append_vars) mesomonth[[v]]<-do.call(merge, lapply(mesoclimate_tiles,`[[`, v))
      for(v in append_vars) names(mesomonth[[v]])<-terra::time(mesomonth[[v]])

      allmonths[[length(allmonths)+1]] <-mesomonth
      #print(paste("Time for downscaling ALL tiles for ALL months of year",yr,"=",format(now()-t0)))
    } # month
  } # years

  # Reconstruct final output combining all montly timesteps
  mesoout<-list()
  mesoout$dtm<-allmonths[[1]]$dtm
  mesoout$windheight_m<-allmonths[[1]]$windheight_m
  mesoout$tempheight_m<-allmonths[[1]]$tempheight_m
  mesoout$tme<-do.call(c,(lapply(allmonths,"[[","tme")))
  append_vars<-names(allmonths[[1]])[c(5:length(allmonths[[1]]))]
  for(v in append_vars) mesoout[[v]]<-rast(unlist(lapply(allmonths,"[",v)))
  for(v in append_vars) names(mesoout[[v]])<-terra::time(mesoout[[v]])

  return(mesoout)
} # function

#' Create overlapping tile set
#' @param template.r = 1km land mask spatraster to inform whether land present
#' @param overlap = overlap in metres of each tile (x & y)
#' @param sz = base size of tiles in metres - some output tiles may be smaller
#'
#' @returns named list where "tile_extents" holds the terra extent of each tile and
#' "tile_land" is BOOLEAN of whether land cells within tile (based on template.r).
#' Returned list ordered by cols.
#' @export
#'
#' @examples
#' r<-terra::rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))
#' testtiles<-create_overlapping_tiles(r,overlap=200,sz=2000)
#' terra::plot(r)
#' for(t in testtiles$tile_extents) terra::plot(terra::vect(terra::ext(t)),add=TRUE)
#' testtiles<-create_overlapping_tiles(r,overlap=200,sz=20000)
create_overlapping_tiles<-function(template.r,overlap=1000,sz=10000){
  # If tile size > template.r return single tile with message
  if(sz>ncol(template.r)*res(template.r)[1] & sz>nrow(template.r)*res(template.r)[2]){
    warning("Requested tile size larger than input area - returning a single tile of whole area!")
    tileset<-list("tile_extents"=list(ext(template.r)),"tile_land"='y')
    return(tileset)
  }
  xmax<-ext(template.r)[2]
  xmin<-ext(template.r)[1]
  if(sz%%res(template.r)[1]!=0) warning("Choice of tile size is NOT divisible by resolution of template.r!!")

  xtiles<-(xmax-xmin)%/%sz
  xstart<-seq(xmin,xmin+(sz*xtiles-1),sz-overlap)
  xend<-(xstart+sz)
  xrem<-xmax-xend[length(xend)]
  if(xrem< -(0.5*sz)){
    xstart<-xstart[1:(length(xstart)-1)]
    xend<-xend[1:(length(xend)-1)]
    xend[length(xend)]<-xmax
  }
  if(xrem<=0.5*sz & xrem>-(0.5*sz)){
    xend[length(xend)]<-xmax
  }
  if(xrem>0.5*sz){
    xstart<-c(xstart,xend[length(xend)]-overlap)
    xend<-c(xend,xmax)
  }

  ymax<-ext(template.r)[4]
  ymin<-ext(template.r)[3]
  ytiles<-(ymax-ymin)%/%sz
  ystart<-seq(ymin,ymin+(sz*ytiles-1),sz-overlap)
  yend<-(ystart+sz)
  yrem<-ymax-yend[length(yend)]
  if(yrem< -(0.5*sz)){
    ystart<-ystart[1:(length(ystart)-1)]
    yend<-yend[1:(length(yend)-1)]
    yend[length(yend)]<-ymax
  }
  if(yrem<=0.5*sz & yrem> -(0.5*sz)){
    yend[length(yend)]<-ymax
  }
  if(yrem>(0.5*sz)){
    ystart<-c(ystart,yend[length(yend)]-overlap)
    yend<-c(yend,ymax)
  }

  elist<-list()
  etype<-c()
  for(x in 1:length(xstart)){
    for(y in 1:length(ystart)){
      e<-ext(xstart[x],xend[x],ystart[y],yend[y])
      r<-crop(template.r,e)
      elist<-c(elist,e)
      if(!all(is.na(values(r)))){
        etype<-c(etype,'y')
      } else etype<-c(etype,'n')
    }
  }
  tileset<-list("tile_extents"=elist,"tile_land"=etype)
  return(tileset)
}
