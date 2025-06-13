
#' @title Spatially downscale all climate variables by tiles
#' @description Spatially downscales coarse-resolution climate data over large areas by tiling.
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
  if(is.logical(basins)) basins<-basindelin(dtmf, boundary = 2)
  if(is.logical(skyview) | is.logical(horizon)){
    results<-calculate_terrain_shading(dtmf,steps=24,toArrays=FALSE)
    skyview<-results$skyview
    horizon<-results$horizon
  }

  # Calculate overlapping tile set (nb size can vary)
  tileset<-create_overlapping_tiles(dtmf,overlap,sz)
  tiles<-tileset$tile_extents[which(tileset$tile_land=="y")]

  # Calculate year/months of data supplied
  start<-climdata$tme[1]
  end<-climdata$tme[length(climdata$tme)]
  yrs<-unique(year(climdata$tme))
  message(paste("Running for",length(tiles),"tiles from",format(start),"to",format(end)))

  allmonths<-list()
  t0<-now()
  for(y in yrs){
    yrstart<-as.POSIXlt(paste0(y,'/01/01'),tz=tz(climdata$tme))
    yrend<-as.POSIXlt(paste0(y,'/12/31'),tz=tz(climdata$tme))
    climdata_y<-subset_climdata(climdata,yrstart,yrend)
    mnths<-unique(month(climdata_y$tme))
    for(m in mnths){
      message(paste0("Downscaling ",sprintf("%02d", m),"/",y))
      # Get this yr/month data
      sdatetime<-as.POSIXlt(paste0(yr,'/',m,'/01'),tz=tz(climdata$tme))
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
        basins_tile<-crop(basins,t)
        wca_tile<-.is(crop(.rast(wca,dtmf),t))
        sky_tile<-crop(skyview,t)
        hor_tile<-crop(horizon,t)

        mesoclimate<-spatialdownscale(climdata_m, sst, dtmf_tile, dtmm,
                                      basins = basins_tile, wca=wca_tile, skyview=sky_tile, horizon=hor_tile,
                                      cad, coastal, thgto, whgto, include_tmean,
                                      rhmin, pksealevel, patchsim,
                                      terrainshade, precipmethod, fast, noraincut)

        mesoclimate_tiles[[length(mesoclimate_tiles)+1]]<-mesoclimate
      } # downscale tiles

      # Recreate single month data of merged tiles
      mesomonth<-list()
      mesomonth$dtm<-dtmf
      mesomonth$tme<-mesoclimate_tiles[[1]]$tme
      mesomonth$windheight_m<-mesoclimate_tiles[[1]]$windheight_m
      mesomonth$tempheight_m<-mesoclimate_tiles[[1]]$tempheight_m
      append_vars<-names(mesoclimate_tiles[[1]])[c(5:length(mesoclimate_tiles[[1]]))]
      for(v in append_vars) mesomonth[[v]]<-do.call(merge, lapply(mesoclimate_tiles,`[[`, v))
      for(v in append_vars) names(mesomonth[[v]])<-terra::time(mesomonth[[v]])

      allmonths[[length(allmonths)+1]] <-mesomonth
      print(paste("Time for downscaling ALL tiles for ALL months of year",yr,"=",format(now()-t0)))
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
#' r<-rast(extent=ext(246700, 276700, 35500, 65400),res=c(50,50))
#' testtiles<-create_overlapping_tiles(r,overlap=500,sz=10000)
#' plot(vect(ext(dtmf)))
#' for(t in testtiles$tile_extents) plot(vect(ext(t)),add=TRUE)
create_overlapping_tiles<-function(template.r,overlap=1000,sz=10000){
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

#' Blend tile list
#' @param input_list - list of overlapping spatraster tiles
#' @return spatRaster of merged tiles
#' @details
#' Call function `mosaicblend` to within then across columns.
#'
#' @export
#'
#' @examples
blend_tile_lists<-function(input_list){
  tileext<-lapply(input_list,ext)
  # Blend cols
  rws<-unlist(unique(lapply(tileext,"ymax")))
  cls<-unlist(unique(lapply(tileext,"xmax")))

  col_list<-list()
  for(c in cls){
    sel<-which(lapply(tileext,"xmax")==c)
    blend.r<-mosaicblend(rlist=input_list[sel])
    col_list<-c(col_list,blend.r)
  }
  # Blend rows
  output_r<-rast()
  output_r<-c(output_r,mosaicblend(col_list))
  return(output_r)
}

##################### FOLLOWING MIGHT NOT BE REQUIRED ################

#' @title Mosaics a list of overlapping SpatRasters blending overlap areas
#' @description Mosaics a list of overlapping SpatRasters blending
#' the areas of overlap using a distance weighting to eliminate tiling effects
#' @param rlist a list of SpatRasters
#' @details
#' If rlist contains SpatRasters that are not overlapping the conventional terra::moasic function is used.
#' If rlist contains SpatRasters that do overlap, they should comprise a list of adjacent rasters in a single row or column.
#' Overlapping cells calulated using weighted distances to reduce tile effects
#' Function derived from: https://github.com/ilyamaclean/microclimf/tree/main
#' Called by: `blend_tile_lists`
#' @import terra
#' @export
mosaicblend <- function(rlist) {
  # order by row and then by column
  xmn<-0
  ymn<-0
  for (i in 1:length(rlist)) {
    e<-ext(rlist[[i]])
    xmn[i]<-e$xmin
    ymn[i]<-e$ymin
  }
  xmn1<-unique(xmn)
  ymn1<-unique(ymn)
  le<-min(length(xmn1),length(ymn1))
  if (le > 1) warning("rlist not a row or column. Blended mosaicing may not work")
  if (length(xmn1) > length(ymn1)) {
    o<-order(xmn)
  } else o<-order(ymn)
  rlist2<-list()
  for (i in 1:length(o)) rlist2[[i]]<-rlist[[o[i]]]
  rlist<-NULL
  rma<-rlist2[[1]]
  for (i in 2:length(rlist2)) {
    r<-rlist2[[i]]
    it<-terra::intersect(ext(rma),ext(r))
    a<-as.numeric((it$xmax-it$xmin)*(it$ymax-it$ymin))
    if (a>0) {
      rma<-.blendmosaic(rma, r)
    } else rma<-mosaic(rma,r)
  }
  return(rma)
}
#' blend two adjacent rasters that have overlap
#' called by `mosaicblend`
.blendmosaic<-function(r1, r2) {
  # run checks
  reso1<-res(r1)
  reso2<-res(r2)
  if (reso1[1] != reso2[1]) stop("resolutions must match")
  if (reso1[2] != reso2[2]) stop("resolutions must match")
  # Find whether r2 is TT, TR, RR, BR, BB, BL, LL or TL
  e1<-ext(r1)
  e2<-ext(r2)
  corner<-"ID"
  if (e2$ymax > e1$ymax & e2$xmax == e1$xmax) corner<-"TT"
  if (e2$ymax > e1$ymax & e2$xmax > e1$xmax) corner<-"TR"
  if (e2$ymax == e1$ymax & e2$xmax > e1$xmax) corner<-"RR"
  if (e2$ymax < e1$ymax & e2$xmax > e1$xmax) corner<-"BR"
  if (e2$ymax < e1$ymax & e2$xmax == e1$xmax) corner<-"BB"
  if (e2$ymax < e1$ymax & e2$xmax < e1$xmax) corner<-"BL"
  if (e2$ymax == e1$ymax & e2$xmax < e1$xmax) corner<-"LL"
  if (e2$ymax > e1$ymax & e2$xmax < e1$xmax) corner<-"TL"
  if (corner == "ID") {
    ro<-mosaic(r1,r2,fun="mean")
  } else {
    # Calculate overlap area
    if (corner == "TR" || corner == "TT") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    if (corner == "BR" || corner == "RR") eo<-ext(e2$xmin,e1$xmax,e1$ymin,e2$ymax)
    if (corner == "BL" || corner == "BB") eo<-ext(e1$xmin,e2$xmax,e1$ymin,e2$ymax)
    if (corner == "TL" || corner == "LL") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    # Calculate weights
    nx<-as.numeric((eo$xmax-eo$xmin)/res(r1)[1])
    ny<-as.numeric((eo$ymax-eo$ymin)/res(r1)[2])
    wx<-matrix(rep(seq(0,1,length.out=nx),each=ny),ncol=nx,nrow=ny)
    wy<-matrix(rep(seq(0,1,length.out=ny),nx),ncol=nx,nrow=ny)
    # Create a SpatRast of the blended area
    if (corner == "TT") w1<-wy
    if (corner == "TR") w1<-sqrt((1-wx)^2+wy^2)/sqrt(2)
    if (corner == "RR") w1<-1-wx
    if (corner == "BR") w1<-sqrt((1-wx)^2+(1-wy)^2)/sqrt(2)
    if (corner == "BB") w1<-1-wy
    if (corner == "BL") w1<-sqrt(wx^2+(1-wy)^2)/sqrt(2)
    if (corner == "LL") w1<-wx
    if (corner == "TL") w1<-sqrt(wx^2+wy^2)/sqrt(2)
    # Apply weights to raster
    nn<-dim(r1)[3]
    r1c<-crop(r1,eo)
    r2c<-crop(r2,eo)
    w1<-.rast(.rta(w1,nn),r1c[[1]])
    rb<-r1c*w1+r2c*(1-w1)
    # Clip out the overlap area
    ro<-mosaic(r1,r2)
    re<-w1*0-9999
    ro<-mosaic(ro,re,fun="min")
    ro[ro == -9999]<-NA
    # Mosaic with belnded data
    ro<-mosaic(ro,rb,fun="mean")
  }
  return(ro)
}
