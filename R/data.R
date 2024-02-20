#' A 100m grid resolution digital elevation dataset
#'
#' A spatial dataset of elevations in (m) for the Lizard Peninsula, COrnwall, UK,
#' covering the area bounded by 160000, 181400, 11300, 30000  (xmin, xmax, ymin, ymax)
#' with the coordinate reference system OSGB36 / British National Grid (EPSG:27700)
#'
#' @format A PackedSpatRaster object with 197 rows and 214 columns
#' @source \url{http://www.tellusgb.ac.uk/}
"dtmf"
#' A 1km grid resolution digital elevation dataset
#'
#' A spatial dataset of elevations in (m) for West COrnwall, UK,
#' covering the area bounded by 20000, 220000, -34000, 66000  (xmin, xmax, ymin, ymax)
#' with the coordinate reference system OSGB36 / British National Grid (EPSG:27700)
#'
#' @format A PackedSpatRaster object with 100 rows and 100 columns
#' @source \url{http://www.tellusgb.ac.uk/}
"dtmm"


#' A list of ERA5 climate data
#'
#' A list of hourly ERA5 climata data for May 2018 for Cornwall, UK (lat, long boundary) as
#' returned by [era5toclimarray()]
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{climarray}{a list of arrays of climate variables with the following variables:}
#'  \describe{
#'    \item{temp}{Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2)}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'  }
#'  \item{dtmc}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#' }
#' @source \url{https://cds.climate.copernicus.eu//}
"era5data"

