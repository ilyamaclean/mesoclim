#' A list of ERA5 climate data
#'
#' A list of hourly ERA5 climata data for May 2018 for Cornwall, UK (lat, long boundary) as
#' returned by [era5toclimarray()]
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{dtm}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \item{temp}{Temperature (deg C)}
#'  \item{relhum}{Relative humidity (Percentage)}
#'  \item{spechum}{Specific humidity}
#'  \item{pres}{Sea-level atmospheric pressure (kPa)}
#'  \item{swnet}{Net downward shortwave radiation (W/m^2)}
#'  \item{lwnet}{Net downward longwave radiation (W/m^2)}
#'  \item{windspeed}{at 2m (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{prec}{Precipitation (mm)}
#' }
#' @source \url{https://cds.climate.copernicus.eu//}
"era5input"
#' A 0.25 degree grid resolution dataset of sea-surface temperature data
#'
#' A spatial dataset of hourly sea-surface temperatures for May 2018 for sea around
#' West Cornwall, UK, covering the area bounded by  -6.125, -4.125, 49.125, 51.125
#' (xmin, xmax, ymin, ymax) with the WGS84 lat long coordinate reference system (EPSG:4326)
#'
#' @format A PackedSpatRaster object with 8 rows, 8 columns and 744 layers
#' @source \url{https://cds.climate.copernicus.eu//}
"era5sst"
#' A 100m grid resolution landsea mask
#'
#' A spatial landsea mask for West Cornwall, UK, covering the area bounded by  145000,
#' 195000, -9000, 41000 (xmin, xmax, ymin, ymax) with the with the coordinate reference
#' system OSGB36 / British National Grid (EPSG:27700)
#'
#' @format A PackedSpatRaster object with 500 rows and 500 columns
#' @source \url{https://cds.climate.copernicus.eu//}
"landsea"
#' A list of 12km resolution UKCP18 regional climate data
#'
#' A list of hourly UKCP18 regional climata data for May 2018 for Southwest UK as
#' returned by [ukcp18toclimarray()] where collection='land-gcm' & domain='uk'
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{dtm}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \item{tmax}{Maximum Daily Temperature (deg C)}
#'  \item{tmin}{Minimum Daily Temperature (deg C)}
#'  \item{relhum}{Relative humidity (Percentage)}
#'  \item{spechum}{Specific humidity}
#'  \item{pres}{Sea-level atmospheric pressure (kPa)}
#'  \item{swnet}{Net downward shortwave radiation (W/m^2)}
#'  \item{lwnet}{Net downward longwave radiation (W/m^2)}
#'  \item{windspeed}{at 2m (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{prec}{Precipitation (mm)}
#' }
#' @source \url{ftp.ceda.ac.uk}
"ukcpinput"
#' A ~0.1 degree grid resolution dataset of UKCP18 sea-surface temperature data for NW Europe
#'
#' A spatial dataset of monthly sea-surface temperatures for April to June 2018 for sea around
#' UK & North West Europe, covering the area bounded by  -19.9, 13.1, 40.0, 65.0
#' (xmin, xmax, ymin, ymax) with the WGS84 lat long coordinate reference system (EPSG:4326)
#'
#' @format A PackedSpatRaster object with 375 rows, 297 columns and 3 layers
#' @source \url{ftp.ceda.ac.uk}
"ukcp18sst"
#' A model member lookup table for different UKCP18 collections
#'
#' Matches member number with model name and in which collections / domains it is avaiolble
#' See: https://www.metoffice.gov.uk/binaries/content/assets/metofficegovuk/pdf/research/ukcp/ukcp18-guidance-data-availability-access-and-formats.pdf
#' @format A Dataframe object with 28 rows, 6 columns
"ukcp18lookup"

