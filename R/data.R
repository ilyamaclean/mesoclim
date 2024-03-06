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
#'    \item{prec}{Precipitation (mm)}
#'  }
#'  \item{dtmc}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#' }
#' @source \url{https://cds.climate.copernicus.eu//}
"era5data"
#' A 0.25 degree grid resolution dataset of sea-surface temperature data
#'
#' A spatial dataset of hourly sea-surface temperatures for May 2018 for sea around
#' West Cornwall, UK, covering the area bounded by  -6.125, -4.125, 49.125, 51.125
#' (xmin, xmax, ymin, ymax) with the WGS84 lat long coordinate reference system (EPSG:4326)
#'
#' @format A PackedSpatRaster object with 8 rows, 8 columns and 744 layers
#' @source \urlhttps://cds.climate.copernicus.eu//}
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
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{climarray}{a list of arrays of climate variables with the following variables:}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \describe{
#'    \item{tmax}{Maximum Daily Temperature (deg C)}
#'    \item{tmin}{Minimum Daily Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{spechum}{Specific humidity}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swnet}{Net downward shortwave radiation (W/m^2)}
#'    \item{lwnet}{Net downward longwave radiation (W/m^2)}
#'    \item{windspeed at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#'  \item{dtmc}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#' }
#' @source \url{ftp.ceda.ac.uk}
"ukcp18rcm"

#' A list of 60km resolution UKCP18 global climate data (UK domain)
#'
#' A list of hourly UKCP18 regional climata data for May 2018 for Southwest UK as
#' returned by [ukcp18toclimarray]  where collection='land-gcm'& domain='uk'
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{climarray}{a list of arrays of climate variables with the following variables:}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \describe{
#'    \item{tmax}{Maximum Daily Temperature (deg C)}
#'    \item{tmin}{Minimum Daily Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{spechum}{Specific humidity}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swnet}{Net downward shortwave radiation (W/m^2)}
#'    \item{lwnet}{Net downward longwave radiation (W/m^2)}
#'    \item{windspeed at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#'  \item{dtmc}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#' }
#' @source \url{ftp.ceda.ac.uk}
"ukcp18gcm"

#' A list of ~60km resolution UKCP18 global climate data (global domain)
#'
#' A list of hourly UKCP18 regional climata data for May 2018 for East Australia UK as
#' returned by [ukcp18toclimarray]  where collection='land-gcm'& domain='global'.
#' Lon res = 0.833 deg; Lat res = 0.55 deg
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{tme}{POSIXlt object of dates and times}
#'  \item{climarray}{a list of arrays of climate variables with the following variables:}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \describe{
#'    \item{tmax}{Maximum Daily Temperature (deg C)}
#'    \item{tmin}{Minimum Daily Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{spechum}{Specific humidity}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swnet}{Net downward shortwave radiation (W/m^2)}
#'    \item{lwnet}{Net downward longwave radiation (W/m^2)}
#'    \item{windspeed at 2m (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#'  \item{dtmc}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
#' }
#' @source \url{ftp.ceda.ac.uk}
"ukcp18aus"
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
#'
#' @format A Dataframe object with 28 rows, 6 columns
#' @source \url{}
"ukcp18lookup"

