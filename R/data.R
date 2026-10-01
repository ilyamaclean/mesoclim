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
#' A list of daily UKCP18 regional climate model data (rcp85, member 01) for May 2018 for
#' Southwest UK, as returned by [ukcp18toclimarray()] with collection='land-rcm' and domain='uk'.
#' Covers a grid of 4 x 3 cells of 12 km (x 144000-192000, y 0-36000, OSGB), which includes
#' the example DTMs `dtmf.tif` and `lizard50m.tif` in `inst/extdata/dtms`.
#'
#' @format a list with the following elements (climate variables as 3D arrays):
#' \describe{
#'  \item{dtm}{a wrapped SpatRast of UKCP18 land elevations (m) matching the climate data; NA where UKCP18 classes a cell as sea}
#'  \item{tme}{POSIXlt object of daily dates}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \item{cloud}{Cloud cover (Percentage)}
#'  \item{relhum}{Relative humidity (Percentage)}
#'  \item{prec}{Precipitation (mm/day)}
#'  \item{pres}{Sea-level atmospheric pressure (kPa)}
#'  \item{lwrad}{Downward longwave radiation (W/m^2)}
#'  \item{swrad}{Downward shortwave radiation (W/m^2)}
#'  \item{tmax}{Maximum daily temperature (deg C)}
#'  \item{tmin}{Minimum daily temperature (deg C)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
#' @source \url{https://catalogue.ceda.ac.uk}. Created by `data-raw/ukcp_sample_data.R`.
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
#' A list of 12km resolution UKCP18 model future climate data
#'
#' A list of daily UKCP18 regional climate model data (rcp85, member 01) for May 2030 for
#' Southwest UK, as returned by [ukcp18toclimarray()] with collection='land-rcm' and domain='uk'.
#' Covers the same 4 x 3 grid of 12 km cells as [ukcpinput].
#'
#' @format a list with the following elements (climate variables as PackedSpatRasters):
#' \describe{
#'  \item{dtm}{a wrapped SpatRast of UKCP18 land elevations (m) matching the climate data; NA where UKCP18 classes a cell as sea}
#'  \item{tme}{POSIXlt object of daily dates}
#'  \item{windheight_m}{numeric value in metres of wind height above ground}
#'  \item{tempheight_m}{numeric value in metres of temperature height above ground}
#'  \item{cloud}{Cloud cover (Percentage)}
#'  \item{relhum}{Relative humidity (Percentage)}
#'  \item{prec}{Precipitation (mm/day)}
#'  \item{pres}{Sea-level atmospheric pressure (kPa)}
#'  \item{lwrad}{Downward longwave radiation (W/m^2)}
#'  \item{swrad}{Downward shortwave radiation (W/m^2)}
#'  \item{tmax}{Maximum daily temperature (deg C)}
#'  \item{tmin}{Minimum daily temperature (deg C)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
#' @source \url{https://catalogue.ceda.ac.uk}. Created by `data-raw/ukcp_sample_data.R`.
"ukcpfuture"
#' Sample bias-correction model list for UKCP18 RCM member 01
#'
#' A list of bias-correction models for UKCP18 regional climate model member 01 for
#' Southwest UK, fitted for the land cells of [ukcpinput] (May 2018) against HadUK-Grid 1 km
#' observations aggregated to 12 km, as in the bias correction vignette. Models are
#' available for minimum temperature, maximum temperature and precipitation. Other climdata variables (e.g. relhum, pres, swrad, lwrad,
#' windspeed) are not included; [biascorrect_climdata()] will apply correction only
#' to variables present in the list and warn about any that are missing.
#'
#' @format A named list with three elements:
#' \describe{
#'  \item{tmin}{Bias-correction model for minimum daily temperature (deg C), of class \code{biascorrectmodels}}
#'  \item{tmax}{Bias-correction model for maximum daily temperature (deg C), of class \code{biascorrectmodels}}
#'  \item{prec}{Bias-correction model for daily precipitation (mm): a list of wrapped SpatRasters \code{mu_tot} and \code{mu_frac} as returned by [precipcorrect()] with \code{mod_out = TRUE}}
#' }
#' @source Derived from UKCP18 RCM and HadUK-Grid data by `data-raw/ukcp_sample_data.R`. See [biascorrect()] and [precipcorrect()] for model-fitting functions.
"model_list"
