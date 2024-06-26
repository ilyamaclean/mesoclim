% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{era5input}
\alias{era5input}
\title{A list of ERA5 climate data}
\format{
a list with the following elements:
\describe{
\item{dtm}{a wrapped SpatRast object of elevations (m) matching the extent and resolution of climate data}
\item{tme}{POSIXlt object of dates and times}
\item{windheight_m}{numeric value in metres of wind height above ground}
\item{tempheight_m}{numeric value in metres of temperature height above ground}
\item{temp}{Temperature (deg C)}
\item{relhum}{Relative humidity (Percentage)}
\item{spechum}{Specific humidity}
\item{pres}{Sea-level atmospheric pressure (kPa)}
\item{swnet}{Net downward shortwave radiation (W/m^2)}
\item{lwnet}{Net downward longwave radiation (W/m^2)}
\item{windspeed}{at 2m (m/s)}
\item{winddir}{Wind direction (decimal degrees)}
\item{prec}{Precipitation (mm)}
}
}
\source{
\url{https://cds.climate.copernicus.eu//}
}
\usage{
era5input
}
\description{
A list of hourly ERA5 climata data for May 2018 for Cornwall, UK (lat, long boundary) as
returned by \code{\link[=era5toclimarray]{era5toclimarray()}}
}
\keyword{datasets}
