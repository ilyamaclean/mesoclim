% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporaldownscale.R
\name{blendtemp_hadukera5}
\alias{blendtemp_hadukera5}
\title{Blends Met Office and ERA5 data to produce hourly 1km resolution temperature data}
\usage{
blendtemp_hadukera5(tasmin, tasmax, era5t2m)
}
\arguments{
\item{tasmin}{a stacked SpatRaster of haduk daily minimum temperatures (deg C)}

\item{tasmax}{a stacked SpatRaster of haduk daily maximum temperatures (deg C)}

\item{era5t2m}{a stacked SpatRaster of hourly ERA5 temperatures (deg C or K)}
}
\description{
The function \code{blendtemp_hadukera5} ERA5 data to 1 km grid resoltuion,
calculates the diurnal cycle in each grid cell, and then adjusts this by the
maximum and minimum daily temperatures in the one km met office data.
}
\keyword{preprocess}
\keyword{temporal}
