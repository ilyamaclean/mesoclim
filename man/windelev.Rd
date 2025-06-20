% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climatefunctions.R
\name{windelev}
\alias{windelev}
\title{derive wind terrain adjustment coefficient}
\usage{
windelev(dtmf, dtmm, dtmc, wdir, uz = 2)
}
\arguments{
\item{dtmf}{a high-resolution SpatRast of elevations}

\item{dtmm}{a medium-resolution SpatRast of elevations covering a larger area
than dtmf (see details)}

\item{dtmc}{a coarse-resolution SpatRast of elevations usually matching
the resolution of climate data used for downscaling (see details)}

\item{wdir}{wind direction (from, decimal degrees).}

\item{uz}{height above ground (m) of wind speed measurement}
}
\value{
a SpatRast of wind adjustment coefficients matching the resoltuion,
coordinate reference system and extent of \code{dtmf}.
}
\description{
The function \code{windelev} is used to spatially downscale wind, and
adjusts wind speed for elevation and applies a terrain shelter coefficient for
a specified wind direction.
}
\details{
Elevation effects are derived by sampling the dtms at intervals in
an upwind direction, determining the elevation difference form each focal cell and
performing a standard wind-height adjustment. Terrain sheltering is computed
from horizon angles following the method detailed in Maclean et al (2019) Methods
Ecol Evol 10: 280-290. By supplying three dtms, the algorithm is able to account for
elevation differences outside the boundaries of \code{dtmf}. The area covered by \code{dtmm} was
extend at least one \code{dtmc} grid cell beyond \code{dtmf}. Elevations must be in metres.
The coordinate reference system of \code{dtmf} must be such that x and y are also in metres.
\code{dtmm} and \code{dtmc} are reprojected to match the coordinate reference system of \code{dtmf}.
}
\examples{
dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
climdata<-read_climdata(mesoclim::ukcpinput)
wc <- windelev(dtmf, dtmm, climdata$dtm, wdir = 270)
terra::plot(wc)
}
\seealso{
\code{\link[=winddownscale]{winddownscale()}}
}
\keyword{spatial}
