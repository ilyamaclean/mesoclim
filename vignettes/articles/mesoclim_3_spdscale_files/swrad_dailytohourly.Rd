% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporaldownscale.R
\name{swrad_dailytohourly}
\alias{swrad_dailytohourly}
\title{Daily to hourly shortwave radiation}
\usage{
swrad_dailytohourly(radsw, tme, clearsky = NA, r = r, adjust = TRUE)
}
\arguments{
\item{radsw}{\itemize{
\item an array of daily mean radiation values (W/m**2)
}}

\item{tme}{\itemize{
\item POSIXlt object of dates corresponding to radsw
}}

\item{clearsky}{\itemize{
\item optionally an array with dimensions matching radsw of daily clearsky
radiation as returned by clearskyraddaily(). Calculated if not supplied
}}

\item{r}{\itemize{
\item a terra::SpatRaster object giving the extent of radsw used for deriving the lat and long of grid cells
}}

\item{adjust}{\itemize{
\item optional logical which if TRUE ensures that, after interpolation, returned
hourly values, when averaged to daily, match the input
}}
}
\value{
an array of hourly radiation values (W/m**2)
}
\description{
Derives an array of hourly radiation values from an array of daily radiation values.
}
\details{
NB - this assumes input radiation is downward flux, not net radiation (as assumed in UKCP)
to get from net to downward flux we need to recognise that rswnet = (1-alb)*radsw, so
radsw = rswnet/(1-alb), where alb is white sky albedo. White-sky albedo changes as a function
of solar angle, but in a manner dependent on ground reflectance, leaf area, leaf inclination
angles and leaf transmittance and the ratio of diffuse and direct. There are too
many vegetation parameter unknowns to reverse engineer, so suggest ignoring this.
discrepancies probably quite minor expect in areas with very low cover and will be handled
mostly by bias correction anyway
~~ * Need to spline interpolate clear-sky fraction (bounding by 1 and 0) and
~~   then calculate clear-sky radiation
~~ * Need to spline interpolate sky emissvity (bounding by 1 and 0) and
~~   then calculate longwave
}
\examples{
climdata<- read_climdata(system.file('data/ukcpinput.rds',package='mesoclim'))
swh<-swrad_dailytohourly(climdata$swrad,climdata$tme,r=climdata$dtm)
}
\keyword{temporal}
