% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporaldownscale.R
\name{temporaldownscale}
\alias{temporaldownscale}
\title{Downscale daily climate data to hourly}
\usage{
temporaldownscale(
  climdaily,
  adjust = TRUE,
  clearsky = NA,
  srte = 0.09,
  relmin = 10,
  noraincut = 0,
  toArrays = FALSE
)
}
\arguments{
\item{climdaily}{\itemize{
\item list of daily climate variables, dtm etc as output by \code{ukcp18toclimarray}
}}

\item{adjust}{\itemize{
\item to adjust hourly values fo mean = daily mean
}}

\item{clearsky}{\itemize{
\item clearsky values for terrain (see \code{swrad_dailytohourly})
}}

\item{srte}{\itemize{
\item a parameter controlling speed of decay of night time temperatures (see \code{temp_dailytohourly})
}}

\item{relmin}{\itemize{
\item minimum relative humidity value (see \code{hum_dailytohourly})
}}

\item{noraincut}{\itemize{
\item single numeric value indicating rainfall amounts that should
be considered as no rain (see)
}}

\item{toArrays}{\itemize{
\item if FALSE outputs SpatRasters
}}
}
\value{
a list of spatRasters or 3D arrays of hourly climate variables in same format as input \code{climdaily}
}
\description{
Downscale daily climate data to hourly
}
\examples{
climdaily<-read_climdata(mesoclim::ukcpinput)
allhrly<-temporaldownscale(climdaily, adjust = TRUE, clearsky=NA, srte = 0.09, relmin = 10, noraincut = 0)
for(n in 5:length(allhrly)) terra::plot(allhrly[[n]][[12]],main=names(allhrly)[n])
}
