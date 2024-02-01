\name{era5data}
\alias{era5data}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("era5data")}
\format{
  The format is:
List of 3
 $ tme      : POSIXlt[1:744], format: "2018-05-01 00:00:00" "2018-05-01 01:00:00" "2018-05-01 02:00:00" ...
 $ climarray:List of 9
  ..$ temp     : num [1:3, 1:3, 1:744] 8.86 8.9 9.23 7.6 7.62 ...
  ..$ relhum   : num [1:3, 1:3, 1:744] 73.6 72.6 70.5 77.3 76.3 ...
  ..$ pres     : num [1:3, 1:3, 1:744] 101 101 101 101 101 ...
  ..$ swrad    : num [1:3, 1:3, 1:744] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ difrad   : num [1:3, 1:3, 1:744] -6.47e-14 -6.47e-14 -6.47e-14 -6.47e-14 -6.47e-14 ...
  ..$ lwrad    : num [1:3, 1:3, 1:744] 358 359 362 348 349 ...
  ..$ windspeed: num [1:3, 1:3, 1:744] 4.66 4.6 4.73 3.71 3.7 ...
  ..$ winddir  : num [1:3, 1:3, 1:744] 265 265 267 266 268 ...
  ..$ prech    : num [1:3, 1:3, 1:744] 0.003416 0.001677 0.000497 0.003688 0.001107 ...
 $ dtmc     :Formal class 'PackedSpatRaster' [package "terra"] with 3 slots
  .. ..@ definition: chr "rast(ncols=3, nrows=3, nlyrs=1, xmin=130022.665503526, xmax=191016.344812007, ymin=-973.379523697156, ymax=6002"| __truncated__
  .. ..@ values    : num [1:9, 1] -28.25 -5.06 66.05 43.2 54.77 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "dtm_z10"
  .. ..@ attributes: list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(era5data)
## maybe str(era5data) ; plot(era5data) ...
}
\keyword{datasets}
