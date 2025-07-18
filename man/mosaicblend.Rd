% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile_functions.R
\name{mosaicblend}
\alias{mosaicblend}
\title{Mosaics a list of overlapping SpatRasters blending overlap areas}
\usage{
mosaicblend(rlist)
}
\arguments{
\item{rlist}{a list of SpatRasters}
}
\description{
Mosaics a list of overlapping SpatRasters blending
the areas of overlap using a distance weighting to eliminate tiling effects
}
\details{
If rlist contains SpatRasters that are not overlapping the conventional terra::moasic function is used.
If rlist contains SpatRasters that do overlap, they should comprise a list of adjacent rasters in a single row or column.
Overlapping cells calulated using weighted distances to reduce tile effects
Function derived from: https://github.com/ilyamaclean/microclimf/tree/main
Called by: \code{blend_tile_lists}
}
