% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile_functions.R
\name{blend_tile_lists}
\alias{blend_tile_lists}
\title{Blend tile list}
\usage{
blend_tile_lists(input_list)
}
\arguments{
\item{input_list}{\itemize{
\item list of overlapping spatraster tiles
}}
}
\value{
spatRaster of merged tiles
}
\description{
Blend tile list
}
\details{
Call function \code{mosaicblend} to within then across columns.
}
