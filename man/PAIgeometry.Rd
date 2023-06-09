% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microclimatemodel.R
\name{PAIgeometry}
\alias{PAIgeometry}
\title{Generates plant area index profile}
\usage{
PAIgeometry(PAI, skew, spread, n = 100)
}
\arguments{
\item{PAI}{Total plant area index of canopy}

\item{skew}{number between 0 and 10 indicating the degree of skew towards top of
canopy in canopy foliage (see details)}

\item{spread}{positive non-zero number less than 100 indicating the degree of spread in
canopy foliage (see details)}

\item{n}{Number of plant area index values to generate. Default: 100 (see details)}
}
\value{
a vector of length `n` of plant area index values the sum of which equals `PAI`
}
\description{
Generates a vector of length `n` of plausible plant area index values
}
\details{
when specifying `skew`, lower numbers indicate greater skew towards top of
canopy (5 = symmetrical). In specifying `spread` a value of one indicates almost
all the foliage in concentrated in one canopy layer, whereas a value of 100 indicates
completely evenly spread.
}
\examples{
pai <- PAIgeometry(3, 7, 70)
plot(pai, type = "l")
}
