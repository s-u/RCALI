\name{generVois}
\alias{generVois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate the neighbors of each polygon of a regular grid }
\description{
  For  each polygon of a regular grid, generate the indexes
  of its neighbour polygons }
\usage{
generVois(np = 10, nvois = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{np}{ Number of polygons on each axis of the grid. }
  \item{nvois}{ Number of polygons required in the neighbourhood,
    on each direction, for each polygon}
}
\details{
  For each polygon P of a regular grid made of \code{np}
  polygons along each axis, the indexes of its neighbours
  are generated: the neighbours are the polygons included
  in the square  made of
  \code{nvois} polygons 
  below, above, at the
  right and left sides of P,  in the limit
  of the grid bounds.
}
\value{
  A two-columns matrix, of all pairs of neighbourgs.
}
\seealso{ \code{\link{generPoly}},
  \code{\link{califlopp}}  }
\examples{
generVois(np=3)
}
\keyword{ data }
