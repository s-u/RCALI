\name{generPoly}
\alias{generPoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Generate a regular grid of polygons }
\description{
  Generate a regular grid of polygons and optionally write
  them on a    file for input to \code{\link{califlopp}}
}
\usage{
generPoly(step = 5, np = 10, file = "data", plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{step}{ Step of grid, in meter. }
  \item{np}{ Number of polygons on each axis. }
  \item{file}{ Pathname of the polygons file to be created. If
  NULL, no file created.}
  \item{plot}{ If TRUE, plot of the  polygons. }
}
\details{
  When \code{file} is not NULL, a polygons file is created. It
  is in the  format 2 of \code{\link{califlopp}}
  and the values separator is tabulate.
}
\note{
  The polygons are numbered from 1, from
  the left to the right, then from bottom to top.
  }
\value{
  An object of class \code{"\link[=listpoly-class]{listpoly}"}
}
\seealso{ \code{"\link[=listpoly-class]{listpoly-class}"}, \code{\link{plot.listpoly}}}
\examples{
a <- generPoly(np=3, file=NULL)
}
\keyword{ data }
\keyword{ hplot}
