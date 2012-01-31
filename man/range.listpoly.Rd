\name{range.listpoly}
\alias{range.listpoly}
\title{ Range of the coordinates of an object
  of class 'listpoly' }
\description{
  Return the minimum and maximum coordinates over all the polygones
  of an object of class \code{"\link[=listpoly-class]{listpoly}"}.
}
\usage{
\method{range}{listpoly}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object  of class \code{"\link[=listpoly-class]{listpoly}"} }
  \item{\dots}{ arguments passed as-it to the R function \code{range} }
}
\value{
  Matrix ; its two columns are \dQuote{xrange} and \dQuote{yrange} and
  its two lines  are \dQuote{lower} and \dQuote{upper}.
}
\keyword{manip}
