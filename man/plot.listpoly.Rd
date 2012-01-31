\name{plot.listpoly}
\alias{plot.listpoly}
\title{ Plot of an object of class 'listpoly' }
\description{
  Plot all the polygons of an  object of class \code{"\link[=listpoly-class]{listpoly}"}, optionnaly with
  colors.
}
\usage{
\method{plot}{listpoly}(x, add = F, color = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of class \code{"\link[=listpoly-class]{listpoly}"}.}
  \item{add}{ TRUE, to add to a previous plot, for example after a
    zoom. }
  \item{color}{ TRUE, to fill each polygon in a different color.}
  \item{\dots}{ List of arguments passed as-is to \link{plot.poly}}
}
\value{
  None.
}
\note{
  To make a zoom, type in: zoom() (see 'splancs') then use this function
  with argument \dQuote{add=TRUE}.}

\keyword{hplot}
