\name{export.listpoly}
\alias{export.listpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a polygons file  in format 1. }
\description{
  From an object of class \code{"\link[=listpoly-class]{listpoly}"},
  create a polygons file for input to the software \code{CaliFloPP}
  or to the function
  \link{califlopp}. The format of this file corresponds to the
  format \code{1} of them, i.e you should set the param \code{input} to
  \code{1} when used as input to  \link{califlopp}.
}
\usage{
export.listpoly(x, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of class \code{"\link[=listpoly-class]{listpoly}"}}
  \item{filename}{ Name of the file to be created.}
}
\value{
  None.
}
\keyword{IO}
