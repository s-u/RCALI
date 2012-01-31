\name{readpoly1}
\alias{readpoly1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read a  polygons file in format 1}
\description{
  Create an object of  class \code{"\link[=listpoly-class]{listpoly}"} from the values read on a
   polygons file in format 1.
}
\usage{
readpoly1(filename, delim=" ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ Name of the polygons file}
  \item{delim}{ Separator character on the file. }
}
\value{
  An  object of  class \code{"\link[=listpoly-class]{listpoly}"}, i.e a list, each component of which is an
  object of class \code{"\link[=poly-class]{poly}"}}

\seealso{ \link{readpoly2} }
\keyword{IO}
