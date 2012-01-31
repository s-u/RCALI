\name{as.poly}
\alias{as.poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create an object  of class 'poly' }
\description{
 Create an object of class \code{"\link[=poly-class]{poly}"} from two vectors or from a matrix.
}
\usage{
as.poly(x, y=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of x-coordinates or two-columns matrix  }
  \item{y}{ vector of y-coordinates when 'x' is a vector. }
}
\value{
  An object of class \code{"\link[=poly-class]{poly}"}, i.e a two-columns matrix  labelled "xcoord", "ycoord".
}
\seealso{
  \code{\link{crpoly}}
}
\examples{
# A triangle
a <- as.poly(matrix(c(2,2,2,3,3,3), ncol=2, byrow=TRUE))
}

\keyword{manip}

