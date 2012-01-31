\name{crlistpoly}
\alias{crlistpoly}
\title{ Create an object of class 'listpoly' from objects of class 'poly' }
\description{
 Create an object of class \code{"\link[=listpoly-class]{listpoly}"} from objects of class \code{"\link[=poly-class]{poly}"}
}
\usage{
crlistpoly(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ objetcs of class \code{"\link[=poly-class]{poly}"} }
}

\value{
  An object of class \code{"\link[=listpoly-class]{listpoly}"} : a list
  where each component is a 'poly' object (see \code{\link{poly-class}}).
}
\seealso{  \code{\link{poly-class}} }
\examples{
# A triangle:
a <- as.poly(matrix(c(2,2,2,3,3,3), ncol=2, byrow=TRUE))
# A square:
b <- as.poly(matrix(c(2.5,2,2.5,2.5,3,2.5,3,2), ncol=2, byrow=TRUE))
# The both:
 z <- crlistpoly(a,b)
}
\keyword{ data }
