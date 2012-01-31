\name{crpoly}
\alias{crpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create un object of class 'poly' by clicking on points }
\description{
  Create un object of class \code{"\link[=poly-class]{poly}"}   by clicking on
  a graphic: the locations of the clicks will be
  the vertices of  the polygon.
}
\usage{
crpoly()
}
%- maybe also 'usage' for other objects documented here.
\details{
The system prompts the user to enter points on the current graphic
 using the mouse or other pointing device. The points are joined
on the screen with the current line symbol. A polygon of the points
entered is drawn on the current graphics device and returned as a 'poly'
object. 
}
\value{
  An object of class \code{"\link[=poly-class]{poly}"} i.e a two-columns matrix with the coordinates.
}
\section{Warning }{A plot should be drawn on the current graphics device
  before.}
\seealso{ \code{\link{as.poly}}, \code{getpoly} of the package \code{splancs} }
\examples{
\dontrun{plot(x=c(1,10), y=c(1,10), type='n')}
\dontrun{a<-crpoly()}
#    Enter points with button 1
#    Finish with button 2

}
\keyword{hplot }
\keyword{ data}
