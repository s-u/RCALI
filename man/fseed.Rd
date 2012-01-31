\name{fseed}
\alias{fseed}

\title{ Individual seed dispersion function}
\description{
  From a vector of distances, calculate the individual seed dispersion
   function used in  \emph{Genesys}
}
\usage{
fseed(point)
}
\arguments{
  \item{point}{ scalar or vector of length 2, whose first element is
    a distance expressed in meter.}
}
\value{
  The calculated dispersions.
}
\seealso{ \link{fpollen} }
\references{ Colbach, N. and Clermont-Dauphin, C. and Meynard, J.M. \emph{Genesys: a model of the
influence of cropping system on gene escape from herbicide tolerant rapeseed crops
to rape volunteers. i. temporal evolution of a population of rapeseed volunteers in a
field}. Agriculture, Ecosystems and Environnement, 83:235-253,
2001.
}

\examples{

distance = seq(1,1.5, by=0.05)
a=matrix(distance, ncol=1)
b= apply(a,1,fseed)
par(pty="s")
plot(x=distance, y =b)
}
\keyword{manip}
