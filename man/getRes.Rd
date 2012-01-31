\name{getRes}
\alias{getRes}
\title{ Read a result-file of 'califlopp' }
\description{
  Return the results stored on a file created by \code{califlopp}.
}
\usage{getRes(ficres)}
\arguments{
  \item{ficres}{ Pathname of a result-file created by \code{\link{califlopp}}. }
  }
\value{
  A data.frame with as many rows as pairs of polygons on the file.
  The columns are:
  \item{poly1, poly2 }{Identifiers of the  polygons.}
  \item{mean.flow/area }{Integrated flow divided by the area of the second polygon.
    If several dispersion functions have been studied, as many columns
    as
    functions. The columns labels are then
 \code{mean.flow.f1/area}, \code{mean.flow.f2/area}, etc\ldots}
  \item{area1, area2 }{Areas of the polygons in squared meters.}
  \item{mean.flow, conf.int.lower, conf.int.upper, abs.err, n.eval}{
when the method is "cubature", only: integrated flow, lower and upper bounds
of the confidence interval, absolute error, number of evaluations.
If several dispersion functions have been studied, the columns labels
are suffixed with ".f1", "f2", etc\ldots}
  \item{mean.flow, std}{
    when the method is "grid", only: mean of the integrated flow,
    standard
    deviation.
    If several dispersion functions have been studied, the columns labels
    are suffixed with ".f1", "f2", etc\ldots}
  }

\note{ This function  works  when \code{RCALI}      has been configured with
  \code{OUTPUT_FILE_FORMAT = LIGHT} only (this is the default; see the
  file \code{src/caliconfig.h} to be sure your configuration is compatible).
  
Details about the returned values can be found in the Reference
  Manual.
  }
  \seealso{  \code{\link{califlopp}}}
\keyword{ data }
