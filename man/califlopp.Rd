\name{califlopp}
\alias{califlopp}
\title{ Calculation of the Integrated Flow of Particles
    between Polygons}
\description{
  Calculation of the flow of particles between polygons by two integration methods: integration by a cubature method and
integration on a grid of points.
}
\usage{
califlopp(file, dispf=c(1,2), param=NULL, resfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ Pathname of the polygons-file. See details. }
  \item{dispf}{ The required dispersion functions. Vector of integers  or vector of R functions. The maximum length of this vector is 5.
    \itemize{
      \item{    If vector of integers, the dispersion functions are
    then compiled. Five are provided.  To modify them, see
    Details.
    By default, 
    1 is for dispersal of oilseed rape pollen, 2 for dispersal
    of oilseed rape seed (dispersals of oilseed rape are the ones
    defined in \emph{GeneSys} - see References),
    3 for the constant function, 4 for an
    anisotropic version of the dispersal of yellow rust of wheat defined
    in \emph{Soubeyrand and all}, 5 for a discontinuous function.\cr
    These functions are viewable in the C file \code{src/functions.cc}.
    To modify them, see Details.}
   \item{    If vector of functions, the dispersion functions are coded in R
    (more time consuming than compiled version). Two R dispersal functions are provided,  \code{\link{fpollen}} and
    \code{\link{fseed}}, the functions used in Genesys.
    To specify your own function, see details.}}}

  \item{param}{ Optional list of parameters. Valid components are
    \code{input}, \code{output},
    \code{verbose}, \code{warn.poly}, \code{warn.conv},\code{delim},
    \code{poly}, \code{send.and.receive}, \code{method}, \code{dz}, \code{dp}.
    In addition, when \code{method} is \dQuote{\code{cub}}: \code{maxpts},
    \code{reler}, \code{abser}, \code{tz}.
    When \code{method} is \dQuote{\code{grid}}: \code{seed},
    \code{step}, \code{nr}.
    See details. }
  \item{resfile}{  Optional pathname of a result-file. When set,
    the results are written on it. This file can be read
     by using  function
    \code{\link{getRes}} or \code{\link{read.table}}. See details, as to the content of the file.}
}
\value{
  Nothing. To store the results, set the argument \code{resfile},
  then use the function
  \code{\link{getRes}} or \code{\link{read.table}}, with option
      \code{skip=1}}

\details{
  \bold{The polygons-file}
  
  The coordinates of the polygons should be provided in an ASCII file,
  denoted here "polygons-file". The unit is the meter. The vertices should
  be ordered clockwise. The polygons can be closed
  or not, but  without holes. The first line contains the number of polygons.
  The following lines depend of the \code{input}
  parameter:
  \itemize{
  \item{\code{input=1}}{ Two lines per polygon: on the first one,
    an identifier (a positive
    integer),   followed by the x-coordinates,
    on the second one, the same identifier
    followed by the y-coordinates. The function
    \link{export.listpoly} generates such a file from R structures}
  \item{\code{input=2}}{ Three lines per polygon: on the first one,
    an identifier (a positive
    integer),   followed by a name for the polygon and
    by the number of its vertices, on the second one, the x-coordinates,
    and on the third one, the y-coordinates.}}

  \bold{The individual dispersion functions}
  
  The individual dispersion functions can  either be compiled
  or R functions.
\itemize{

\item{
  \emph{Compiled function:}}{ Five compiled dispersal functions are
  provided (see argument \code{dispf}). To  replace them by yours,
  you have to download the source of \code{RCALI},
  modify and compile it.\cr 
  
  Suppose that you have download the tar-archive
  in the directory \code{MyDir}.
     The steps to customize the dispersion functions are:

     1/ Replace one or several functions in
     \code{MyDir/RCALI/src/functions.cc} by yours:
     The dispersion function has one argument,
    the current point, \code{p}, of class \code{Point}.
    You can use \code{p.getX()}
    and \code{p.getY()} to get the coordinates of the current point
    (in meters*\code{SCALE}, where \code{SCALE} is the  rescaling
    parameter defined in the file \code{src/calicinfig.h}),
    \code{p.dist0()}, the distance of \code{p}
    from the origin (in meters*\code{SCALE})
    and \code{p.angle0()}, the angle (in degrees, in \eqn{[-\pi, +\pi]}{[-PI, +PI]}) between the line  (0,\code{p})
    and the horizontal   line.

    2/ Create a directory \code{MyDir/RCALI/libs},
    place you in \code{MyDir/RCALI/src} and  type in:
    \preformatted{R CMD SHLIB -o ../libs/RCALI.so  *.cc}
    to create the compiled shared library.

3/ To use in a R-session:
    \preformatted{
      source("MyDir/RCALI/R/sourceDir.R")
      sourceDir("MyDir/RCALI/R")
      dyn.load("MyDir/RCALI/libs/RCALI.so")
  }

  The help-files are viewable by opening in a browser
  \code{MyDir/RCALI/inst/doc/html/00Index.html}
}
You can also build  the tar.gz file again, after modifications,
 and install it as a library
by using the standard R commands \code{R CMD build} and
\code{R CMD INSTALL --html}. 
  \item{\emph{R function:}}{  Two R dispersal functions are provided,  \code{\link{fpollen}} and
    \code{\link{fseed}}, the functions used in Genesys.\cr
    You can define your own R dispersal function:
    it should have one vector
    argument,  the localization of the current point, \code{p}.
    The first element of this vector is the distance of \code{p} from the origin
  (in meters) and the second one is the angle (in degrees, in \eqn{[-\pi, +\pi]}{[-PI, +PI]}) between the line  (0,\code{p})
  and the horizontal   line (i.e, stating \code{x} and \code{y} are the coordinates
  of \code{p}, the angle is
  \eqn{atan2(y,x)*\frac{180}{\pi}}{atan2(y,x)*180/PI})
}}
  \bold{The parameters}

  The argument \code{param} is a list which valid components are:
  \itemize{
  \item{\code{input}}{ format of the polygons-file. 1 or 2 (see above).
  Default 2}
  \item{\code{output}}{ output required  on the screen:
 0 nothing, 1: all  results, 2:  progression numbers,
 3: the integrated flows and their means per squared meter.
Default 1}
    \item{\code{verbose}}{ TRUE, if output is required about polygons convexity
      and landscape translation. Default FALSE}
        \item{\code{warn.poly}}{ TRUE, if output is required about polygons simplification. Default FALSE}
        \item{\code{warn.conv}}{ TRUE, if output is required when
	  cubature convergence is not reached. Default TRUE}
 \item{\code{delim}}{ separator character between values in the polygons-file. Default: tabulate}
\item{\code{send.and.receive}}{ TRUE, if results are required from
  sending polygons to target polygons and from target polygons to
  sending polygons (case of anisotropic functions). Default FALSE}
 \item{\code{poly}}{ required pairs of polygons.
   List of vectors of length 2, or two-columns  matrix.
   If only one pair is required, it may be a vector of length 2.
 Default: all pairs of polygons.}
   \item{\code{method}}{ string equal to \code{cub} for cubature method, \code{grid} for the grid
     method. Default: \code{cub}}
   \item{\code{dz}}{ integer vector, whose length is greater or equal to the number of
     required dispersion functions. \code{dz[i]} is the distance in meters
     beyond which the \code{i}st dispersion function is considered as
     nul. Default in a standard configuration: 0,21,0,1000,0
   for functions number 1 to 5, respectively.}
 \item{\code{dp}}{ integer vector, whose length is greater or equal
    to the number of
     required dispersion functions. \code{dp[i]} is the distance in meters
     beyond which the \code{i}st dispersion function is calculated between
     centroids only. Default in a standard configuration: 100, 0, 0,
     500, 0
   for functions number 1 to 5, respectively.}}
   
   In addition, when \code{method}  is \code{cub}:
   \itemize{
     \item{\code{maxpts}}{ maximal
     number of evaluation points required for each function.
     Vector of length equal to the number of required functions.
   Default in a standard configuration: 100000 }
    \item{\code{reler}}{ relative error required for each function.
      Vector of length equal to the number of required functions.
      Should be positive when method is cubature.
    Default in a standard configuration: 1.0e-3 }
    \item{\code{abser}}{ absolute error required for each function.
      Vector of length equal to the number of required functions.
      Should be positive when method is cubature.
      Default in a standard configuration: 1.0e-3}
    \item{\code{tz}}{ integer vector, whose length is greater
      or equal       to the number of
     required dispersion functions. Mode of triangulation for the
     cubature method.
     \code{tz[i]} should be 1, if, for the \code{i}st dispersion
     function,
     triangulation from (0,0) has to be done
when  (0,0) is included in the integration area and, 0 if not.
1 is recommended  when the dispersion function is very "sharp"
at the origin. Default in a standard configuration: 0,1,0,0,0
   for functions number 1 to 5, respectively.}}
    
    When  \code{method} is \code{grid}:
\itemize{
    \item{\code{seed}}{ seed of the random generator.}
    \item{\code{step}}{ step of the grid on the x-axis and on the y-axis
      in meter. Vector of length 2.}
    \item{\code{nr}}{ maximal number of replications or grids.}}

    \bold{The result-file}

    When the argument \code{resfile} is set, a file is created.
    On this file, the values are separated by tabulates.

    
    Its contains, when the method is \code{cub},
    \itemize{
      \item{ on the first line:}{ "npoly:",	"input-file:",
      "nfunc:",	"method:", each of these identifiers
      followed by the actual values.
    }
    \item{ on each of the following lines, the results for a couple of
      polygons:}{ the identifiers of both polygons;
	the integrated flow divided by the area of the second
	polygon, for
	each dispersal function;
	the areas of both polygons; then, for
	each dispersal function,
	the integrated flow, the lower and upper bounds
	of the confidence interval, the absolute error, and the
	number of evaluations.}}

      Its contains, when the method is \code{grid},
      \itemize{
	\item{ on the first line:}{ "npoly:",	"input-file:",
      "nfunc:",	"method:", "stepx:", "stepy:", each of these identifiers
      followed by the actual values.
    }
    \item{ on each of the following lines, the results for a couple of
      polygons:}{ the identifiers of both polygons;
	the integrated flow divided by the area of the second
	polygon, for
	each dispersal function;
	the areas of both polygons; then, for
	each dispersal function,
	the integrated flow, and the standard deviation.
      }}

      This file can be read in a R-session by using the function
      \code{\link{getRes}} or \code{\link{read.table}}, with option
      \code{skip=1}.
      
  }
\section{Side effect}{This function creates a temporary file to store the
  parameters, usually in the directory \code{tmp} of the user. This file
  is destroyed at the end of execution.} 
\references{  - The CaliFloPP software: \url{http://genome.jouy.inra.fr/logiciels/califlopp} \cr
- Main reference paper: A. Bouvier, K. Kieu, K. Adamczyk, and H. Monod. Computation of
  integrated flow of particles between polygons. Environmental Modelling
  & Software, 24:843--849, 2009.\cr
- N. Colbach, and all.Genesys: a model of the influence of cropping
system on gene escape from herbicide tolerant rapeseed crops to rape
volunteers. Agriculture, Ecosystems and Environnement, 83:235--270, 2001. }
\author{ A. Bouvier }
  \seealso{  \code{\link{getRes}}}
  \examples{
# Grid method with compiled constant and seed dispersion functions:
param <- list(method="grid",  grid=list(step=c(50,50)))
\dontrun{califlopp("MyPolygonsFile",dispf=c(3,1), param=param)}

# Cubature method with a R dispersion function:
param <- list( output=1, input=2, dz=0, dp=100, tz=0)
\dontrun{califlopp("MyPolygonsFile", dispf=fpollen, param=param)}
}
\keyword{ optimize }
