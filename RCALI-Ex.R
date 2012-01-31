pkgname <- "RCALI"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('RCALI')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("as.poly")
### * as.poly

flush(stderr()); flush(stdout())

### Name: as.poly
### Title: Create an object of class 'poly'
### Aliases: as.poly
### Keywords: manip

### ** Examples

# A triangle
a <- as.poly(matrix(c(2,2,2,3,3,3), ncol=2, byrow=TRUE))



cleanEx()
nameEx("califlopp")
### * califlopp

flush(stderr()); flush(stdout())

### Name: califlopp
### Title: Calculation of the Integrated Flow of Particles between Polygons
### Aliases: califlopp
### Keywords: optimize

### ** Examples

# Grid method with compiled constant and seed dispersion functions:
param <- list(method="grid",  grid=list(step=c(50,50)))
## Not run: califlopp("MyPolygonsFile",dispf=c(3,1), param=param)

# Cubature method with a R dispersion function:
param <- list( output=1, input=2, dz=0, dp=100, tz=0)
## Not run: califlopp("MyPolygonsFile", dispf=fpollen, param=param)



cleanEx()
nameEx("crlistpoly")
### * crlistpoly

flush(stderr()); flush(stdout())

### Name: crlistpoly
### Title: Create an object of class 'listpoly' from objects of class
###   'poly'
### Aliases: crlistpoly
### Keywords: data

### ** Examples

# A triangle:
a <- as.poly(matrix(c(2,2,2,3,3,3), ncol=2, byrow=TRUE))
# A square:
b <- as.poly(matrix(c(2.5,2,2.5,2.5,3,2.5,3,2), ncol=2, byrow=TRUE))
# The both:
 z <- crlistpoly(a,b)



cleanEx()
nameEx("crpoly")
### * crpoly

flush(stderr()); flush(stdout())

### Name: crpoly
### Title: Create un object of class 'poly' by clicking on points
### Aliases: crpoly
### Keywords: hplot data

### ** Examples

## Not run: plot(x=c(1,10), y=c(1,10), type='n')
## Not run: a<-crpoly()
#    Enter points with button 1
#    Finish with button 2




cleanEx()
nameEx("fpollen")
### * fpollen

flush(stderr()); flush(stdout())

### Name: fpollen
### Title: Individual pollen dispersion function
### Aliases: fpollen
### Keywords: manip

### ** Examples

distance = seq(1,1.5, by=0.05)
a=matrix(distance, ncol=1)
b= apply(a,1,fpollen)
par(pty="s")
plot(x=distance, y =b)
lines(x=distance, y = apply(a,1,fseed))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fseed")
### * fseed

flush(stderr()); flush(stdout())

### Name: fseed
### Title: Individual seed dispersion function
### Aliases: fseed
### Keywords: manip

### ** Examples


distance = seq(1,1.5, by=0.05)
a=matrix(distance, ncol=1)
b= apply(a,1,fseed)
par(pty="s")
plot(x=distance, y =b)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("generPoly")
### * generPoly

flush(stderr()); flush(stdout())

### Name: generPoly
### Title: Generate a regular grid of polygons
### Aliases: generPoly
### Keywords: data hplot

### ** Examples

a <- generPoly(np=3, file=NULL)



cleanEx()
nameEx("generVois")
### * generVois

flush(stderr()); flush(stdout())

### Name: generVois
### Title: Generate the neighbors of each polygon of a regular grid
### Aliases: generVois
### Keywords: data

### ** Examples

generVois(np=3)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
