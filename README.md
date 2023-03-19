# SphericalEllipse

<!-- badges: start -->
<!-- badges: end -->


A spherical ellipse is the locus of points on sphere for which the sum of the
great circle distances to two focal points on the sphere is a constant.  The
SphericalEllipse package provides basic functionality for producing spherical
ellipses on the surface of a (spherical) Earth.

## Installation

You can install the development version of SphericalEllipse like so:

``` r
remotes::install_github("SWotherspoon/SpehericalEllipse")
```

## Example

To plot a spherical ellipse with foci (140E, 45S) and (143E, 40S) for which the
sum of the distances from the foci to any point on the ellipse is 900km
``` r
library(SphericalEllipse)
p1 <- c(140,-45)
p2 <- c(143,-40)
ps <- sphericalEllipse(p1,p2,900000)
plot(ps,type="l",xlab="Lon",ylab="Lat")
points(rbind(p1,p2),col="red",pch=16)
```

