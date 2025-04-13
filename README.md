# SphericalEllipse

<!-- badges: start -->
<!-- badges: end -->


A spherical or ellipsoidal ellipse is the locus of points on a sphere or ellipsoid for which the sum
of the geodesic distances to two focal points on the surface is a constant.  The SphericalEllipse
package provides basic functionality for producing ellipses on the surface of a spherical or
ellipsoidal Earth.

## Installation

The current version of Combin8R can be installed from GitHub using the remotes package.
``` r
remotes::install_github("SWotherspoon/SphericalEllipse")
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

