##' Find the points on a unit sphere that are specified great circle
##' distances from two focal points.
##'
##' Given two focal points `p1` and `p2` on the surface of the unit
##' sphere, this function calculates the two points that are distance
##' `d1` from `p1` and distance `d2` from `p2`. The vectors of
##' distances `d1` and `d2` must be the same length, and two alternate
##' points are given for each pairs of great circle distances.
##' 
##' @title Focal Intersection
##' @param p1 a vector of length 2 giving the longitude and latitude
##'   of the first focal point
##' @param p2 a vector of length 2 giving the longitude and latitude
##'   of the second focal point
##' @param d1 a vector of great circle distances to the first focal
##'   point.
##' @param d2 a vector of great circle distances to the second focal
##'   point.
##' @return a list of two, two-column matrices representing the
##'   longitude and latitude of the two alternate points.
##' @export
focalIntersection <- function(p1,p2,d1,d2) {

  ## Convert to radians
  rad <- pi/180
  p1 <- rad*p1
  p2 <- rad*p2
  
  ## Convert gc distance to chord length on unit sphere
  d1 <- sqrt((cos(d1)-1)^2+sin(d1)^2)
  d2 <- sqrt((cos(d2)-1)^2+sin(d2)^2)

  ## Cartesian coords on the unit sphere
  p1 <- c(cos(p1[2])*cos(p1[1]),cos(p1[2])*sin(p1[1]),sin(p1[2]))
  p2 <- c(cos(p2[2])*cos(p2[1]),cos(p2[2])*sin(p2[1]),sin(p2[2]))

  ## Trilateration
  ex <- p1
  i <- sum(ex*p2)
  ey <- p2-i*ex
  ey <- ey/sqrt(sum(ey^2))
  j <- sum(ey*p2)
  ez <- c(ex[2]*ey[3]-ex[3]*ey[2],ex[3]*ey[1]-ex[1]*ey[3], ex[1]*ey[2]-ex[2]*ey[1])
  x <- (2-d1^2)/2
  y <- (1-d2^2-2*i*x+i^2+j^2)/(2*j)
  z <- sqrt(pmax(0,1-x^2-y^2))

  ## The two points of intersection
  p1 <- outer(x,ex)+outer(y,ey)+outer(z,ez)
  p2 <- outer(x,ex)+outer(y,ey)-outer(z,ez)

  ## Convert to lon,lat
  list(p1=cbind(atan2(p1[,2],p1[,1])/rad,asin(p1[,3])/rad),
       p2=cbind(atan2(p2[,2],p2[,1])/rad,asin(p2[,3])/rad))
}
  

##' Calculate points on a spherical ellipse on the surface of the
##' Earth.
##'
##' A spherical ellipse is the locus of points on sphere for which the
##' sum of the great circle distances to two focal points on the
##' sphere is a constant.
##' 
##' Given the latitude and longitude of two focal points `p1` and
##' `p2`, this function calculates the longitude and latitude of a
##' sequence of points of the spherical ellipse for which the sum of
##' the great circle distances to the focal points is `d`.
##'
##' The final point is identical to the first. A cosine approximation
##' is used to distribute the points evenly around the ellipse.  For
##' highly eccentric ellipses it may be necessary to call
##' [focalIntersection()] directly to calculate a more even
##' distribution of points.
##' 
##' @title Spherical Ellipse
##' @param p1 a vector of length 2 giving the longitude and latitude
##'   of the first focal point
##' @param p2 a vector of length 2 giving the longitude and latitude
##'   of the second focal point
##' @param d the combined distance from the two focal points to a
##'   point on the ellipse.  point.
##' @param n 2n-1 points will be generated.
##' @param R radius of the Earth
##' @examples
##' p1 <- c(140,-45)
##' p2 <- c(143,-40)
##' ps <- sphericalEllipse(p1,p2,900000)
##' plot(ps,type="l",xlab="Lon",ylab="Lat")
##' points(rbind(p1,p2),col="red",pch=16) 
##'
##' @return a (2n-1) x 2 matrix of longitudes and latitudes of points
##'   on the ellipse.
##' @export
sphericalEllipse <- function(p1,p2,d,n=25,R=6378137) {

  ## Radian conversion
  rad <- pi/180

  ## Convert to unit sphere
  d <- d/R
  
  ## Interfocal distance on unit sphere
  f <- acos(sin(rad*p1[2])*sin(rad*p2[2])+
              cos(rad*p1[2])*cos(rad*p2[2])*cos(rad*(p1[1]-p2[1])))
  if(f<d) {
    d1 <- (d-f*cos(seq(0,pi,length.out=n)))/2
    ell <- focalIntersection(p1,p2,d1,d-d1)
    rbind(ell$p1,ell$p2[rev(seq_len(n-1)),])
  }
}
