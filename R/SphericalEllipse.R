##' Find the centroid of n points on a sphere.
##'
##' Calculate the centroid of n points on a sphere by an iterative
##' method.
##' @title Spherical Centroid
##' @param ps a two column matrix giving the longitude and latitude of
##'   the points.
##' @param n.iters number of iterations.
##' @return A vector of length 2 giving the longitude and latitude of
##'   the centroid.
##' @export
sphericalCentroid <- function(ps,n.iters=20) {
  ps <- (pi/180)*ps
  p3 <- cbind(cos(ps[,2L])*cos(ps[,1L]),cos(ps[,2L])*sin(ps[,1L]),sin(ps[,2L]))
  pc <- c(0,0,0)
  for(k in seq_len(n.iters)) {
    pc <- colMeans(p3/sqrt(1-drop(p3%*%pc)^2))
    pc <- pc/sqrt(sum(pc^2))
  }
  (180/pi)*c(atan2(pc[2L],pc[1L]),asin(pc[3L]))
}




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
##'   of the first focal point.
##' @param p2 a vector of length 2 giving the longitude and latitude
##'   of the second focal point.
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
  p1 <- c(cos(p1[2L])*cos(p1[1L]),cos(p1[2L])*sin(p1[1L]),sin(p1[2L]))
  p2 <- c(cos(p2[2L])*cos(p2[1L]),cos(p2[2L])*sin(p2[1L]),sin(p2[2L]))

  ## Trilateration
  ex <- p1
  i <- sum(ex*p2)
  ey <- p2-i*ex
  ey <- ey/sqrt(sum(ey^2))
  j <- sum(ey*p2)
  ez <- c(ex[2L]*ey[3L]-ex[3L]*ey[2L],
          ex[3L]*ey[1L]-ex[1L]*ey[3L],
          ex[1L]*ey[2L]-ex[2L]*ey[1L])
  x <- (2-d1^2)/2
  y <- (1-d2^2-2*i*x+i^2+j^2)/(2*j)
  z <- sqrt(pmax(0,1-x^2-y^2))

  ## The two points of intersection
  p1 <- outer(x,ex)+outer(y,ey)+outer(z,ez)
  p2 <- outer(x,ex)+outer(y,ey)-outer(z,ez)

  ## Convert to lon,lat
  list(p1=cbind(atan2(p1[,2L],p1[,1L])/rad,asin(p1[,3L])/rad),
       p2=cbind(atan2(p2[,2L],p2[,1L])/rad,asin(p2[,3L])/rad))
}


##' Calculate points on a spherical ellipse on the surface of the
##' Earth.
##'
##' A spherical ellipse is the locus of points on a sphere for which
##' the sum of the great circle distances to two focal points on the
##' sphere is a constant.
##'
##' Given the latitude and longitude of two focal points `p1` and
##' `p2`, this function calculates the longitude and latitude of a
##' sequence of points of the spherical ellipse for which the sum of
##' the great circle distances to the focal points is `d`.
##'
##' The final point is identical to the first. If `d` is less than the
##' inter-focal distance, the function returns `NULL`.
##'
##' A cosine approximation is used to distribute the points evenly
##' around the ellipse.  For highly eccentric ellipses it may be
##' necessary to call [focalIntersection()] directly to calculate a
##' more even distribution of points.
##'
##' @title Spherical Ellipse
##' @param p1 a vector of length 2 giving the longitude and latitude
##'   of the first focal point.
##' @param p2 a vector of length 2 giving the longitude and latitude
##'   of the second focal point.
##' @param d the combined distance from the two focal points to a
##'   point on the ellipse.
##' @param n.pts 2n.pts-1 points will be generated.
##' @param R radius of the Earth.
##' @return a (2n.pts-1) x 2 matrix of longitudes and latitudes of points
##'   on the ellipse.
##' @examples
##' p1 <- c(140,-45)
##' p2 <- c(143,-40)
##' ps <- sphericalEllipse(p1,p2,900000)
##' plot(ps,type="l",xlab="Lon",ylab="Lat")
##' points(rbind(p1,p2),col="red",pch=16)
##' @seealso [ellipsoidalEllipse()]
##' @export
sphericalEllipse <- function(p1,p2,d,n.pts=25,R=6378137) {

  ## Radian conversion
  rad <- pi/180

  ## Convert to unit sphere
  d <- d/R

  ## Interfocal distance on unit sphere
  d0 <- acos(sin(rad*p1[2L])*sin(rad*p2[2L])+
               cos(rad*p1[2L])*cos(rad*p2[2L])*cos(rad*(p1[1L]-p2[1L])))
  if(d0<d) {
    d1 <- (d-d0*cos(seq(0,pi,length.out=n.pts)))/2
    ell <- focalIntersection(p1,p2,d1,d-d1)
    ## Ensure result is closed
    rbind(ell$p1,ell$p2[rev(1L+seq_len(n.pts-2L)),],ell$p1[1L,])
  }
}



##' Calculate points on a ellipse on the surface of the Earth.
##'
##' An ellipsoidal ellipse is the locus of points on an ellipsoid for
##' which the sum of geodesic distances to two focal points on the
##' ellipsoid is a constant.
##'
##' Given the latitude and longitude of two focal points `p1` and
##' `p2`, this function calculates the longitude and latitude of a
##' sequence of points of the Earth for which the sum of the geodesic
##' distances to the focal points is `d`.
##'
##' The final point is identical to the first. If `d` is less than the
##' inter-focal distance, the function returns `NULL`.
##'
##' This implementation constructs points on the ellipse by computing
##' the distance to the ellipse along n evenly spaced bearings from
##' the midpoint of the foci using a bisection procedure. The accuracy
##' of the procedure is approximately d/2^iters.
##'
##' This function is based upon [geosphere::distGeo()] and an
##' explanation of the parameters `a` and `f` can be found there.
##'
##' For most applications, [sphericalEllipse()] provides an adequate
##' approximation that is more efficient to compute.
##' 
##' @title Ellipsoidal Ellipse
##' @param p1 a vector of length 2 giving the longitude and latitude
##'   of the first focal point.
##' @param p2 a vector of length 2 giving the longitude and latitude
##'   of the second focal point.
##' @param d the combined distance from the two focal points to a
##'   point on the ellipse.
##' @param n.pts number of points to generate.
##' @param a major (equatorial) radius of the ellipsoid (default is
##'   for WGS84).
##' @param f ellipsoid flattening (default is for WGS84).
##' @param n.iters number of bisection iterations.
##' @return an n x 2 matrix of longitudes and latitudes of points on
##'   the ellipse.
##' @examples
##' p1 <- c(140,-45)
##' p2 <- c(143,-40)
##' ps <- ellipsoidalEllipse(p1,p2,900000)
##' plot(ps,type="l",xlab="Lon",ylab="Lat")
##' points(rbind(p1,p2),col="red",pch=16)
##' @importFrom geosphere midPoint distGeo destPoint
##' @seealso [sphericalEllipse()]
##' @export
ellipsoidalEllipse <- function(p1,p2,d,n.pts=50,
                               a=6378137, f=1/298.257223563,
                               n.iters=2+max(6,ceiling(log2(d)))) {

  p <- midPoint(p1,p2,a=a,f=f)
  d0 <- distGeo(p1,p2,a=a,f=f)
  if(d0 < d) {
    b <- seq(-180,180,length.out=n.pts+1L)[-1L]
    ## Inner boundary
    r.lwr <- rep(0,n.pts)
    p.lwr <- destPoint(p,b,r.lwr,a=a,f=f)
    d.lwr <- distGeo(p.lwr,p1,a=a,f=f)+distGeo(p.lwr,p2,a=a,f=f)
    ## Outer boundary
    r.upr <- rep((d+d0)/2,n.pts)
    p.upr <- destPoint(p,b,r.upr,a=a,f=f)
    d.upr <- distGeo(p.upr,p1,a=a,f=f)+distGeo(p.upr,p2,a=a,f=f)
    
    ## Bisection method
    for(k in seq_len(n.iters)) {
      r.mid <- (r.lwr+r.upr)/2
      p.mid <- destPoint(p,b,r.mid,a=a,f=f)
      d.mid <- distGeo(p.mid,p1,a=a,f=f)+distGeo(p.mid,p2,a=a,f=f)
      r.upr <- ifelse(d.mid >= d,r.mid,r.upr)
      r.lwr <- ifelse(d.mid <= d,r.mid,r.lwr)
    }
    
    r <- (r.lwr+r.upr)/2
    ps <- destPoint(p,b,r,a=a,f=f)
    rbind(ps,ps[1L,])
  }
}




##' Calculate points on an n-ellipse on the surface of the Earth.
##'
##' An ellipsoidal n-ellipse is the locus of points on an ellipsoid
##' for which the sum of geodesic distances to n focal points on the
##' ellipsoid is a constant.
##'
##' Given the latitude and longitude of n focal points `ps`, this
##' function calculates the longitude and latitude of a sequence of
##' points of the Earth for which the sum of the geodesic distances to
##' the focal points is `d`.
##'
##' The final point is identical to the first.
##'
##' This implementation constructs points on the n-ellipse by
##' computing the distance to the n-ellipse along n evenly spaced
##' bearings from the centroid of the foci using a bisection
##' procedure. The accuracy of the procedure is approximately
##' d/2^n.iters. The centroid is calculated assuming the Earth is
##' spherical and so this function may fail for small `d`.
##'
##' This function is based upon [geosphere::distGeo()] and an
##' explanation of the parameters `a` and `f` can be found there.
##'
##' @title Ellipsoidal N-Ellipse
##' @param ps a two column matrix giving the longitude and latitude
##'   of the focal points.
##' @param d the combined distance from the two focal points to a
##'   point on the ellipse.
##' @param n.pts minimum number of points to generate.
##' @param a major (equatorial) radius of the ellipsoid (default is
##'   for WGS84).
##' @param f ellipsoid flattening (default is for WGS84).
##' @param n.iters number of bisection iterations.
##' @return an n x 2 matrix of longitudes and latitudes of points on
##'   the ellipse.
##' @examples
##' p1 <- c(140,-45)
##' p2 <- c(143,-40)
##' p3 <- c(138,-43)
##' plot(rbind(p1,p2,p3),
##'      xlab="Lon",ylab="Lat",xlim=c(135,145),ylim=c(-47,-38),
##'      pch=16,col="red")
##' points(t(sphericalCentroid(rbind(p1,p2,p3))),pch=16,col="blue")
##' ps <- ellipsoidalNEllipse(rbind(p1,p2,p3),1400000)
##' points(ps,type="l")
##' ps <- ellipsoidalNEllipse(rbind(p1,p2,p3),1000000)
##' lines(ps,type="l")
##' ps <- ellipsoidalNEllipse(rbind(p1,p2,p3),800000)
##' lines(ps,type="l")
##' @importFrom geosphere distGeo destPoint bearing
##' @seealso [ellipsoidalEllipse()]
##' @export
ellipsoidalNEllipse <- function(ps,d,n.pts=50,
                                a=6378137, f=1/298.257223563,
                                n.iters=2+max(6,ceiling(log2(d)))) {

  n <- nrow(ps)
  p <- sphericalCentroid(ps)
  b <- sort(union(bearing(p,ps),seq(-180,180,length.out=n.pts+1L)))[-1L]

  ## Inner boundary
  r.lwr <- rep(0,length(b))
  p.lwr <- destPoint(p,b,r.lwr,a=a,f=f)
  d.lwr <- 0
  for(k in seq_len(n)) d.lwr <- d.lwr + distGeo(p.lwr,ps[k,],a=a,f=f)

  if(all(d.lwr <= d)) { 

    ## Outer boundary
    r.upr <- rep(2*d,length(b))
    p.upr <- destPoint(p,b,r.upr,a=a,f=f)
    d.upr <- 0
    for(k in seq_len(n)) d.upr <- d.upr + distGeo(p.upr,ps[k,],a=a,f=f)
    
    ## Bisection method
    for(k in seq_len(n.iters)) {
      r.mid <- (r.lwr+r.upr)/2
      p.mid <- destPoint(p,b,r.mid,a=a,f=f)
      d.mid <- 0
      for(k in seq_len(n)) d.mid <- d.mid + distGeo(p.mid,ps[k,],a=a,f=f)
      r.upr <- ifelse(d.mid >= d,r.mid,r.upr)
      r.lwr <- ifelse(d.mid <= d,r.mid,r.lwr)
    }
    
    r <- (r.lwr+r.upr)/2
    ps <- destPoint(p,b,r,a=a,f=f)
    rbind(ps,ps[1L,])
  }
}





