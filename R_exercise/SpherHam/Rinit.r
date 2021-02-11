## source("~/R/KlimDyn/SpherHam/Rinit.r")
##
##
##  Rinit.r loads functions Rplmhin and Rplmbck
##  and calculates:
##   Gaussian grid
##   Associtated Legendre polynomials (Ralp)
##   Eigenvalue of laplace operator on sphere
##   Weights for integration over Gaussian latitudes
##   Indices l and m of Ralp
##   Deviation of associtated Legendre polynomials
##
 
  source("~/R/KlimDyn/SpherHam/Rplmhin.r")
  source("~/R/KlimDyn/SpherHam/Rplmbck.r")

  source("~/R/KlimDyn/SpherHam/Ruvpsi.r")
  source("~/R/KlimDyn/SpherHam/Ruvchi.r")

  source("~/R/KlimDyn/SpherHam/gauleg.r")
  source("~/R/KlimDyn/SpherHam/lgndre.r")

## STARTN gives PLM.init 
##
  if ( truncation == 42 ){
    mmm <<- 42 ; mmjj <<- (mmm+1)*(mmm+2)/2
    nlong <<- 128 ; nlatg <<- 64
  }else if ( truncation == 21 ){
    mmm <<- 21 ; mmjj <<- (mmm+1)*(mmm+2)/2
    nlong <<- 64 ; nlatg <<- 32
  }else if ( truncation == 63){
    mmm <<- 63 ; mmjj <<- (mmm+1)*(mmm+2)/2
    nlong <<- 192 ; nlatg <<- 96
  }else{ print("ERROR: Choose T21, T42, or T63 only!") }

## Calculate abscissas and weights of the Gauss-Legendre N-points
## quadrature formula (latitudes)
#
  GL <- gauleg(nlatg)
  slatg  <<- GL$x
  gauwch <<- GL$w
  latg   <<- asin(slatg)*180/pi     # latitudes
  clatg  <- cos(latg*pi/180)        # cosine of latitudes

## Calculate longitudes
#
  long <<- seq(0,(360-360/nlong),by=360/nlong)

## Compute index matrix for total wavenumber l and
## zonal wavenumer m
#
  iml <<- array(0,c((mmm+1),(mmm+1)))
  lm <- 0
  lplm <- NULL
  lplm1 <- NULL
  for ( m in 0:mmm ){
   for ( l in m:mmm ){
     lm <- lm+1
     iml[(m+1),(l+1)] <- lm
     lplm[lm] <- l
     lplm1[lm] <- max(l+1,0)
   }
  }

## Complute associated Legendre polynomials using IMSL routine LGNDRE.f,
## which was adapted to R
#
  LP <- lgndre(mmm,mmjj,latg,slatg,clatg)
  Ralp <<- sqrt(2) * LP$alp
