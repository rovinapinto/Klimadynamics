##    source("~/R/SpherHam/Ruvchi.r")
##
##
#
 uvchi <- function(plm)
 {

  mplm <- NULL ; lplm <- NULL 
  ufak1 <- NULL ; ufak2 <- NULL
  plm1 <- NULL ; plm2 <- NULL

  eps <- function(l,j){ sqrt( ( l*l - m*m )/( 4 *l*l - 1 ) ) }
  
  lm <- 0
  for ( m in 0:mmm ){
   for ( l in m:mmm ){

     lm <- lm+1

     mplm[lm] <- m
     lplm[lm] <- l

     ufak1[lm] <- (l-1) * eps(l,m)
     ufak2[lm] <- -(l+2) * eps((l+1),m)

     plm1[lm] <- plm[ min((lm+1),mmjj) ]
     if ( l >= mmm ) plm1[lm] <- 0+0i

     plm2[lm] <- plm[ max((lm-1),1) ]
     if ( l == m ) plm2[lm] <- 0+0i

   }
  }

 # V component
 #
  vplm <- -ufak1*plm2 - ufak2*plm1 

 # U component
 #
  uplm <- +1i * mplm * plm

  return( list( uplm=uplm , vplm=vplm ) )
 }
