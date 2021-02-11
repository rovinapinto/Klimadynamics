## source("~/R/SpherHam/plmbck.r")
##
##
# 
#  Tranformation of sperical hamonics plm onto Gaussian field 
#  for triangular truncation T(mmm) 
#
 plmbck <- function(plm,factor=1)
 {

   dv <- array(as.complex(0),c(nlong,nlatg))
   
   for ( j in 1:nlatg )
   {
     for ( m in 0:mmm )
     {
       lm <- iml[(m+1),(m+1):(mmm+1)]
       dv[(m+1),j] <- sum( plm[lm] * Ralp[lm,j] )
       if ( m > 0 ) dv[(nlong-m+1),j] <- Conj( sum( plm[lm] * Ralp[lm,j] ) )
     }
   }

   gauss <- mvfft(dv,inverse = TRUE)
   gauss <- Re( gauss * factor / nlong )
   return( gauss )
 }
