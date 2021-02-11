## source("~/R/SpherHam/plmhin.r")
##
##
# 
#  Tranformation of Gaussian field in sperical hamonics plm
#  for triangular truncation T(mmm) 
#
 plmhin <- function(gauss)
 {

   plm <- NULL
   dv <- mvfft(gauss)
   
   for ( m in 0:mmm )
   {
     for ( l in m:mmm )
     {
       lm =  iml[(m+1),(l+1)]
       plm[lm] <- 0.5 * sum( gauwch[] * Ralp[lm,] * dv[(m+1),] )
     }
   }

  return( plm )
 }
