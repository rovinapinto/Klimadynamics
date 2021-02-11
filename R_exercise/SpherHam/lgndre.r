## Computes associated Legendre polynomials.
## Adapted from LGNDR function (IMSL) 
## Only triangular trunction
## Petra Friederichs 2012
#
 lgndre <- function(mmm,mmjj,latg,slatg,clatg)
 {
#                                                                                 
#     Calculates legendre polynomials (ALP) and their derivatives                 
#     (DALP) at the JL'th latitude (JL is in LEGAU) using                         
#     recurrence relationships.                                                  
#                                                                                
      alp  <- array(0,c(mmjj,nlatg)       )                                       
      dalp <- array(0,c(mmjj,nlatg)       )                                       
                                                                                 
## Set P(0,0) and P(0,1)                                                      
#                                                                                
      alp[1,] <- sqrt(.5)                                                         
      F1M=sqrt(1.5)                                                              
      alp[2,] <- F1M * slatg                                                        
      dalp[1,] <- 0

## Loop over wavenumbers                                                      
#                                                                                
      LM=2                                                                       
      mmo <- 0

      for ( m1 in 1:mmm ) {                                                              
         m <- m1-1                                                               
         M2M <- m+m                                                                 

         E2 <- sqrt( m+m + 3 )                                                         

         if ( m > 0 & m == mmo ){
            F2M <- -F1M * clatg/sqrt(m+m)                                               
            F1M <-  F2M * E2                                                           
         } else if ( m > 0 ) {
            F2M <- -F1M * clatg/sqrt(m+m)                                               
            F1M <-  F2M * E2                                                           
            LM <- LM+1                                                              
            alp[LM,] <- F2M                                                       
            LM <- LM+1                                                              
            alp[LM,] <- F1M * slatg                                                 
            dalp[(LM-1),] <- -m * alp[LM,]/E2                                      
         }

         m2 <- m+2                                                                  
         MMO <- m+1                                                              
         JFM <- mmm

         if ( mmm >= m2 ) {                                                     
            K = LM - m2 + 1                                                            
                                                                                
##          Loop over degree N                                                   
#                                                                                
            for ( l in m2:mmm )
            {                                                   
               AN=l                                                              
               AN2=l*l                                                           
               ANM2=(l-1)*(l-1)                                                  

               E1 <- sqrt( ((l-1)*(l-1) - m*m)/(4*(l-1)*(l-1) - 1) )                                 
               E2=sqrt((4*l*l - 1)/(l*l - m*m))                                   
               alp[K+l,] <- E2 * ( slatg*alp[K+l-1,] - E1*alp[K+l-2,] )               
               dalp[K+l-1,] <- (1-l)*alp[K+l,]/E2 + l*E1*alp[K+l-2,]        
            }                                                            
            LM=LM+mmm-m2+1                                                      
         }                                                  
         dalp[LM,] <- -l* slatg*alp[LM,]+(l+l+1)*alp[LM-1,]/E2             
      }
  return( list( alp = alp , dalp = dalp ) )
 }
