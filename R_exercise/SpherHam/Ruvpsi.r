##    source("~/R/SpherHam/uvpsi.r")
##
##
#
 uvpsi <- function(plm)
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
  vplm <- 0+1i * mplm * plm

 # U component
 #
  uplm <- ufak1*plm2 + ufak2*plm1 

  return( list( uplm=uplm , vplm=vplm ) )
 }

#     EPS(L,M) = SQRT((FLOAT(L*L)-FLOAT(M*M))/(4.*FLOAT(L*L)-1.))   
#
#     DO M=0,MMM
#      DO L=M,MMM
#      LM = IML(M,L)
#      UFAK1(LM) = FLOAT(L-1) * EPS(L,M)
#      UFAK2(LM) = - FLOAT(L+2) * EPS(L+1,M)
#      ENDDO
#      ENDDO


#BROUTINE UVPSI(PSHAKO,UPHAKO,VPHAKO)
#      include 'para.h'
#      include 'comkfk.h'
#      COMPLEX PSHAKO(MMJJ),UPHAKO(MMJJ),VPHAKO(MMJJ)                    00001780
#      COMPLEX PSH1,PSH2                                                 00001790
#C.                                                                      00001810
#C.    UVPSI BERECHNET AUS DEN SPHERICAL HARMONICS DER STROMFKT          00001820
#C.    DIE SPHERICAL HARMONICS DES U-WINDES UND V-WINDES                 00001830
#C.                                                                      00001840
#C.    V - KOMPONENTE                                                    00001850
#C.                                                                      00001860
#      DO 5 M=0,MMM                                                      00001880
#      DO 10 L=M,MMM
#     LM = IML(M,L)
#      VPHAKO(LM) = CMPLX(0.,1.)*FLOAT(M)*PSHAKO(LM)                     00001920
#  10  CONTINUE                                                          00001930
#  5   CONTINUE                                                          00001940
#C.                                                                      00001950
#C.    U - KOMPONENTE                                                    00001960
#C.                                                                      00001970
#      DO 15 M=0,MMM                                                     00001990
#      DO 20 L=M,MMM
#
#      LM = IML(M,L)
#
#      LMM1 = MAX0(LM-1,1)                                               00002030
#      PSH2 = PSHAKO(LMM1)
#      IF(L.EQ.M) PSH2 = CMPLX(0.,0.)

#      LMP1 = MIN0(LM+1,MMJJ)                                            00002040
#      PSH1 = PSHAKO(LMP1)                                               00002050
#      IF(L.GE.MMM) PSH1 = CMPLX(0.,0.)
#
#      UPHAKO(LM)=UFAK1(LM)*PSH2+UFAK2(LM)*PSH1
#
#  20  CONTINUE                                                          00002090
#  15  CONTINUE                                                          00002100
#      RETURN                                                            00002110
#      END                                                               00002120
