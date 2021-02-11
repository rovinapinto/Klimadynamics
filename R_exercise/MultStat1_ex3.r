## source("~/R/KlimDyn/MultStat1_ex3.r")
##
## met731 - Klimadynamik und Statistik - WS2020/2021
## 
## Exercise III - 19.11.2020
##
##

 library(fields)
 library(ncdf4)
 library("maps")
 library("date")

 # inputfile
 #
 infile = "/home/nrisse/PEA/Lectures/13_ClimDyn/data/HadSLP2_slp-mnmean-real.nc"

## I.
##
 
 # Read netcdf file
 # 

 ncId <- nc_open(infile)
 print(ncId)

 lat  <- ncvar_get(ncId,varid="lat") ; nlat <- length(lat)
 lon  <- ncvar_get(ncId,varid="lon") ; nlon <- length(lon)
 time <- ncvar_get(ncId,varid="time") ; nt <- length(time) # units: days since 1800-1-1 00:00:00
 dtime <- as.Date(time, origin = "1800-01-01")

 Var <- ncvar_get(ncId,varid="slp") 
 VarName <- "Mean Sea Level Pressure" # units: mb or hPa

 # Flip dimension latitude to be ascending
 #
 lat <- lat[nlat:1]
 Var <- Var[,nlat:1,]

 # Select December months
 #
 index.dec <- date.mdy(dtime)$month==12
 VarDec <- Var[,,index.dec]
 dtimeDec <- dtime[index.dec]
 ntdec <- length(dtimeDec)
 rm(Var) ; rm(dtime)

 # Calculate zonal mean averages for december months
 #
  VarDecZ <- apply(VarDec,c(2,3),"mean")
  
 # Plot averaged zonal mean MSLP
 #
  plot(lat,apply(VarDec,c(2),"mean"),xlab="Latitude",ylab="MSLP zonal mean")
  lines(lat,apply(VarDec,c(2),"mean"),col="blue")

## 
##

 # Calculate matrix of anomalies (D)
 #
  D <- t(scale(t(VarDecZ),scale=F))
  
 # Weight variance with cos(pi/180*lat) according to area the values represent
 #
  W <-diag(sqrt(cos(lat*2*pi/360)))
  
 # Hovmoeller diagram (e.g., wikipedia)
 #
  image.plot(lat,dtimeDec,W%*%D,xlab="Latitude",ylab="MSLP anomalies - area weighted")
  
## 4. 
##

 # Estimate covariance matrix S of zonal mean MSLP and plot it
 #
  S <- 1/ntdec* W%*% D%*% t(D) %*% W
  
## 5. 
##

 # Derive and plot first eigenvector of S
 #
  PC <- eigen(S)
  E <- PC$vectors
  L <- PC$values
  A <- t(E) %*% W%*%D

 # Plot Eigenvalues and explained variance spectra
 #

 # Plot 1st and 2nd EOF
 #
  par(mfrow=(c(1,2)))
  plot(lat,E[,1],xlab="Latitude",main="1. EOF")
  lines(lat,E[,1],col="blue")
  plot(lat,E[,2],xlab="Latitude",main="2. EOF")
  lines(lat,E[,2],col="blue")

 # Plot principal components
 #
  par(mfrow=(c(2,1)))
  plot(dtimeDec,A[1,],type='l',
       ylab='PC', xlab='Time', main='PC of 1. EOF')
  plot(dtimeDec,A[2,],type='l',
       ylab='PC', xlab='Time', main='PC of 2. EOF')



 # Derive expansion coefficients for 1st and 2nd EOF, and plot them
 #
 
 # plot filtered D for qt variable
 #
  qt <- 2
  G <- E[,1:qt]
  Df <- G %*% t(G) %*%W %*%D

  par(mfrow=(c(1,2)))
  image.plot(lat,dtimeDec,Df,xlab="Latitude",ylab="Filtered anomalies")
  image.plot(lat,dtimeDec,Df-W%*%D,xlab="Latitude",ylab="Difference")
  
