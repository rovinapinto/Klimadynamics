## source("~/R/KlimDyn/KlimDyn1.r")
##
## Klimadynamik und Statistik 1
## January 2017
## Petra Friederichs
##
##
   graphics.off()
   rm(list=ls())
   library("fields")
   library("maps")

## 
 # Initialisation for spherical harminics
 #
   a <<-  6371000  # earth radius
   truncation <- 21         # choose triangular truncation T21 or T42
   cT <- sprintf("T%d",truncation)

   source("/data/KlimDyn1/SpherHam/Rinit.r")

## Read ISCCP data
 #
   lon.isccp <- seq(-179.5,179.5)
   lat.isccp <- seq(-89.5,89.5)
   nlat.isccp <- 180 ; nlon.isccp <- 360

   cVar <- 'TNR'
   fname <- sprintf("/data/KlimDyn1/ISCCP/%s_9195_ANN.txt", cVar)
   ascfile <- file(fname,"r")
   print(readLines(ascfile,n=7))
   tnr <- scan(ascfile)
   tnr <- matrix(tnr,nrow=nlon.isccp,ncol=nlat.isccp)
   close(ascfile)
  
## Plot original TOA Net Radiation on original grid
 #
 # X11()
   par(mfrow=c(2,1),mar=c(2,2,3,1),las=1,mgp=c(3, 0.6, 0))
   image.plot(lon.isccp,lat.isccp,tnr,xlab="",ylab="",main="TOA Net Radiation (W/m2)",tcl=-0.3,zlim=c(-130,130))
   map("world", col="black",xaxt="s", yaxt="s",add=T)

## Interpolation on Gaussian grid
 #
   lon.isccp2 <- c(lon.isccp,lon.isccp+360)
   tnr2       <- rbind(tnr,tnr)
   obj       <- list(x=lon.isccp2, y=lat.isccp, z=data.matrix(tnr2))
   loc       <- make.surface.grid(list(x=long, y=latg))
   tnr.gauss <- array(interp.surface(obj, loc),c(nlong,nlatg))

 # Plot original TOA Net Radiation on Gaussian grid
 #
   par(mfrow=c(2,1),mar=c(2,2,3,1),las=1,mgp=c(3, 0.6, 0))
   image.plot(long,latg,tnr.gauss,xlab="",ylab="",main=sprintf("TOA Net Radiation (%s)",cT),zlim=c(-130,130))
   map("world2", col="black",xaxt="s", yaxt="s",add=T)

## Transform to sperical harmonics
 #
   tnr.plm <- plmhin(tnr.gauss)

 # Transform back to Gaussian grid
 #
   tnrTest.gauss <- plmbck(tnr.plm)

 # Plot backtransformed field on Gaussian grid and differences to original field
 #
# X11()
   par(mfrow=c(2,1),mar=c(2,2,3,1),las=1,mgp=c(3, 0.6, 0))
   image.plot(long,latg,tnrTest.gauss,xlab="",ylab="",main=sprintf("TOA Net Radiation (%s)",cT),zlim=c(-130,130))
   map("world2", col="black",xaxt="s", yaxt="s",add=T)
   image.plot(long,latg,tnrTest.gauss-tnr.gauss,xlab="",ylab="",main=sprintf("TOA Net Radiation (%s)",cT))
   map("world2", col="black",xaxt="s", yaxt="s",add=T)

## 1. Calculate and plot zonal mean by setting coefficients of all non-zonal harmonics to zero and transform back
 #
## 2. Calculate and plot zonal eddies by setting coefficients of all zonal harmonics to zero and transform back
 #
## 3. Calculate and plot transport potential
 #
