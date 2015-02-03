##' Function do_DATA1
##' 
##' Function do_DATA1 is the method used with Agilent Chemstation data 
##' and is set as default. It is called from read_cdf.R
##' This function is used to arrange the initial data table after import
##' @param INT    intensity values from CDF import
##' @param MZ     mass values from CDF import
##' @param TIME   scantime values from CDF import
##' @param count  point count from CDF import
##' @param MZmin  min value mass range from CDF import
##' @param MZmax  max value mass range from CDF import
##' @param MZP    mass precision from settings file
##' @return DATA  rawdata matrix with rows for timepoints
##' and columns for mz values. The mz values range from 
##' zero to maxMZ.
do_DATA1<-function(INT,MZ,TIME,count,MZmin,MZmax,MZP){
  
  ## prepare/process rawdata
  MZ        <- round(MZ,digits = MZP)*10^MZP
  INT       <- round(INT)
  lenINT    <- length(INT)
  lencount  <- length(count)
  lenTIME   <- length(TIME)
  
  ## check if all mass values from min to max are abundant. This is usually not 
  ## the case in Agilent datasets. Therefore the dataset has to be constructed
  ## by looping through the individual datavalues. Otherwise, a matrix assignment
  ## can be done directly
  if(length(abs(diff(MZ))[abs(diff(MZ)) != 1 & abs(diff(MZ)) != (MZmax-MZmin)])){ #, !all(count == (MZmax-MZmin+1))))
    cat("Using alternative method to arrange data matrix (slower)\n")

    DATA  <-  matrix(0,lenTIME,max(MZ))
    scans <-  min(which(count != 0))

    for(i in 1:lenINT){
      DATA[scans,MZ[i]]  <- DATA[scans,MZ[i]] + INT[i]
      if(i != lenINT)
        if(MZ[i+1] < MZ[i])
          scans  <- scans + min(which(count[(scans+1):lencount] != 0))
      }
    }else
      DATA <- cbind(matrix(0,nrow=lenTIME,ncol=max((min(MZ)-1),0)),matrix(INT,nrow=lenTIME,ncol=(max(MZ)-min(MZ)+1),byrow=TRUE))
  return(DATA)
}
