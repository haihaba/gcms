##' Function do_DATA1
##' 
##' Function do_DATA1
##' This funcition is used to arrange the initial data table after import
##' @param INT
##' @param MZ
##' @param TIME
##' @param count
##' @param MZmin
##' @param MZmax
##' @param MZP
##' @return DATA
do_DATA1<-function(INT,MZ,TIME,count,MZmin,MZmax,MZP){
  MZ        <- round(MZ,digits = MZP)*10^MZP
  INT       <- round(INT)
  lenINT    <- length(INT)
  lencount  <- length(count)
  lenTIME   <- length(TIME)
  
  if(length(abs(diff(MZ))[abs(diff(MZ)) != 1 & abs(diff(MZ)) != (MZmax-MZmin)])){#, !all(count == (MZmax-MZmin+1))))
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
