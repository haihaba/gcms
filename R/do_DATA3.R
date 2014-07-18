do_DATA3<-function(INT,MZ,TIME,count,MZmin,MZmax,MZP){
  MZ        <- round(MZ,digits = MZP)*10^MZP
  INT       <- round(INT)
  lenINT    <- length(INT)
  lencount  <- length(count)
  lenTIME   <- length(TIME)

  DATA	    <- matrix(0,lenTIME,max(MZ))
  scans	    <- 0

  for(i in 1:lenTIME){
    if(count[i] != 0){
      mz<-MZ[(sum(count[0:(i-1)])+1):sum(count[1:i])]
      int<-INT[(sum(count[0:(i-1)])+1):sum(count[1:i])]

      DATA[i,mz] <- int 
  

    }
  }
  return(DATA)
}

