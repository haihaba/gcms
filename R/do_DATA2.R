do_DATA2<-function(INT,MZ,TIME,count,MZmin,MZmax,MZP){
  MZ        <- round(MZ,digits = MZP)*10^MZP
  INT       <- round(INT)
  lenINT    <- length(INT)
  lencount  <- length(count)
  lenTIME   <- length(TIME)

  DATA	    <- matrix(0,lenTIME,max(MZ))
  scans	    <- 0

  for(i in 1:lenTIME){
    if(count[i] != 0){
      for(k in 1:count[i]){
        mz<-MZ[sum(count[1:i])-count[i]+k]
        DATA[i,mz]	<-	DATA[i,mz] + INT[sum(count[1:i])-count[i]+k]
        }
      }
    }
  return(DATA)
}

