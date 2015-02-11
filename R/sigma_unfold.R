sigma_unfold  <-  function(predpath){
  
  
	if(!file.exists(file.path(predpath,"Aligned","shift.Rdata"))){
    cat("Error! You need to align data first!\n")
		return(NULL)
	
  }else{
    load(file.path(predpath,"Aligned","files.Rdata"))
		load(file.path(predpath,"Aligned","COLUMNID1.Rdata"))
		load(file.path(predpath,"Aligned","shift.Rdata"))
		OBSID1	<-	COLUMNID1
		load(file.path(predpath,"SETTINGS.Rdata"))
		DOBL	<-	SETTINGS$BC2
		
    if(file.exists(file.path(predpath,"Edges","edges.Rdata"))){
			load(file.path(predpath,"Edges","edges.Rdata"))
			load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
			ROWID1	<-	varnamegen(length(edges)-1,length(SCAN_RANGE),min(SCAN_RANGE))
			load(file.path(predpath,"Aligned","NUM_scans.Rdata"))
			
      for(i in 1:length(files)){
		    DATA_temp	<- numeric()
		    cat("Loading ",basename(files[i]),paste("(",i,"/",length(files),")",sep=""),"\n")
		    load(files[i])

		    if(shift[i]+tail(edges,1)>nrow(Xbc))    # Avoids "subscript out of bounds" if the shift is too big (outliers)
		        Xbc <-  rbind(Xbc,matrix(0,shift[i]+tail(edges,1)-nrow(Xbc),ncol(Xbc)))

		    for(k in 1:(length(edges)-1)){
		      X_temp <- Xbc[edges[k]:edges[k+1]+shift[i],]
		      
          if(DOBL){
            BL               <- apply(X_temp,2,function(X_temp) approx(c(1,length(X_temp)),X_temp[c(1,length(X_temp))],1:length(X_temp))$y)       # method = 'linear' by default
		        X_temp           <- X_temp-BL
		        X_temp[X_temp<0] <- 0
					}
          
		      DATA_temp <- cbind(DATA_temp, matrix(colSums(X_temp),nrow=1))
				}
        
		    if(i == 1)
		    	DATA  <-  matrix(0,length(files),length(DATA_temp))
				
        if(ncol(DATA) < ncol(DATA_temp))
					DATA  <-  cbind(DATA,matrix(0,nrow=nrow(DATA),ncol=ncol(DATA_temp)-ncol(DATA)))
			  
        DATA[i,1:ncol(DATA_temp)]	<-	DATA_temp
			}
      
			dir.create(file.path(predpath,"HMCR","REG"),showWarnings=FALSE)
			save(DATA,OBSID1,file=file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))
			save(edges,shift,ROWID1,files,file=file.path(predpath,"info.Rdata"))
			output <- DATA
		}else{
		  cat("Error! You need to set edges first!\n")
		  return(NULL)
		}
	}
}
