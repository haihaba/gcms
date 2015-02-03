##' Function Baseline
##' 
##' Function Baseline
##' @param X              raw imported Data. Time as rows and mz as columns
##' @param projectpath    baseproject directory
##' @param FL             Filterlength
##' @return Xbc           moving average filtered data table, one file
baseline <- function(X,projectpath,FL = 0){
	
  N <-  nrow(X)
	K <-  ncol(X)
	#X	<-	sweep(X,2,apply(X,2,min))
	
  if(FL){
		Xbc								<-	X*0
		Xbc[1:FL,]						<-	X[1:FL,]
		Xbc[(nrow(X)-FL+1):nrow(X),] 	<-	X[(nrow(X)-FL+1):nrow(X),]

		
		start<-FL+1
		end<-dim(X)[1]-FL
    
    ## moving average filter along rows (scantimes)
		Xbc[start:end,]<-filter(X,filter=rep(1/(FL*2+1),(FL*2+1)))[start:end,]
		
				
		if(any(Xbc<0))
			Xbc[Xbc<0]	<-	0
	}else
		Xbc <-  X
	return(Xbc)
}

