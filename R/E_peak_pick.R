##' Function E_peak_pick
##' 
##' Function E_peak_pick
##' @param x
##' @param NL
##' @param min_win
##' @return x 
E_peak_pick <- function(x,NL,min_win){
	require(signal)
	#xd1 <- -sav.gol(x,11,3,1)
	xd1 <-  sgolayfilt(x,3,11,1)
	N1	<-	which(x < NL)
	N2	<-  numeric(length(x))
	
	
	for(i in 3:(length(x)-2)){
		if(xd1[i-2] < 0  & xd1[i-1]<0 & xd1[i+1]>0 & xd1[i+2]>0 &  sum(N1 == i)==1)
    		N2[i]	<- 1
	}
	
	N	<-	intersect(N1,which(N2 == 1))
	
	if(length(N)==0)
		return()
	
	if(length(N) != 1)
		while(min(diff(N)) < min_win){
			p1	<- 	which.min(diff(N))
			p2	<-	p1+1
			
 			if(x[N[p1]] < x[N[p2]])
 				N <- N[-p1]
			else
				N <- N[-p2]

 			if(length(N) == 1)
 				break
		}
		
	xpeak		<-	numeric(length(x))
	xpeak[N]	<-	x[N]
	x			<-	list(peak = xpeak, out = x)

}

