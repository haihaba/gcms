##' Function unimodal2
##' 
##' Function unimodal2
##' @param C
##' @param R
##' @param var
##' @param obs
##' @param RT_LIMIT
##' @return C, CP
unimodal2<-function(C,R,var,obs,RT_LIMIT){
	
	N 	<-  nrow(C)
	CP  <-  numeric(R)
	for(i in 1:R){
		c       <-  matrix(C[,i],var,obs)
		out		<-	unimodal3(c,obs,RT_LIMIT)
		C[,i]	<-	matrix(out$C,N,1)
		CP[i]	<-	out$cp
	}
	out <-  list(C=C,CP=CP)
}

