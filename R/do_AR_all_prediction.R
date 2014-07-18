##' Function do_AR_all_prediction
##' 
##' Function do_AR_all_prediction
##' Funciton to do Alternate Regression on a prediction set
##' @param x
##' @param S
##' @param CP
##' @param RT_LIMIT
##' @return Cnew
do_AR_all_prediction <- function(x,S,CP,RT_LIMIT){
	
	var 		<-  nrow(x)
	mz  		<-  ncol(x)
  	C			<-	x%*%S%*%ginv(t(S)%*%S)
	C[C<0]	<-	0
	Cnew		<-	C*0
 	
 	for(i in 1:ncol(C)){
 		c		<-	as.matrix(C[,i],nrow=var)
    	Cnew[,i]	<-	unimodal_3(c,1,CP[i],RT_LIMIT)
 	}
	
	C	<-	Cnew
}

