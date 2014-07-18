##' Function wmat
##' 
##' Function wmat
##' @param tempmat
##' @param imp
##' @param irank
##' @param jvar
##' @return dm 
wmat<-function(tempmat,imp,irank,jvar){
	dm								<-	matrix(0,irank,irank)
 	dm[1,1]							<-	tempmat[jvar,jvar]
 	dm[1,2:irank]					<-	tempmat[jvar,imp[1:(irank-1)]]
	dm[2:irank,1]					<-	tempmat[imp[1:(irank-1)],jvar]
	dm[2:irank,2:irank] <-  tempmat[imp[1:(irank-1)],imp[1:(irank-1)]]

	dm
  
	
}

