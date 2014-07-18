##' Function pca
##' 
##' Function pca
##' @param X
##' @param comp
##' @return vec, p
pca<-function(X,comp){
	
	#   Calculates Principal componets of X by calc eigvectors of X'*X or X*X'
	#   Depending on whats easiest to calculate....

	if(nrow(as.matrix(X)) < ncol(as.matrix(X))){
		tX  <-  t(X)
		p 	<-  eigen(X%*%tX)$vectors[,1:comp]
		
		if(comp>1){
			p <-  apply(p,2,function(p) tX%*%p)
			p <-  apply(p,2,function(p) p/sqrt(sum(p^2)))
		
		}else{
			p	<-	tX%*%p
			p 	<-	p/sqrt(sum(p^2))
		}

	}else
 		p		<-	eigen(t(X)%*%X)$vectors[,1:comp]

	vec		<-	X%*%p
	vecp	<-	list(vec=vec,p=p)
}

