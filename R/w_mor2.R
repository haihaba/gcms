w_mor2 <-
function(X){
	
	X								<-	sweep(X,2,colMeans(X))
	ind								<-  apply(X,2,function(X) sqrt(sum(X^2))/sqrt(sum(diff(X)^2)))
	ind[!is.finite(ind)]  <-  0
	ind
}

