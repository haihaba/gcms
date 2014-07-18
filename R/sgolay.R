##' Function sgolay
##' 
##' Function sgolay
##' @param p
##' @param n
##' @param ts
##' @return F 
sgolay<-function(p, n, m = 0, ts = 1){ 
	
	library(MASS)
	if (n %% 2 != 1)
		stop("sgolay needs an odd filter length n")
	
	if (p >= n)
		stop("sgolay needs filter length n larger than polynomial order p")

	F = matrix(0., n, n)
	k = floor(n/2)
	for (row  in  1:(k+1)) {
    	C = ( ((1:n)-row) %*% matrix(1, 1, p+1) ) ^ ( matrix(1, n) %*% (0:p) )
    	A = ginv(C)
    	F[row,] = A[1+m,]
	}
  
	F[(k+2):n,] = (-1)^m * F[k:1,n:1]
  
  	if (m > 0)
  		F = F * prod(1:m) / (ts^m)
	
	class(F) = "sgolayFilter"
  	F
}

