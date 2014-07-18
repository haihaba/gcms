##' Function sgolayfilt
##' 
##' Function sgolayfilt
##' @param x
##' @param p
##' @param n
##' @param m
##' @param ts
sgolayfilt<-function(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
	
	len = length(x)
	if (class(p) == "sgolayFilter" || (!is.null(dim(p)) && dim(p) > 1)){
		F = p
		n = nrow(F)
  	}else
		F = sgolay(p, n, m, ts)
		k = floor(n/2)
		
		#z = filter(F[k+1,n:1], 1, x)
		
		filt  <-  F[k+1,n:1]
		z 		<-	na.omit(stats:::filter(c(rep(0,length(filt) - 1), x), filt, sides = 1))
  c(F[1:k,] %*% x[1:n], z[n:len], F[(k+2):n,] %*% x[(len-n+1):len])
}

