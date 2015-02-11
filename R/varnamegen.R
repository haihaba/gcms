varnamegen <-
function(intervals,peks,start)
{
		#A <- cbind(paste("int ",sort(rep(1:intervals,peks)),"  M/Z ", rep(start:(start+peks-1),intervals)," ",sep=""))
		A <-  matrix(c(sort(rep(1:intervals,peks)),rep(start:(start+peks-1),intervals)),peks*intervals,2)
}

