##' A_E Function, Automatic Edges
##' 
##' Probably the funciton to set automatic edges for the processing windows
##' @param DATA
##' @param Nmin
##' @param Nmax
##' @return edges 
A_E <- function(DATA,Nmin,Nmax){
  edges	<-	c(1,length(DATA))
	N		<-	0.5
	while(max(diff(edges)) > Nmax){
		n		<-	which.max(diff(edges)) 
		range	<-	(edges[n] + Nmin):(edges[n+1] - Nmin)
		x		<-	E_peak_pick(DATA[range],median(DATA[range]),Nmin)
    	
	    if(!any(x$peak>0)){
	    	N	<-	N*1.5
	    	if(N > 7){
	    		n	<-	which.min(DATA[range])
	    		
	        	if(length(n) > 1){
	        		k	<-	round(length(n)/2)
	        		n	<-	n[k]
	        	}
				
				N			<-	0.5
				x$peak[n]	<-	1
			}
		}else
	    	N	<-	0.5
		edges	<-	sort(c(edges,range[which(x$peak>0)]))
	}
	
	edges <-  edges[!(edges == 1 | edges == length(DATA) | edges > max(which(DATA>0)) | edges < min(which(DATA>0)))]
}

