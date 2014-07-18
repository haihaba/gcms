peak_pick <-
function(x){
	
	xout				<-	x	<-	as.matrix(x)

	xd1         <-  sgolayfilt(xout,3,11,1)
 	NOISE_LIMIT	<-	median(x)
	N1					<-	which(x>NOISE_LIMIT)
	N2					<-	matrix(0,nrow(x),ncol(x))

	for (i in 3:(ncol(x)-2))
  		if(xd1[i-2] > 0)
			if(xd1[i-1] > 0)
				if(xd1[i+1] < 0)
					if(xd1[i+2] < 0)
						if(sum(N1 == i) == 1)
				   			N2[i]	<-	TRUE

	N	<-	which(as.logical(N2))
	
	if(length(N) > 1){
		
		while(min(diff(N)) < 3){
			p1	<-	min(which(diff(N) < 3))
			p2	<-	p1+1
			
			if(xout[N[p1]] < xout[N[p2]])
				N	<-	N[-p1]
			else
				N	<-	N[-p2]
			if(length(N) == 1)
				break
		}
 	}
 	xpeak			<-	matrix(0,nrow(x),ncol(x))
	xpeak[N]	<-	xout[N]
 	out <-  list(xpeak = xpeak, xout = xout)
}

