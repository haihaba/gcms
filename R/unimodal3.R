unimodal3 <-
function(C,obs,RT_LIMIT){
	
	cp	<-	excel_round(median(which.max(colMeans(t(C)))))
	if(!length(cp)){
		C		<-	C*0
		cp		<-	0

	}else{
		for(i in 1:obs){
			xpeak	<-	peak_pick(t(C[,i]))$xpeak
			mp		<-	which(as.logical(xpeak))
			if(length(mp)){
				if(which.min(abs(mp-cp)))
					mp	<-	mp[which.min(abs(mp-cp))]
				if(abs(mp-cp)<RT_LIMIT){
					D			<-	diff(C[1:mp,i])
					if(any(D<0))
						poss	<-	max(which(D<0))
					
					else
						poss	<-	1
					
					for(j in poss:1)
						C[j,i]	<-	min(C[j:mp,i])
					
					D			<-	diff(C[mp:nrow(C),i])
					
					if(length(which(D>0)))
						poss	<-	min(which(D>0))
					
					else
						poss	<-	FALSE

					if(poss)
						for(j in (mp-1+poss):nrow(C))
							C[j,i]	<-	min(C[mp:j,i])

					# Correction of strange peaks
					knew						<-  C[-which(diff(C[,i]) == 0),i]
					mpnew						<-	which.max(knew)
					k							<-	C[,i]*0
					k[1:length(knew)+mp-mpnew]	<-	knew
        			C[,i]	<-	k
        		
        		}else
        			C[,i]	<-	C[,i]*0
        	}else
        		C[,i]	<-	C[,i]*0
        }
	}
	out <-  list(C=C,cp=cp)
}

