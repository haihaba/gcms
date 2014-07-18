##' Function do_AR_all
##' 
##' Function do_AR_all
##' Funciton to do Alternate Regression
##' @param x
##' @param projectpath
##' @param NL
##' @param RP
##' @param RT_LIMIT
##' @param SCAN_RANGE
##' @return C, S, INDEX, R2, noise, CP
do_AR_all <- function(x,projectpath,NL,RP,RT_LIMIT,SCAN_RANGE){
	
	var <-  dim(x)[1]
	mz  <-  dim(x)[2]
	obs <-  dim(x)[3]
	cat("\n[TIMEPOINTS,M/Z-Channels,Observations] = [",var,",",mz,",",obs,"]\n\n")
	R					<-	OK	<-	1
	R2					<-	0
	C					<-	S	<-	numeric()
	OK_out				<-	3
	X 					<-  t(matrix(aperm(x,c(2,1,3)),nrow=mz,ncol=var*obs))
	ind					<-	w_mor2(X)
	noise				<-	which(ind<NL)
	x[,noise,]			<-	0
	ssX					<-  sum(apply(x,3,ss))
	Xmean				<-	apply(x,2,function(x) apply(x,1,sum))
	rm(x)
	X[,noise]		<-	0
	max_R				<-	qr(X[,setdiff(1:length(ind),noise)])$rank
	S						<-	abs(pca(X,1)$p)
	vars				<-	which(as.logical(SCAN_RANGE))
	cat("-----------------------------\n")
	gc()
	
	while(OK_out>0){
		iter	<-	0
		C 		<-  numeric()
		if(max_R >= R){
			if(R > 1){
				p		<-	pure(Xmean[,vars],R,0.01)
				sp  <-  p$sp
				imp <-  p$imp
      			S	<-	matrix(0,mz,R)
        
       		for(i in 1:R)
       			S[vars[imp[i]],i]	<-	1
     		}
     		
     		R2			<-	0.001
     		I			<-	numeric()
     		dif			<-	1
     		C			<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)
     		
			
			
			C[C<0]		<-	0
			out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
			C	  		<-  out$C
			CP			<-  out$CP
			
			while(dif > 1e-6){
				r2			<-	R2
				iter		<-	iter+1
				tC  		<-  t(C)
				S			<-	t(ginv(tC%*%C,tol=.Machine$double.eps^2)%*%(tC%*%X))
				S[S<0]		<-	0
				S 			<-  apply(S,2,function(S) S/sum(S))
				
				if(any(is.na(S)))
					break
				
				C	<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)
					
				if(any(is.na(C)))
					break
				C[C<0]	<-	0
				out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
				C       <-  out$C
				CP      <-  out$CP
				R2		<-	1-ss(X-C%*%t(S))/ssX
				#				cat(iter,": ",R2,"| sum(C) = ",sum(C),"\n")
				
				if(iter == 50)
					dif <-  0
				
				else
					dif <-  abs(R2-r2)
			}
				
			if(any(is.na(S)) | any(is.na(C))){
					
				cat("Warning: NAs found.\n")
				C		<-	matrix(1,nrow(C),ncol(C))
				S		<-	matrix(0,nrow(S),ncol(S))
				R2		<-  0
			}
				
			INDEX	<-	matrix(0,obs,R)
				
			for(i in 1:R){
				c	<-	matrix(C[,i],var,obs)
				for(j in 1:obs){
					if(max(c[,j])>0)
						INDEX[j,i]  <-  min(which.max(c[,j]))
					else
						INDEX[j,i]  <-  -99
				}
			}
				
			if(sum(INDEX == -99) > 0){
				I <-  row(INDEX)[INDEX == -99]
				J <-  col(INDEX)[INDEX == -99]
	
				for(i in 1:length(J)){
					index				<-	INDEX[,J[i]]
					index				<-	index[-(index == -99)]
					INDEX[I[i],J[i]]	<-	median(index)
				}
			}
				
			I		<-  order(colMeans(INDEX))
			Y		<-  sort(colMeans(INDEX))
			C		<-	C[,I]
			S		<-	S[,I]
			CP		<-	CP[I]
			INDEX	<-	INDEX[,I]
			cat("Rank: ", R,"\n")
				
			if(R>1)
				PERMUTED_ORDER	<-	sum(diff(t(INDEX)) < -RP)
			else
				PERMUTED_ORDER	<-	0
					
			cat("PERMUTED_ORDER = ",PERMUTED_ORDER,"\n")
			cat("Number of iterations: ", iter,"\n")
			cat("R2X = ",round(R2,5),"\n")
				
		}else{
			R2				<-	0
			S				<- 	C 	<-  numeric()
			PERMUTED_ORDER	<-	99
		}
			
		s		<-	S
		c		<-	C
		r2		<-	R2
		
		if(r2>0){
			if(!PERMUTED_ORDER){
				C_out		<-	c
				S_out		<-	s
				CP_out		<-	CP	
				R2_out		<-	r2
				OK_out		<-	3
				OK2_out		<-	0
				cat("OK\n")
				
			}else{
				OK2_out	<-	1
				cat("NOT OK [1]\n")
			}
			
		}else{
			OK2_out	<-	1
			cat("NOT OK [2]\n")
		}
			
		OK_out	<-	OK_out-OK2_out
		R		<-	R+1
			
		cat("Timestamp: ",format(Sys.time(), "%X"),"\n")
		cat("-----------------------------\n")
	}
		
	if(!exists("C_out")){
		C_out		<- numeric()
		S_out		<-	numeric()
		R2_out	<-	numeric()
		CP_out	<-	numeric()
	}
		
	C		<-	as.matrix(C_out)
	S		<-	as.matrix(S_out)
	R2		<-	R2_out
	CP		<-	CP_out
	R		<-	ncol(C)
		
	cat("FINAL MODEL\n")
	cat("Rank:\t",R,"\n")
	cat("R2X:\t",R2,"\n")
	cat("-----------------------------")
		
			
	if(!is.null(R)){ #produces errors when R=1
	#if(R>1){
		INDEX	<-	matrix(0,obs,R)
			
		for(i in 1:R){
			c	<-	matrix(C[,i],var,obs)
			for(j in 1:obs)
				INDEX[j,i]	<-	min(which.max(c[,j]))
	
		}
				
		I		<-	order(round(colMeans(INDEX)))
		INDEX 	<-  sort(round(colMeans(INDEX)))
	}else
		INDEX	<-	numeric()

	gc()
	out <-  list(C=C,S=S,INDEX=INDEX,R2=R2,noise=noise,CP=CP)
}

