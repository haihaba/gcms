##' Function sfd_par
##' 
##' Function sfd_par
##' @param X
##' @param num
##' @param exclude
##' @return object, j, md
sfd_par<-function(X,num,exclude = nrow(as.matrix(X))+1){
	
	N				<-  nrow(as.matrix(X))
	K				<-  ncol(as.matrix(X))
	
	if(num == (N-length(exclude))){
		sel <-  (1:N)[-exclude]
		return(list(object=sel))
  	}

	
	D			<-	obj_eu(X)
	diag(D)		<-	2*max(D)
 	MD			<-	j	<-	0
	
	for(run in 1:500){
		object	<-	sample((1:N)[-exclude])[1:num]
    	OK			<-	3
    	md			<-	0
    	number	<-	1
    
    	while(OK){
    		XD	<-	D[object,object]
    		OK	<-	OK-1
    		nk	<-	which.min(XD)
			k	<-	col(XD)[nk]
			
			if(OK)
				n	<-	row(XD)[nk]*(OK != 2) + k*(OK == 2)
			else
				n	<-	number
			
			number 		<-  number + (OK==0)
			
			if(!OK)
				OK	<-	(number <= num)
      		
      		obs			<-	object[-n]
      		ex_object	<-	setdiff((1:N)[-exclude],obs)
      		I			<-	which.max(apply(D[obs,ex_object],2,min))
      		
      		for(jj in I){
      			obs[num]	<-	ex_object[jj]
      			
      			if((md <= min(D[obs,obs])) & sum(abs(sort(obs)-sort(object))) > 0){
      				XD				<-	D[obs,obs]
      				diag(XD)		<-	0
          
					if(md == min(D[obs,obs])){
					
						if(j < sum(XD)/2){
							md			<-	min(D[obs,obs])
							object		<-	sort(obs)
							j			<-	sum(XD)/2
							number		<-	number*(OK != 1) + (OK == 1)
						
							if(OK>1)
								OK	<-	3
						}
					}else{
						md		<-	min(D[obs,obs])
						object	<-	sort(obs)	
						j		<-	sum(XD)/2
						number	<-	number*(OK != 1) + (OK == 1)
					
						if(OK>1)
							OK	<-	3
					}
				}
			}
    	}
    
    	XD				<-	D[object,object]
    	md				<-	min(D[object,object])
    
    	diag(XD)	<-	0
		j			<-	sum(XD)/2
    	object		<-	sort(object)
    
    	if(MD < md){
   		 	J	<-	j
   		 	sel	<-	object
    		MD	<-	md
    	}
	}
	list(object=sel,j=J,md=MD)
}

