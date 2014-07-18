baseline <-
function(X,projectpath,FL = 0){
	#dyn.load("filterlength.so")
	N <-  nrow(X)
	K <-  ncol(X)
	#X	<-	sweep(X,2,apply(X,2,min))
	if(FL){
		Xbc								<-	X*0
		Xbc[1:FL,]						<-	X[1:FL,]
		Xbc[(nrow(X)-FL+1):nrow(X),] 	<-	X[(nrow(X)-FL+1):nrow(X),]

		#out	<-	.C("filterlength",as.integer(FL),as.integer(nrow(X)),as.integer(ncol(X)),as.double(X),results = as.double(Xbc))
		#Xbc	<-  matrix(out[[5]],N,K)
		
		start<-FL+1
		end<-dim(X)[1]-FL
		Xbc[start:end,]<-filter(X,filter=rep(1/(FL*2+1),(FL*2+1)))[start:end,]
		
				
		if(any(Xbc<0))
			Xbc[Xbc<0]	<-	0
	}else
		Xbc <-  X
	return(Xbc)
}

