ARmethod1 <-
function(X,projectpath)
{
	require(MASS)
	X[X<0]	<-	0
	XX			<-	sweep(X,2,colMeans(X))
	vec    <-  pca(XX,min(c((dim(X)-1),15)))$vec
	SSX			<-	sum(XX^2)
	load(file.path(projectpath,"SETTINGS.Rdata"))
	limit		<-	SETTINGS$R2Xc
	R				<-	sum(cumsum(diag(t(vec)%*%vec)/SSX)<limit)+1
	# R f?r ej vara 0

	p		<-	pca(X,R)$p

	S				<-	abs(p)
	Sstart	<-	S
	diff		<-	1
	R2			<-	0
	z				<-	0
	zz			<-	100
	C				<-	X%*%S%*%ginv(t(S)%*%S)
	C[C<0]	<-	0
	Cstart	<-	C
	while(diff > 1e-6)
	{
  	z				<-	z+1
    S				<-	t(X)%*%C%*%ginv(t(C)%*%C)
    S[S<=0]	<-	0
		if(R == 1) # Applyfunktionen nedan fungerar endast f?r R>1
		{
		  S[,1]/sqrt(sum(S[,1]^2))
 		}else
			S[,1:R] 	<-  apply(S[,1:R],2,function(S) S/sqrt(sum(S^2)))

    CC	<-	t(S)%*%S
    if(sum(CC>0.95) > R)
    {
    	zz	<-	z
    	z		<-	100
		}
   	C				<-	X%*%S%*%ginv(t(S)%*%S)
    C[C<0]	<-	0
    R2new		<-	1-ss(X-C%*%t(S))/ss(X)
    diff		<-	abs(R2-R2new)
    R2			<-	R2new
    if(any(is.na(R2)))
        z	<-	100
    if(z == 100)
    {
        cat("Rank reduced after: ",zz, " iterations.\n")
        z			<-	0
        R			<-	R-1
        C			<-	Cstart[,1:R]
        diff	<-	1
        R2		<-	0
        zz  	<- 100
		}
	}
	cat("Rank:       ",R,"\n")
	cat("R2X:        ",R2,"\n")
	cat("Iterations: ",z,"\n")
	#do_log(projectpath,paste("Rank:       ",R,"\n"));
	#do_log(projectpath,paste("R2X:        ",R2,"\n"));
	#do_log(projectpath,paste("Iterations: ",z,"\n"));
	#do_log(projectpath,"----------------------")
	output	<-	list(C=C,S=S)
}

