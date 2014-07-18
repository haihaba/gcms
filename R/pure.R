pure <-
function(d,nr,f){
	
	# [sp,imp]=pure(d,nr,f)
	# sp purest row/column profiles
	# imp indexes of purest variables
	# d data matrix; nr (rank) number of pure components to search
	# if d(nspectra,nwave) imp gives purest nwave => sp are conc. profiles (nr,nspectra)
	# if d(nwave,nspectra) imp gives purest nspectra => sp are spectra profiles (nr,nwave)
	# f percent of noise allowed respect maximum of the average spectrum given in % (i.e. 1% or 0.1%))

	#[nrow,ncol]=size(d);
	# calculation of the purity spectrum
 	w   		<-  p	<-	s	<-	matrix(0,nrow=max(nr,1),ncol=ncol(d))
	f			<-	f/100
	s[1,]		<-	apply(d,2,sd)
	m			<-	colMeans(d)
	ll			<-	s[1,]^2+m^2
	f			<-	max(m)*f
	p[1,]		<-	s[1,]/(m+f)

	imp 		<-  numeric()
	mp 			<-  max(p[1,])
	imp[1]  <-  which.max(p[1,])
	# calculation of the correlation matrix

	l	<-	sqrt(s[1,]^2+(m+f)^2)
	for(i in 1:nrow(d))
		d[i,]	<-	d[i,]/l
		
	c	<-	(t(d)%*%d)/nrow(d)

	# calculation of the weights
	# first weight
	w[1,]	<-	ll/l^2
	p[1,]	<-	w[1,]*p[1,]
	s[1,]	<-	w[1,]*s[1,]

	# next weights
	if(nr > 1){
		for(i in 2:nr){
			for(j in 1:ncol(d)){
				dm			<-	wmat(c,imp,i,j)
				w[i,j]	<-	det(dm)
				p[i,j]	<-	p[1,j]*w[i,j]
				s[i,j]	<-	s[1,j]*w[i,j]
			}
			
			mp[i] 	<-  max(p[i,])
			imp[i]  <-  which.max(p[i,])
		}
	}
	
	sp	<-	normv2(t(d[,imp[1:nr]]))
	purelist  <-  list(sp=sp,imp=imp)
}

