##' Function scaling
##' 
##' Function scaling
##' @param Xin
##' @param K
##' @return Xout, scale_text
scaling <-function(Xin,K){
	
	if(missing(K))
		K	<-	menu(c("Centering","UV-Scaling","Pareto","None"),title="Select scaling function for X-matrix")
	if(K == 1){
		#Xout		<-	cent(Xin)   #sweep(Xin,2,colMeans(Xin),check.margin=FALSE)
    Xout <- sweep(Xin,2,colMeans(Xin),check.margin=FALSE)
		scale_text	<-	"DATA is mean centered"
	}else if( K==2){
		#Xout <-  uv_sc(Xin)   #sweep(Xin,2,colMeans(Xtr),check.margin=FALSE)/(apply(Xin,2,sd) + (apply(Xin,2,sd) == 0))  ## Last term is an adjustment to prevent division by zero
    Xout <- sweep(Xin,2,colMeans(Xin),check.margin=FALSE)/(apply(Xin,2,sd) + (apply(Xin,2,sd) == 0))
    	scale_text	<-	"DATA is centered and UV-scaled"
	
	}else if( K==3){
		#Xout <-  par_sc(Xin)  #sweep(Xin,2,colMeans(Xin),check.margin=FALSE)/sqrt((apply(Xin,2,sd) + (apply(Xin,2,sd) == 0) ))  ## Last term is an adjustment to prevent division by zero
    Xout <-  sweep(Xin,2,colMeans(Xin),check.margin=FALSE)/sqrt((apply(Xin,2,sd) + (apply(Xin,2,sd) == 0)))
		scale_text	<-	"DATA is centered and Pareto-scaled"
  	}else if(K == 4 || K == 0){
  		Xout	<-	Xin
		scale_text	<-	"DATA neither centered nor scaled"
	}
	
	scale <-  list(Xout=Xout,scale_text=scale_text)
}

