##' Function replot
##' 
##' Function replot
##' @param edges
##' @param DATA
##' @param zoomwidth
##' @param maintext
replot<-function(edges,DATA,zoomwidth,maintext){
	
  ## setting zoomwidth borders for no zoom
	if(sum(zoomwidth) == 0){
		zoomwidth[1]  <-  1
		zoomwidth[2]  <-  ncol(DATA)
	
  ## setting zoomwidth borders for zoom
	}else{
		zoomwidth[1]	<-	max(1,zoomwidth[1])
		zoomwidth[2]  <-  min(zoomwidth[2],ncol(DATA))
	}
  
   	
 	ymax	<-	max(DATA[,zoomwidth[1]:zoomwidth[2]])*1.05
 	#par(new=FALSE)
  
 	#for(i in 1:nrow(DATA)){
  averagedData<-apply(DATA,2,mean)
  plot(zoomwidth[1]:zoomwidth[2],averagedData[zoomwidth[1]:zoomwidth[2]], xlim=c(zoomwidth[1],zoomwidth[2]),ylim=c(0,ymax), type="l", main=maintext,xlab="",ylab="")
 	#	plot(zoomwidth[1]:zoomwidth[2],DATA[i,zoomwidth[1]:zoomwidth[2]],xlim=c(zoomwidth[1],zoomwidth[2]),ylim=c(0,ymax),type="l",col=i,main=maintext,xlab="",ylab="")
  #par(new=TRUE)
	#}
 	abline(v=edges,col="red",lwd=2)
}

