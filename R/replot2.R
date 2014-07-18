replot2 <-
function(edges,DATA,zoomwidth,maintext){
	
	if(sum(zoomwidth) == 0){
		zoomwidth[1]  <-  1
		zoomwidth[2]  <-  ncol(DATA)
	
	}else{
		zoomwidth[1]	<-	max(1,zoomwidth[1])
		zoomwidth[2]  <-  min(zoomwidth[2],ncol(DATA))
	}
 	
 	ymax	<-	max(DATA[,zoomwidth[1]:zoomwidth[2]])*1.05
 	par(new=FALSE)
 	for(i in 1:nrow(DATA)){
 		plot(zoomwidth[1]:zoomwidth[2],DATA[i,zoomwidth[1]:zoomwidth[2]],xlim=c(zoomwidth[1],zoomwidth[2]),ylim=c(0,ymax),type="l",col=i,main=maintext,xlab="",ylab="")
		par(new=TRUE)
	}
 	abline(v=edges,col="red",lwd=2)
}

