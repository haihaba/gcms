##' Function replot
##' 
##' Function replot
##' @param plotlimits
##' @param plotedges
##' @param plotfigure
##' @param ymax
##' @param edges
##' @param DATA
##' @param num
##' @param zoomwidth
##' @param maintext
##' @param min_win
##' @param max_win
replot<-function(plotlimits=FALSE,plotedges=FALSE,plotfigure=FALSE,ymax,edges,DATA,num,zoomwidth,maintext,min_win,max_win){
	
	if(num == 0 & !missing(DATA))
		ymax  <-	max(DATA)*1.05
		
	if(plotfigure){
		zoomwidth	<-	ifelse(num == 0,length(edges)-1,zoomwidth)
		num			<-	max(num,1)
		L			<-	edges[num+zoomwidth]-edges[num]
		xmin		<-	edges[num]-L/100
		xmax		<-	edges[num+zoomwidth]+L/100
		ymax		<-	max(DATA[,edges[num]:edges[num+zoomwidth]])*1.05
		par(new=FALSE)
		
		for(i in 1:nrow(DATA)){
			
			plot(edges[num]:edges[num+zoomwidth],DATA[i,edges[num]:edges[num+zoomwidth]],xlim=c(xmin,xmax),ylim=c(0,ymax),type="l",col=i,main=c(maintext,paste("Window: ",paste("[",num,"-",num+zoomwidth-1,"] (",length(edges)-1,")",sep=""))),xlab="",ylab="")
			par(new=TRUE)
		}
	}

	if(plotlimits){
		m	<-	edges_move_limit(edges,min_win,max_win)
		segments(x0=edges,y0=rep(ymax,length(edges)),x1 = m$right, y1=rep(ymax,length(edges)),lwd=2.5,col="blue")
		segments(x1=m$left,y1=rep(0,length(edges)),x0=edges, y0=0,lwd=2.5,col="green")
	}

	if(plotedges)
		abline(v=edges,col="red",lwd=2)
}

