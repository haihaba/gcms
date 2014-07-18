##' Function replot_align
##' 
##' Function replot_align
##' @param DATA
##' @param ADATA
##' @param start
##' @param stop
##' @param zoomwidth
##' @param datasource
##' @param maintext
##' @param targetfil
replot_align<-function(DATA,ADATA,start,stop,zoomwidth = c(0,0),datasource,maintext,targetfile){
	
	zoomwidth[1]	<-	max(0,zoomwidth[1])
	zoomwidth[2]	<-	min(zoomwidth[2],ncol(DATA))
	if(sum(zoomwidth) == 0){
		screen(1)
		for(i in 1:nrow(DATA)){
	    	plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim=c(0,ncol(DATA)*1.01),ylim=c(0,max(DATA[,min(start):max(stop)])*1.01),main=c(maintext,paste("Target file: ", targetfile)),xlab="",ylab=datasource)
			par(new=TRUE)
		}
		
		if(!missing(ADATA)){
			
			screen(2)
			for(i in 1:nrow(ADATA)){
	    		plot(1:length(ADATA[i,]),ADATA[i,],type="l",col=i,xlim=c(0,ncol(ADATA)*1.01),ylim=c(0,max(DATA[,min(start):max(stop)])*1.01),main="Aligned data",xlab="Time (scans)",ylab=datasource)
				par(new=TRUE)
			}
		}
	}else{
		
		## Plot unaligned data in upper half of graph:
		xlimit  <-  c(max(0,zoomwidth[1]),min(zoomwidth[2],ncol(DATA)))
	  	#ylimit  <-  c(0,max(ADATA[,zoomwidth[1]:zoomwidth[2]],DATA[,min(start):max(stop)])*1.01)  ### did not adjust the ylimit
	  	ylimit  <-  c(0,max(ADATA[,zoomwidth[1]:zoomwidth[2]])*1.01)
	  	screen(1)
	  	
	  	for(i in 1:nrow(DATA)){
	  		plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim = xlimit,ylim = ylimit,main=c(maintext,paste("Target file: ", targetfile)),xlab="",ylab=datasource)
			par(new=TRUE)
	 	}
	 	
    	if(!missing(ADATA)){
			## Plot aligned data in lower half of graph:
			screen(2)
    	for(i in 1:nrow(ADATA)){
    		plot(1:length(ADATA[i,]),ADATA[i,],type="l",col=i,xlim = xlimit,ylim = ylimit,main="Aligned data",xlab="Time (scans)",ylab=datasource)
			par(new=TRUE)
			}
		}
	}
}

