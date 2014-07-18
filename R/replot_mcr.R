replot_mcr <-
function(x,y,maintext,incl,excl,numfiles,set1,set2,seltype=(!missing(set1) & !missing(set2)),plotlabels = TRUE,pred){
	
	if(seltype == 1){
		color <-  2*(1:numfiles %in% set1)+4*(1:numfiles %in% set2)+3*(1:numfiles %in% pred)+(1:numfiles %in% setdiff(1:numfiles,c(set1,set2,pred)))
	
	}else{
		color <-  2*(1:numfiles %in% incl)+(1:numfiles %in% excl)+3*(1:numfiles %in% pred)
	}
	
	plot(x,y,xlim=c(range(x)),ylim=range(y),main=maintext,xlab="t[1]",ylab="t[2]",col=color,cex=0.9,pch=16)

	if(plotlabels)
  		text(x,y,pos=4,cex=0.7)
}

