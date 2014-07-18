replot_checkdata <-
function(x,y,scale_text,maintext,OBS,ssX,s,coords=NULL){
	
	xtext <-  paste("Component 1 (",round(t(x)%*%x/ssX*100,1),"% of variation in X)",sep="")
	ytext <-  paste("Component 2 (",round(t(y)%*%y/ssX*100,1),"% of variation in X)",sep="")

	color <-  s %in% s[grep("[OUT]",s)]+1
        if (missing(coords)){
          plot(x,y,main=c(maintext,scale_text),xlab=xtext,ylab=ytext,col=color,pch=16)
          text(x,y,OBS,pos=4,col=color)
        }else{
          plot(x,y,main=c(maintext,scale_text),xlab=xtext,ylab=ytext,col=color,pch=16,xlim=coords$x,ylim=coords$y)
	  text(x,y,OBS,pos=4,col=color)
        }
}

