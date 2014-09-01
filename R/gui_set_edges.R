##' Function gui_set_edges
##' 
##' Function gui_set_edges
##' @export
##' @param projectpath
gui_set_edges<-function(projectpath){
	
	X11.options(type="Xlib")
	dir.create(file.path(projectpath,"Edges"),showWarnings=FALSE)
	#do_log(projectpath,"Setting edges...")
		
	load(file.path(projectpath,"Aligned","files.Rdata"))
	load(file.path(projectpath,"Aligned","shift.Rdata"))
	load(file.path(projectpath,"maxMZ.Rdata"))
	min_win 			<-	120
	max_win 			<-	450
	edges				<-	numeric()
	deletededges		<-	numeric()
	menuvector  		<-	c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)
	menulist 			<-	c("Load data","Load old edges","Auto edges","Plot controls","Save windows as .png","Quit")
	k  					<-	1
	X11()
	
	while(k){
		cat("Number of edges:\t",length(edges),"\n============================\n")
		k  <-  menu(menulist[menuvector],title="Set edges")
		
    	if(k == sum(menuvector) | !k){
    		k <- 0
    	}else if(menulist[menuvector][k] == "Load data"){
    		datasource  <-  menu(choices=c("Total ion Current","Basepeak Chromatogram"),title = "Choose type of data to load")
			if(datasource){
				if(datasource == 1){
					load(file.path(projectpath,"Aligned","TIC.Rdata"))
		    		DATA  		<-  TIC
		    		maintext  <-  "Total ion current"
		    		rm(TIC)
		    		
		    	}else if(datasource == 2){
		    		load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
		    		DATA		<-	BASEPEAK
		    		maintext  	<-	"Basepeak Chromatogram"
					rm(BASEPEAK)
				}
				
				#graphics.off()
				
				for(i in 1:nrow(DATA)){
					plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim=c(0,ncol(DATA)),ylim=c(0,max(DATA)*1.01),main=maintext,xlab="",ylab="")
		  			par(new=TRUE)
 				}
 				
 				edges <-  c(0,ncol(DATA))
				menuvector[c(2,3,4)]	<-	TRUE
				
			}
			
		}else if(menulist[menuvector][k] == "Load old edges"){
			edgesdata	<-	tk_choose.files(file.path(projectpath,"Edges","edges.Rdata"),caption="Select edges data file",multi = FALSE)
			
			if(!length(grep(pattern="^.+\\.([rR][dD][aA][tT][aA])$",edgesdata))){
				cat("File is not .Rdata!\n\n")
			
			}else{
				load(edgesdata)
				
				if(!length(edges) | !is.numeric(edges)){
					cat("Invalid edges data!\n")
				
				}else{
					abline(v=edges,col="red",lwd=2)
					menuvector[5]	<-	TRUE
				}
			}
		
		}else if(menulist[menuvector][k] == "Auto edges"){
			OK  <-  menu(choices=c("Yes","No"),title="Current edges will be deleted! Continue anyway?")
			if(OK == 1){
				
				while(!is.numeric(min_win	<-	as.numeric(readline(prompt = "Input smallest possible window (scans) (120 default): "))))
        			cat("Invalid input!\n")
		
				if(is.na(min_win))
					min_win <-  120
		
				while(!is.numeric(max_win	<-	as.numeric(readline(prompt = "Input largest possible window (scans) (450 default): "))))
					cat("Invalid input!\n")
				if(is.na(max_win))
					max_win <-  450
				if(length(edges))
					abline(v=edges,col="white",lwd=2)
				edges	<-	A_E(colSums(DATA),min_win,max_win)
				abline(v=edges,col="red",lwd=2)
				menuvector[5]	<-	TRUE
			}
			
		}else if(menulist[menuvector][k] == "Save windows as .png"){
			cat("Storing windows as .png in ",file.path(projectpath,"Edges",maintext,"png"),"\n")
 			dir.create(file.path(projectpath,"Edges",maintext,"png"),showWarnings = FALSE,recursive = TRUE)
			X11.options(type="cairo")
			x11()
			for(i in 1:(length(edges)-1)){
					replot(plotlimits=TRUE,plotedges=TRUE,plotfigure=TRUE,edges=edges,DATA=DATA,num=i,zoomwidth=1,maintext=maintext,min_win=min_win,max_win=max_win)
				savePlot(file.path(projectpath,"Edges",maintext,"png",paste("win",ifelse(i<=99,ifelse(i<=9,paste("00",i,sep=""),paste("0",i,sep="")),i),sep="")),type="png")
			}
   			dev.off()
   			X11.options(type="Xlib")
   			cat("Done!\n\n")
	 	
	 	}else if(menulist[menuvector][k] == "Plot controls"){
	 		zoommenu  	<-  c("Zoom in","Zoom default","Zoom out","Pan left","Pan right","Move edges","Set edges","Delete edges","Back")
			zoom 		<-	1
			zoomwidth 	<-  c(0,0)
		  	
		  	while(zoom){
		  		
		  		if(sum(zoomwidth) != 0)
		  			cat("\n\nCurrent zoom width:\t",zoomwidth[2] - zoomwidth[1]," scans.\n")
		  		cat("============================\n")
 				zoom  <-  menu(choices=zoommenu,title = "Plot Controls v.2")
				
				if(zoom == 9 | !zoom){
					zoom  		<-  0
					zoomwidth <-  c(0,0)
  		  	# Zoom back to default after quit
   					replot2(edges,DATA,zoomwidth,maintext)
				
				}else if(zoom == 1){ #Zoom in
					cat("\n\nSelect two points on the x-axis to zoom in.\n\n")
					zoomwidth			<-	round(sort(locator(2,"p")$x))
					replot2(edges,DATA,zoomwidth,maintext)
				
				}else if(zoom == 2){ #Zoom default
					zoomwidth	<-	c(0,0)
   					replot2(edges,DATA,zoomwidth,maintext)
   					
				}else if(zoom == 3){ #Zoom out
					
					if(sum(zoomwidth) == 0){
				    	cat("You must zoom in before you zoom out!\n")
					
					}else{
						tmp1		<-	zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
						tmp2		<-	zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
					  	zoomwidth	<-	c(tmp1,tmp2)
					  	
						if(diff(zoomwidth) > ncol(DATA))
							zoomwidth <-  c(0,0)
   						replot2(edges,DATA,zoomwidth,maintext)
					}
				
				}else if(zoom == 4){ #pan left
					if(sum(zoomwidth) == 0){
						cat("You must zoom in before you can pan!\n")
					
					}else if(zoomwidth[1] <= 0){
						cat("Beginning of chromatogram.\n")
					
					}else{
						tmp1  		<-  zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
						tmp2  		<-  zoomwidth[2] - (zoomwidth[2]-zoomwidth[1])/3
						zoomwidth <-  c(tmp1,tmp2)
						replot2(edges,DATA,zoomwidth,maintext)
					}
				
				}else if(zoom == 5){ #pan right
					
					if(sum(zoomwidth) == 0){
						cat("You must zoom in before you can pan!\n")
					
					}else if(zoomwidth[2] >= ncol(DATA)){
						cat("End of chromatogram.\n")
						
					}else{
						tmp1  		<-  zoomwidth[1] + (zoomwidth[2]-zoomwidth[1])/3
						tmp2  		<-  zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
						zoomwidth	<-	c(tmp1,tmp2)
						replot2(edges,DATA,zoomwidth,maintext)
					}
				
				}else if(zoom == 6){
					
					more <- "yes"
					while(more == "yes"){
						
						selected_edge	<-  getGraphicsEvent(prompt = "Click near an edge to select it", onMouseDown = function(buttons,x,y) which.min(abs(edges-grconvertX(x,from="nic",to="user"))))
						cat("Selected edge: ",selected_edge,"\n")
						abline(v=edges[selected_edge],col="black",lwd=2)
				 		x	<-	getGraphicsEvent(prompt = paste("Click where you want to move edge ",selected_edge,".",sep=""), onMouseDown = function(buttons,x,y) ifelse(grconvertX(x,from="nic",to="user")<0,1,ifelse(grconvertX(x,from="nic",to="user")>ncol(DATA),ncol(DATA),grconvertX(x,from="nic",to="user"))))
						oldedge					<-	edges[selected_edge]
            			edges[selected_edge]	<-	round(x)
            			edges					<-	sort(edges)
            			replot2(edges,DATA,zoomwidth,maintext)
						abline(v=oldedge,col="grey",lwd=2)
						more  <-  tk_messageBox(type = "yesno",message = "Move more edges?")
					}
				
				}else if(zoom == 7){
					tempzoomwidth <-  zoomwidth
					more  <-  "yes"
					
					while(more == "yes"){
						x			<-	getGraphicsEvent(prompt = "Click where you want to set an edge", onMouseDown = function(buttons,x,y) ifelse(grconvertX(x,from="nic",to="user") < 1,1,ifelse(grconvertX(x,from="nic",to="user")>ncol(DATA),ncol(DATA),grconvertX(x,from="nic",to="user"))))
						edges 		<-	sort(c(edges,round(x)))
   						replot2(edges,DATA,zoomwidth,maintext)
   						more	<-	tk_messageBox(type = "yesno",message = "Set more edges?")
					}
				}else if(zoom == 8){
					
					more <- "yes"
					while(more == "yes"){
						del_edge	<-	getGraphicsEvent(prompt = "Click near an edge to delete it", onMouseDown = function(buttons,x,y) which.min(abs(edges-grconvertX(x,from="nic",to="user"))))
				 		cat("\nSelected edge: ",del_edge,"\n\n")
				 		abline(v=edges[del_edge],col="black",lwd=2)
				 		if(tk_messageBox(type="yesno",message=paste("Delete selected edge? (Edge ",del_edge,")")) == "yes"){
							deletededges <-  c(deletededges,edges[del_edge])
							edges 	<-  edges[-del_edge]
       						replot2(edges,DATA,zoomwidth,maintext)
 				 	  		abline(v=deletededges,col="grey",lwd=2)
						
						}else
							abline(v=edges[del_edge],col="red",lwd=2)
							more  <-  tk_messageBox(type = "yesno",message = "Delete more edges?")
					}
						
					deletededges	<-numeric()
				}
	 		}
		}
	}
	
	if(length(edges)){
		cat("\nWindow information: \n\nNumber of windows: ", length(edges)-1,"\n")
		cat("Largest window:\t",max(diff(edges))," (win",which.max(diff(edges)),")\n")
 		cat("Smallest window: ", min(diff(edges)),"(win",which.min(diff(edges)),")\n")
		cat("Number of files: ", length(files),"\n\n")
 		#do_log(projectpath,paste(paste("Edges: ",edges,sep=""),"\n"))
		load(file.path(projectpath,"Aligned","SCAN_RANGE.Rdata"))
		load(files[which.min(shift)])
 		EDGES_TIME	<-	SCAN_INFO[,2]
 		save(EDGES_TIME,file=file.path(projectpath,"Edges","EDGES_TIME.Rdata"))
 		save(edges,min_win,max_win,file=file.path(projectpath,"Edges","edges.Rdata"))
	
	}else
		cat("No edges set!\n")
	graphics.off()
	
}

