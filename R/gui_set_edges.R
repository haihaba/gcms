##' Function gui_set_edges
##' 
##' Function gui_set_edges
##' @export
##' @param projectpath
gui_set_edges<-function(projectpath){
	
	dir.create(file.path(projectpath,"Edges"),showWarnings=FALSE)
		
	load(file.path(projectpath,"Aligned","files.Rdata"))
	load(file.path(projectpath,"Aligned","shift.Rdata"))
	load(file.path(projectpath,"maxMZ.Rdata"))
	min_win 			<-	120
	max_win 			<-	450
	edges				<-	numeric()
	deletededges		<-	numeric()
	menuvector  		<-	c(TRUE,FALSE,FALSE,FALSE,TRUE)
	menulist 			<-	c("Load data","Load old edges","Auto edges","Plot controls","Quit")
	k  					<-	1
	
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
				
        ### Here the chromatograms are plotted
				#for(i in 1:nrow(DATA)){
          averagedData<-apply(DATA,2,mean)
          plot(averagedData,type='l', xlim=c(0,length(averagedData)),ylim=c(0,max(averagedData)*1.01),main=maintext,xlab="",ylab="")        
					#plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim=c(0,ncol(DATA)),ylim=c(0,max(DATA)*1.01),main=maintext,xlab="",ylab="")
		  		par(new=TRUE)
 				#}
 				
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
				
			}
			
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
   					replot(edges,DATA,zoomwidth,maintext)
				
				}else if(zoom == 1){ #Zoom in
					cat("\n\nSelect two points on the x-axis to zoom in.\n\n")
					zoomwidth			<-	round(sort(locator(2,"p")$x))
					replot(edges,DATA,zoomwidth,maintext)
				
				}else if(zoom == 2){ #Zoom default
					zoomwidth	<-	c(0,0)
   					replot(edges,DATA,zoomwidth,maintext)
   					
				}else if(zoom == 3){ #Zoom out
					
					if(sum(zoomwidth) == 0){
				    	cat("You must zoom in before you zoom out!\n")
					
					}else{
						tmp1		<-	zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
						tmp2		<-	zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
					  	zoomwidth	<-	c(tmp1,tmp2)
					  	
						if(diff(zoomwidth) > ncol(DATA))
							zoomwidth <-  c(0,0)
   						replot(edges,DATA,zoomwidth,maintext)
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
						replot(edges,DATA,zoomwidth,maintext)
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
						replot(edges,DATA,zoomwidth,maintext)
					}
				
				}else if(zoom == 6){
					
					#more <- "yes"
					#while(more == "yes"){
						cat('Click near an edge to select it\n')
            selected_edge <- which.min(abs(edges-locator(1,type='p')[[1]]))
						#selected_edge	<-  getGraphicsEvent(prompt = "Click near an edge to select it", onMouseDown = function(buttons,x,y) which.min(abs(edges-grconvertX(x,from="nic",to="user"))))
						cat("Selected edge: ",selected_edge,"\n")
						
            abline(v=edges[selected_edge],col="black",lwd=2)
				 		cat('click where you want to move the edge')
            x <- locator(1,type='p')[[1]]
						x <- ifelse(x < 0, 1, ifelse(x > ncol(DATA), ncol(DATA), x))
            #x	<-	getGraphicsEvent(prompt = paste("Click where you want to move edge ",selected_edge,".",sep=""), onMouseDown = function(buttons,x,y) ifelse(grconvertX(x,from="nic",to="user")<0,1,ifelse(grconvertX(x,from="nic",to="user")>ncol(DATA),ncol(DATA),grconvertX(x,from="nic",to="user"))))
						      
            oldedge <- edges[selected_edge]
            edges[selected_edge] <- round(x)
            edges <- sort(edges)
            replot(edges,DATA,zoomwidth,maintext)
						abline(v=oldedge,col="grey",lwd=2)
						#more <- tk_messageBox(type = "yesno",message = "Move more edges?")
					#}
				
				}else if(zoom == 7){
					tempzoomwidth <-  zoomwidth
					#more  <-  "yes"
					
					#while(more == "yes"){
            cat('Click where you want to set an edge\n')
            x <- locator(n=1,type='p')[[1]]
            
            ## check if chosen value is in
            ## the range of the chromatorgram
            x<-ifelse(x<0,1,ifelse(x>ncol(DATA),ncol(DATA),x))
            #x			<-	getGraphicsEvent(prompt = "Click where you want to set an edge", onMouseDown = function(buttons,x,y) ifelse(grconvertX(x,from="nic",to="user") < 1,1,ifelse(grconvertX(x,from="nic",to="user")>ncol(DATA),ncol(DATA),grconvertX(x,from="nic",to="user"))))
						edges 		<-	sort(c(edges,round(x)))
   					replot(edges,DATA,zoomwidth,maintext)
   					#more	<-	tk_messageBox(type = "yesno",message = "Set more edges?")
					#}
          
				}else if(zoom == 8){
					
          
					#more <- "yes"
					#while(more == "yes"){
            cat("Click near an edge to delete it\n")
            del_edge <- which.min(abs(edges-locator(1,'p')[[1]]))
						#del_edge	<-	getGraphicsEvent(prompt = "Click near an edge to delete it", onMouseDown = function(buttons,x,y) which.min(abs(edges-grconvertX(x,from="nic",to="user"))))
				 		cat("\nSelected edge: ",del_edge,"\n\n")
				 		abline(v=edges[del_edge],col="black",lwd=2)
				 		#if(tk_messageBox(type="yesno",message=paste("Delete selected edge? (Edge ",del_edge,")")) == "yes"){
							deletededges <-  c(deletededges,edges[del_edge])
							edges 	<-  edges[-del_edge]
       						replot(edges,DATA,zoomwidth,maintext)
 				 	  		abline(v=deletededges,col="grey",lwd=2)
						
						#}else
						#	abline(v=edges[del_edge],col="red",lwd=2)
							#more  <-  tk_messageBox(type = "yesno",message = "Delete more edges?")
					#}
						
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
	#graphics.off()
	
}

##' A_E Function, Automatic Edges
##' 
##' Probably the funciton to set automatic edges for the processing windows
##' @param DATA
##' @param Nmin
##' @param Nmax
##' @return edges 
A_E <- function(DATA,Nmin,Nmax){
  edges  <-	c(1,length(DATA))
  N		<-	0.5
  while(max(diff(edges)) > Nmax){
    n		<-	which.max(diff(edges)) 
    range	<-	(edges[n] + Nmin):(edges[n+1] - Nmin)
    x		<-	E_peak_pick(DATA[range],median(DATA[range]),Nmin)
    
    if(!any(x$peak>0)){
      N	<-	N*1.5
      if(N > 7){
        n	<-	which.min(DATA[range])
        
        if(length(n) > 1){
          k	<-	round(length(n)/2)
          n	<-	n[k]
        }
        
        N			<-	0.5
        x$peak[n]	<-	1
      }
    }else
      N	<-	0.5
    edges	<-	sort(c(edges,range[which(x$peak>0)]))
  }
  
  edges <-  edges[!(edges == 1 | edges == length(DATA) | edges > max(which(DATA>0)) | edges < min(which(DATA>0)))]
}

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

##' Function E_peak_pick
##' 
##' Function E_peak_pick
##' @param x
##' @param NL
##' @param min_win
##' @return x 
E_peak_pick <- function(x,NL,min_win){
  require(signal)
  #xd1 <- -sav.gol(x,11,3,1)
  xd1 <-  sgolayfilt(x,3,11,1)
  N1	<-	which(x < NL)
  N2	<-  numeric(length(x))
  
  
  for(i in 3:(length(x)-2)){
    if(xd1[i-2] < 0  & xd1[i-1]<0 & xd1[i+1]>0 & xd1[i+2]>0 &  sum(N1 == i)==1)
      N2[i]	<- 1
  }
  
  N	<-	intersect(N1,which(N2 == 1))
  
  if(length(N)==0)
    return()
  
  if(length(N) != 1)
    while(min(diff(N)) < min_win){
      p1	<- 	which.min(diff(N))
      p2	<-	p1+1
      
      if(x[N[p1]] < x[N[p2]])
        N <- N[-p1]
      else
        N <- N[-p2]
      
      if(length(N) == 1)
        break
    }
  
  xpeak		<-	numeric(length(x))
  xpeak[N]	<-	x[N]
  x			<-	list(peak = xpeak, out = x)
  
}


##' Function sgolayfilt
##' 
##' Function sgolayfilt
##' @param x
##' @param p
##' @param n
##' @param m
##' @param ts
sgolayfilt<-function(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
  
  len = length(x)
  if (class(p) == "sgolayFilter" || (!is.null(dim(p)) && dim(p) > 1)){
    F = p
    n = nrow(F)
  }else
    F = sgolay(p, n, m, ts)
  k = floor(n/2)
  
  #z = filter(F[k+1,n:1], 1, x)
  
  filt  <-  F[k+1,n:1]
  z 		<-	na.omit(stats:::filter(c(rep(0,length(filt) - 1), x), filt, sides = 1))
  c(F[1:k,] %*% x[1:n], z[n:len], F[(k+2):n,] %*% x[(len-n+1):len])
}

##' Function sgolay
##' 
##' Function sgolay
##' @param p
##' @param n
##' @param ts
##' @return F 
sgolay<-function(p, n, m = 0, ts = 1){ 
  
  library(MASS)
  if (n %% 2 != 1)
    stop("sgolay needs an odd filter length n")
  
  if (p >= n)
    stop("sgolay needs filter length n larger than polynomial order p")
  
  F = matrix(0., n, n)
  k = floor(n/2)
  for (row  in  1:(k+1)) {
    C = ( ((1:n)-row) %*% matrix(1, 1, p+1) ) ^ ( matrix(1, n) %*% (0:p) )
    A = ginv(C)
    F[row,] = A[1+m,]
  }
  
  F[(k+2):n,] = (-1)^m * F[k:1,n:1]
  
  if (m > 0)
    F = F * prod(1:m) / (ts^m)
  
  class(F) = "sgolayFilter"
  F
}


