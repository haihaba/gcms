##' Function check_data
##' 
##' Function check_data
##' @export
##' @param projectpath
check_data <- function(projectpath){
  
  
  ## Loading the list of files and setting variables
	load(file.path(projectpath,"Aligned","files.Rdata"))
  OBS			<-	1:length(files)
	select	<-	character()
	COMP  		<-	1:5
	# default setting, no scaling
  scalenum	<-	0
	x       	<-  y	<-  numeric()

  ## load the data
	load(file.path(projectpath,"Aligned","TIC.Rdata"))
 	
  ## Use TIC for PCA, scale date, calc PCA
 	Xin				<-	TIC[OBS,]
 	maintext		<-	"Total Ion Current"
  	scale 			<-	scaling(Xin,scalenum)
  	X     			<-	scale$Xout
 	scale_text		<-	scale$scale_text
 	vec				<-	pca(X,2)$vec
	
  ## Set up folder for Check Data
	dir.create(file.path(projectpath,"Check Data"),showWarnings = FALSE)
	
	for(i in 1:length(files)){
		if(i < 100){
			if(i > 9)
				select[i]	<-	paste("[0", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
      		else
				select[i]	<-	paste("[00", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
		}
		else
			select[i]	<-	paste("[", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
  	}
  		
 		selectTemp <-  select

	outliers 	<-  character()
	k			<- 1
        #X11.options(type="cairo")
        #X11()

        replot=FALSE

        while(k!="Quit"){
          if(replot){
            menuselector<-c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
          }else{
            menuselector<-c(TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
          }

          menulist<-c("Data source","Zoom in","Standard zoom","Scaling","Export to excel","List files","Mark outliers","Unmark outliers","Quit")
          k<-menu(choices=menulist[menuselector],title="Check Data")
          k<-menulist[menuselector][k]

          if(k == "Data source"){
  				a	<-	menu(choices=c("Total Ion Current","Base-Peak Chromatogram","Total Mass-Spectrum","Done"),title="Choose data source")
    		
    			if(a == 1){
    				load(file.path(projectpath,"Aligned","TIC.Rdata"))
      				Xin				<-	TIC[OBS,]
      				maintext	<-	"Total Ion Current"
			
			}else if(a == 2){
     				load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
      				Xin			<-	BASEPEAK[OBS,]
				maintext	<-	"Basepeak Chromatogram"
			}else if(a == 3){
				load(file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
		      	Xin			<-	SUM_MZ[OBS,]
      				maintext	<-	"Total Mass-Spectrum"
			}
			
  	  			scale 			<-	scaling(Xin,scalenum)
  	  			X     			<-	scale$Xout
    			scale_text		<-	scale$scale_text
    			vec				<-	pca(X,2)$vec

			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select)
                        replot<-TRUE
                                
		}else if(k == "Scaling"){
			scalenum   <- menu(c("Center","UV-scaled","Pareto-scaled","None"),title="Choose scaling method")
			scale      <- scaling(Xin,scalenum)
			X          <- scale$Xout
			scale_text <- scale$scale_text
			vec        <- pca(X,2)$vec
 			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select)
                        replot<-TRUE
  	
                }else if(k == "Export to excel"){
                  write.table(Xin,file=file.path(projectpath,"Check Data",paste(maintext,".txt",sep="")),sep="\t",row.names = selectTemp)
                  cat("Data exported to: ",file.path(projectpath,"Check Data",paste(maintext,".txt\n",sep=""),fsep="\\"))
	
		}else if(k == "List files"){
			select.list(select,title="Outliers")
		
		}else if(k == "Mark outliers"){
      
      ## the data points in mark are excluded
		  require(sp)
      coordinates <- locator(type = "l")
		  coordinates$x[length(coordinates$x + 1)] <- coordinates$x[1]
		  coordinates$y[length(coordinates$y + 1)] <- coordinates$y[1]
		  new.poly <- Polygon(coordinates)
		  in.out <- as.logical(point.in.polygon(vec[,1], vec[,2], new.poly@coords[, 1], new.poly@coords[, 2]))
		      
			#mark  <-  select.list(select,preselect=select[1],multiple=TRUE,title="Mark outliers to remove.")
			select[in.out] <-  paste(selectTemp[in.out]," [OUT]")
			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select)
                        replot<-TRUE
                        
		}else if(k == "Unmark outliers"){
      
		  require(sp)
		  coordinates <- locator(type = "l")
		  coordinates$x[length(coordinates$x + 1)] <- coordinates$x[1]
		  coordinates$y[length(coordinates$y + 1)] <- coordinates$y[1]
		  new.poly <- Polygon(coordinates)
		  in.out <- as.logical(point.in.polygon(vec[,1], vec[,2], new.poly@coords[, 1], new.poly@coords[, 2]))
      
			#unmark  <-  select.list(select,preselect=s[1],multiple = TRUE,title = "Unmark outliers.")
			select[in.out] <-  selectTemp[in.out]
			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select)
                        replot<-TRUE
                }else if(k == "Zoom in"){
                        coords<-locator(n=2,type="p")
                        replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select,coords)
		}else if(k == "Standard zoom"){
                        replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),select)
                }
	}
	#graphics.off()
	remove_outliers(projectpath, grep("[OUT]",select,fixed=TRUE))
}

##' Function replot_checkdata
##' 
##' Function replot_checkdata
##' @param x
##' @param y
##' @param scale_text
##' @param maintext
##' @param OBS
##' @param ssX
##' @param select
##' @param coords
replot_checkdata<-function(x,y,scale_text,maintext,OBS,ssX,select,coords=NULL){
  
  xtext <-  paste("Component 1 (",round(t(x)%*%x/ssX*100,1),"% of variation in X)",sep="")
  ytext <-  paste("Component 2 (",round(t(y)%*%y/ssX*100,1),"% of variation in X)",sep="")
  
  color <-  select %in% select[grep("[OUT]",select)]+1
  if (missing(coords)){
    plot(x,y,main=c(maintext,scale_text),xlab=xtext,ylab=ytext,col=color,pch=16)
    text(x,y,OBS,pos=4,col=color)
  }else{
    plot(x,y,main=c(maintext,scale_text),xlab=xtext,ylab=ytext,col=color,pch=16,xlim=coords$x,ylim=coords$y)
    text(x,y,OBS,pos=4,col=color)
  }
}




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



##' Function pca
##' 
##' Function pca
##' @export
##' @param X
##' @param comp
##' @return vec, p
pca<-function(X,comp){
  
  #   Calculates Principal componets of X by calc eigvectors of X'*X or X*X'
  #   Depending on whats easiest to calculate....
  
  if(nrow(as.matrix(X)) < ncol(as.matrix(X))){
    tX  <-  t(X)
    p 	<-  eigen(X%*%tX)$vectors[,1:comp]
    
    if(comp>1){
      p <-  apply(p,2,function(p) tX%*%p)
      p <-  apply(p,2,function(p) p/sqrt(sum(p^2)))
      
    }else{
      p	<-	tX%*%p
      p 	<-	p/sqrt(sum(p^2))
    }
    
  }else
    p		<-	eigen(t(X)%*%X)$vectors[,1:comp]
  
  vec		<-	X%*%p
  vecp	<-	list(vec=vec,p=p)
}

##' Function remove_outliers
##' 
##' Function remove_outliers
##' @param projectpath
##' @param outliers
remove_outliers<-function(projectpath,outliers){
  
  OK	<-	menu(c("Yes","No"),title="Remove outliers from further analysis?")
  if(OK == 1 & length(outliers) & is.numeric(outliers)){
    
    cat("Removing outliers...\n\n")
    
    load(file.path(projectpath,"Aligned","TIC.Rdata"))
    load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
    load(file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
    load(file.path(projectpath,"Aligned","shift.Rdata"))
    load(file.path(projectpath,"Aligned","files.Rdata"))
    load(file.path(projectpath,"Aligned","COLUMNID1.Rdata"))
    
    
    files			<-	files[-outliers]
    shift			<-	shift[-outliers]
    BASEPEAK	<-	BASEPEAK[-outliers,]
    TIC				<-	TIC[-outliers,]
    SUM_MZ		<-	SUM_MZ[-outliers,]
    COLUMNID1	<-	COLUMNID1[-outliers]
    
    save(TIC,file=file.path(projectpath,"Aligned","TIC.Rdata"))
    save(BASEPEAK,file=file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
    save(SUM_MZ,file=file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
    save(shift,file=file.path(projectpath,"Aligned","shift.Rdata"))
    save(files,file=file.path(projectpath,"Aligned","files.Rdata"))
    save(COLUMNID1,file=file.path(projectpath,"Aligned","COLUMNID1.Rdata"))
    
  }else
    cat("Aborted by user or no outliers selected!\n")
  
}

##' Function ss
##' 
##' Function ss
##' @param X
##' @return SSX
ss<-function(X)
  SSX	<-	sum(X^2)


