##' Function MCR_settings
##' 
##' Function MCR_settings
##' @export
##' @param projectpath
##' @return MCR
MCR_settings <-function(projectpath){
	
	if(file.exists(file.path(projectpath,"Edges","edges.Rdata"))){
		load(file.path(projectpath,"Edges","edges.Rdata"))
	
	}else{
		cat("Error! You must set edges first!\n")
		Sys.sleep(1)
		return(NULL)	
	}
	
	load(file.path(projectpath,"Aligned","files.Rdata"))
	filestemp       <-  sub("(.+)[.][^.]+$", "\\1", basename(files))
	DATA  <-  list()
	
	if(file.exists(file.path(projectpath,"HMCR","MCR.Rdata"))){  # Load and check old settings
		cat("Loading old settings..\n")
		load(file.path(projectpath,"HMCR","MCR.Rdata"))
		
		if(is.null(MCR$edges))
			MCR$edges <-  0
		
 		if(!all(edges == MCR$edges)){
			cat("Warning! Mismatch between old window edges and current window edges. Selected windows for processing will be reset!\n")
			MCR$windowlist	<-	c(paste("Window[00",1:9,"] [P]",sep=""),paste("Window[0",10:99,"] [P]",sep=""),paste("Window[",100:999,"] [P]",sep=""))[1:(length(edges)-1)]
			unlink(file.path(projectpath,"Edges","dat"),recursive=TRUE)

			if(menu(c("Yes","No"),title="Reset to default settings as well?") == 1)
				file.remove(file.path(projectpath,"HMCR","MCR.Rdata"))
		}
			
		if(length(list.files(file.path(projectpath,"Edges","dat")))){
			cat("Warning! Window data folder is not empty. The content of this folder must be removed if any changes are made regarding the samples used when processing!\n")
 		if(menu(c("Yes","No"),title="Remove content in window data folder?") == 1)
			unlink(file.path(projectpath,"Edges","dat"),recursive=TRUE)
		}
	}
		
	if(!file.exists(file.path(projectpath,"HMCR","MCR.Rdata"))){  ## MCR default values ##
		MCR	<-  list()
		MCR$reg$incl    <-  MCR$cvD$incl    <-  1:length(files)
		MCR$reg$excl	<-	MCR$reg$pred    <- numeric()
		MCR$cp   		<-  0.9
		MCR$sp  		<-  0.95
		MCR$edges 		<-  edges
		MCR$windowlist	<-	c(paste("Window[00",1:9,"] [P]",sep=""),paste("Window[0",10:99,"] [P]",sep=""),paste("Window[",100:999,"] [P]",sep=""))[1:(length(edges)-1)]
	}
	
	samproc			<-	" All samples"
	numfiles		<-	length(files)
	incl  			<-	1:length(files)
	pred			<-	numeric()
	excl			<-	numeric()
	datamenu		<-	character()

	if(file.exists(file.path(projectpath,"Aligned","SUM_MZ.Rdata")))
		datamenu	<-	c(datamenu,"Tot Data")
	if(file.exists(file.path(projectpath,"Aligned","BASEPEAK.Rdata")))
		datamenu	<-	c(datamenu,"BASEPEAK")
	if(file.exists(file.path(projectpath,"Aligned","TIC.Rdata")))
		datamenu	<-	c(datamenu,"TIC")
	
	dsource 		<-	scaling	<-	"None"
	k				<-	1
	menulist		<-	c("Data source","Scaling","Regular method settings","Cross-validation method settings","Select windows to process","Quit")
	menuvector		<-	c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE)

	while(k){
		cat("==========================\n")
		cat("====== MCR Settings ======\n")
		cat("==========================\n\n")
		cat("Windows to process:\t",grep("[P]",MCR$windowlist),"\n")
		cat("Number of files:\t",numfiles,"\n")
		cat("Datasource:\t\t",dsource,"\n")
		cat("Scaling:\t\t",scaling,"\n")

		k	<-	menu(choices=menulist[menuvector],title="")

		if(k == sum(menuvector) | k == 0){
			k <-  0
			
		}else if(menulist[menuvector][k] == "Data source"){
			value <-  menu(datamenu,title="Choose a data source")
				
			if(value){
				dsource	<-	datamenu[value]
					
				if(dsource == "Tot Data"){
					load(file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
					DATA$data	<-	DATA$unalt  <-  SUM_MZ
					
				}else if(dsource == "BASEPEAK"){
					load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
					DATA$data  <-  DATA$unalt  <-  BASEPEAK
					
				}else if(dsource == "TIC"){
					load(file.path(projectpath,"Aligned","TIC.Rdata"))
					DATA$data  <-  DATA$unalt  <-  TIC
					
				}
									
				scaling <-  "None"
				maintext=paste(c("PCA plot",paste("Scaling: ",scaling),paste("Data: ",dsource)))
				pc  <-  pca(DATA$data,2)
					
				replot_mcr(x=pc$vec[,1],y=pc$vec[,2],maintext,incl,excl,numfiles,plotlabels=FALSE,pred=pred)
				
				menuvector[c(2,3,5)] <-  TRUE
					
			}
			
		}else if(menulist[menuvector][k] == "Scaling"){
			value  <-  1
			while(value){
				scalingmenu <-	c("Center","UV-scale","Pareto","None","Done")
				value <-  menu(scalingmenu,title="Choose scaling method")
				if(value != 5 & value){
					scaling <-  scalingmenu[value]
					if(value == 1)
						DATA$data	<-	cent(DATA$unalt)
					if(value == 2)
						DATA$data	<-	uv_sc(DATA$unalt)
					if(value == 3)
						DATA$data	<-	par_sc(DATA$unalt)
					if(value == 4)
						DATA$data <-  DATA$unalt
							
					maintext	<-	paste(c("PCA plot",paste("Scaling: ",scaling),paste("Data: ",dsource)))
					pc  <-  pca(DATA$data,2)
					replot_mcr(x=pc$vec[,1],y=pc$vec[,2],maintext,incl,excl,numfiles,plotlabels=FALSE,pred=pred)
					
						
				}
					
				value <-  !((value == 5) | !value)
			}
		}else if(menulist[menuvector][k] == "Regular method settings"){
			reglist		<-	c("Auto-select","Selection controls","Done")
			regvec		<-	!logical(3)
			value		<-	1
			filelist	<-	MCR$k$files #DATA$k$files
		  		
			while(value){
				cat("===================================\n")
				cat("===== Regular method settings =====\n")
				cat("===================================\n\n")
				cat("Number of files:\t",numfiles,"\n")
				cat("Included files:\t\t",length(MCR$reg$incl),"\n")
				cat("Prediction files:\t",length(MCR$reg$pred),"\n")
				cat("Excluded files:\t\t",length(MCR$reg$excl),"\n\n")
		  			
				value <-  menu(reglist[regvec],title="")
		
				if(value == sum(regvec) | !value){
					value <-  0
				}else if(reglist[regvec][value] == "Auto-select"){
					while(1){
						num_samples	<-	as.integer(readline(prompt = paste("Input number of model samples (integer from 5 to ",numfiles,"): ")))
						ifelse(!is.na(num_samples) & is.integer(num_samples) & num_samples >= 5 & num_samples <= numfiles,break,"")
						cat("Invalid input!\n")
					}
		
					if(num_samples < numfiles){
						sfd		<-	sfd_par(pc$vec,num_samples)
						incl	<-	sort(sfd$object)
						
					}else  #if num_samples == numfiles
						incl  <-  1:numfiles
 							
 					if(num_samples < numfiles){
 						pred  <-  numeric()
 		
 						while(1){
 							num	<-	as.integer(readline(paste("Input number of prediction samples (integer from 0 to: ",(numfiles-num_samples),"): ")))
							ifelse(!is.na(num) & is.integer(num) & num >= 0 & num <= (numfiles-num_samples),break,"")
							cat("Invalid input!\n")
						}
								
						if(num){
							if(num+num_samples == numfiles){
								pred  <-  (1:numfiles)[-incl]
									
							}else{
								sfd   <-  sfd_par(pc$vec,num,exclude=incl)
								pred 	<-  sort(sfd$object)
							}
						}
					}
						
					excl  <-  (1:numfiles)[-c(incl,pred)]
							
					filelist  <-  filestemp   ###MCR$reg$filelist TRYING TO FIX ORIGNAL CODE
					filelist[1:length(filelist) %in% incl]	  <-  paste("[incl] ",filestemp[1:length(filelist) %in% incl],sep="")
					filelist[1:length(filelist) %in% pred]	  <-  paste("[pred] ",filestemp[1:length(filelist) %in% pred],sep="")
					filelist[1:length(filelist) %in% excl]	  <-  paste("[excl] ",filestemp[1:length(filelist) %in% excl],sep="")
						
					MCR$reg$incl	<-	grep("[incl]",filelist,fixed=TRUE)
					MCR$reg$pred	<-	grep("[pred]",filelist,fixed=TRUE)
					MCR$reg$excl	<-	grep("[excl]",filelist,fixed=TRUE)
					MCR$k$files		<-	MCR$reg$filelist  <-  filelist
							
					replot_mcr(pc$vec[,1],pc$vec[,2],maintext,incl=MCR$reg$incl,excl=MCR$reg$excl,numfiles,plotlabels=FALSE,pred=MCR$reg$pred)
						
				}else if(reglist[regvec][value] == "Selection controls"){
					sellist <-  c("Include model files","Include prediction files","Exclude files","List files","Done")
					selvec  <-  !logical(5)
					choice  <-  1
					filelist  <-  MCR$reg$filelist
						
					while(choice){
						cat("Number of files:\t",numfiles,"\n")
						cat("Included files:\t\t",length(MCR$reg$incl),"\n")
						cat("Prediction files:\t",length(MCR$reg$pred),"\n")
						cat("Excluded files:\t\t",length(MCR$reg$excl),"\n\n")
						choice  <-  menu(sellist[selvec],title="")
								
						if(choice == sum(selvec) | choice == 0){
							choice  <-  0
							
						}else if(sellist[selvec][choice] == "Include model files"){
							incl	<-	select.list(filelist,multiple=TRUE,title="Choose model files to include")
							filelist[filelist %in% incl]	<-	paste("[incl] ",filestemp[filelist %in% incl],sep="")
								
						}else if(sellist[selvec][choice] == "Include prediction files"){
							pred	<-	select.list(filelist,multiple=TRUE,title="Choose prediction files to include")
							filelist[filelist %in% pred]	<-	paste("[pred] ",filestemp[filelist %in% pred],sep="")
						}else if(sellist[selvec][choice] == "Exclude files"){
							excl	<-	select.list(filelist,multiple=TRUE,title="Choose files to exclude")
							filelist[filelist %in% excl]	<-	paste("[excl] ",filestemp[filelist %in% excl],sep="")
						}else if(sellist[selvec][choice] == "List files")
							select.list(filelist,title = "Samples")
								
						replot_mcr(pc$vec[,1],pc$vec[,2],maintext,incl=grep("incl",filelist),excl=grep("excl",filelist),numfiles,plotlabels=FALSE,pred=grep("pred",filelist))

						MCR$reg$incl		<-	grep("[incl]",filelist,fixed=TRUE)
						MCR$reg$pred		<-	grep("[pred]",filelist,fixed=TRUE)
						MCR$reg$excl		<-	grep("[excl]",filelist,fixed=TRUE)
						MCR$reg$filelist	<-	MCR$k$filelist	<-  filelist
					}
				}
			}
				
		}
			
		else if(menulist[menuvector][k] == "Select windows to process"){
				
			## Select which windows to process
			winmenu 		<-  c("Add windows","Remove windows","Next window","Previous window","Add this window","Remove this window","Jump to window","Change Data Source","Back")
			choice		<-	num	<-	1
			windowlist	<-	MCR$windowlist
			tmp			<-	substr(MCR$windowlist,1,11)
		
			load(file.path(projectpath,"Aligned","TIC.Rdata"))
			windowDATA	<-	TIC
			maintext		<-	"Total ion current"
			temp			<-	0.1*(edges[num+1]-edges[num])
			replot2(edges,windowDATA,zoomwidth=c(edges[num]-temp,edges[num+1]+temp),maintext)
			title(c("","",paste(windowlist[num])))
		
			while(choice){
				temp	<-	0.1*(edges[num+1]-edges[num])
				replot2(edges,windowDATA,zoomwidth=c(edges[num]-temp,edges[num+1]+temp),maintext)
				title(c("","",paste(windowlist[num])))
				choice  <-  menu(winmenu,title="Select windows to process")
			
				if(!choice | choice == 9){
					choice  <-  0
			
				}else if(choice == 1){ #Add windows
					mark  <-  select.list(windowlist,preselect=windowlist[1],multiple=TRUE,title="Select windows to process")
					windowlist[windowlist %in% mark] <-  paste(tmp[windowlist %in% mark]," [P]")
					cat("Window added!\n")
			
				}else if(choice == 2){ #Remove windows
					unmark  <-  select.list(windowlist,preselect=windowlist[1],multiple = TRUE,title = "Remove windows")
					windowlist[windowlist %in% unmark] <-  tmp[windowlist %in% unmark]
					cat("Window removed!\n")
				
				}else if(choice == 3){ #Next window
					num <-  min(length(edges)-1,num+1)
				
				}else if(choice == 4){ #Previous window
					num <-  max(1,num-1)
				
				}else if(choice == 5){ #Add this window
					windowlist[1:(length(edges)-1) %in% num] <-  paste(tmp[1:(length(edges)-1) %in% num]," [P]")
				
				}else if(choice == 6){ #Remove this window
					windowlist[1:(length(edges)-1) %in% num] <-  tmp[1:(length(edges)-1) %in% num]
				
				}else if(choice == 7){ #Jump to window
					numtemp  <-  as.numeric(substr(select.list(windowlist,preselect=windowlist[1],multiple = FALSE,title = "Jump to window"),8,10))
					
					if(is.numeric(numtemp))
						num <-  numtemp
				
				}else if(choice == 8){ #Change Data Source
					datasource  <-  menu(choices=c("Total ion Current","Basepeak Chromatogram"),title = "Choose type of data to load")
					if(datasource){
						if(datasource == 1){
							load(file.path(projectpath,"Aligned","TIC.Rdata"))
							windowDATA	<-	TIC
							maintext		<-	"Total ion current"
							rm(TIC)
						}else if(datasource == 2){
							load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
							windowDATA  		<-  BASEPEAK
							maintext  <-  "Basepeak Chromatogram"
							rm(BASEPEAK)
						}
					}
				}
			}
	
			MCR$windowlist  <-  windowlist
		}

	}
	graphics.off()
	dir.create(file.path(projectpath,"HMCR"),showWarnings=FALSE)
	save(MCR=MCR,file=file.path(projectpath,"HMCR","MCR.Rdata"))
	return(MCR)
	
}


##' Function sfd_par
##' 
##' Function sfd_par
##' @param X
##' @param num
##' @param exclude
##' @return object, j, md
sfd_par<-function(X,num,exclude = nrow(as.matrix(X))+1){
  
  N				<-  nrow(as.matrix(X))
  K				<-  ncol(as.matrix(X))
  
  if(num == (N-length(exclude))){
    sel <-  (1:N)[-exclude]
    return(list(object=sel))
  }
  
  
  D			<-	obj_eu(X)
  diag(D)		<-	2*max(D)
  MD			<-	j	<-	0
  
  for(run in 1:500){
    object	<-	sample((1:N)[-exclude])[1:num]
    OK			<-	3
    md			<-	0
    number	<-	1
    
    while(OK){
      XD	<-	D[object,object]
      OK	<-	OK-1
      nk	<-	which.min(XD)
      k	<-	col(XD)[nk]
      
      if(OK)
        n	<-	row(XD)[nk]*(OK != 2) + k*(OK == 2)
      else
        n	<-	number
      
      number 		<-  number + (OK==0)
      
      if(!OK)
        OK	<-	(number <= num)
      
      obs			<-	object[-n]
      ex_object	<-	setdiff((1:N)[-exclude],obs)
      I			<-	which.max(apply(D[obs,ex_object],2,min))
      
      for(jj in I){
        obs[num]	<-	ex_object[jj]
        
        if((md <= min(D[obs,obs])) & sum(abs(sort(obs)-sort(object))) > 0){
          XD				<-	D[obs,obs]
          diag(XD)		<-	0
          
          if(md == min(D[obs,obs])){
            
            if(j < sum(XD)/2){
              md			<-	min(D[obs,obs])
              object		<-	sort(obs)
              j			<-	sum(XD)/2
              number		<-	number*(OK != 1) + (OK == 1)
              
              if(OK>1)
                OK	<-	3
            }
          }else{
            md		<-	min(D[obs,obs])
            object	<-	sort(obs)	
            j		<-	sum(XD)/2
            number	<-	number*(OK != 1) + (OK == 1)
            
            if(OK>1)
              OK	<-	3
          }
        }
      }
    }
    
    XD				<-	D[object,object]
    md				<-	min(D[object,object])
    
    diag(XD)	<-	0
    j			<-	sum(XD)/2
    object		<-	sort(object)
    
    if(MD < md){
      J	<-	j
      sel	<-	object
      MD	<-	md
    }
  }
  list(object=sel,j=J,md=MD)
}

##' Function X
##' 
##' Function X
##' @param X
##' @return D
obj_eu<-function(X)
  D <-  as.matrix(dist(X,upper=TRUE,diag=TRUE))


##' Function replot_mcr
##' 
##' Function replot_mcr
##' @param x
##' @param y
##' @param maintext
##' @param incl
##' @param excl
##' @param numfiles
##' @param set1
##' @param set2
##' @param seltype
##' @param plotlabels
##' @param pred 
replot_mcr<-function(x,y,maintext,incl,excl,numfiles,set1,set2,seltype=(!missing(set1) & !missing(set2)),plotlabels = TRUE,pred){
  
  if(seltype == 1){
    color <-  2*(1:numfiles %in% set1)+4*(1:numfiles %in% set2)+3*(1:numfiles %in% pred)+(1:numfiles %in% setdiff(1:numfiles,c(set1,set2,pred)))
    
  }else{
    color <-  2*(1:numfiles %in% incl)+(1:numfiles %in% excl)+3*(1:numfiles %in% pred)
  }
  
  plot(x,y,xlim=c(range(x)),ylim=range(y),main=maintext,xlab="t[1]",ylab="t[2]",col=color,cex=0.9,pch=16)
  
  if(plotlabels)
    text(x,y,pos=4,cex=0.7)
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



