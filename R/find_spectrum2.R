##' Function find_spectrum2
##' 
##' Function find_spectrum2
##' @param predpath
##' @param projectpath
##' @return type, windwolist
find_spectrum2<-function(predpath,projectpath){
	
	require(MASS)
	#X11.options(type="cairo")
	datamenu  <-  character()
	
	if(file.exists(file.path(projectpath,"HMCR","REG","MVA_DATA.Rdata")))
		datamenu  <-  "REG H-MCR DATA"
	
	
	if(length(datamenu)){
		dataexport  <-  menu(c(datamenu,"Cancel"),title="Select data")
		if(dataexport == 0 | dataexport == length(c(datamenu,"Cancel")))
  			datamenu <- character()
	}
	
	if(!length(datamenu)){
		cat("Error! (Aborted by user or no data found)\n\n")
		return(character())
	
	}else{
		type = ifelse(length(datamenu) == 1,strsplit(datamenu," ")[[1]][1],c("REG","CV")[dataexport])
		load(file.path(predpath,"Aligned","files.Rdata"))
		load(file.path(predpath,"Edges","edges.Rdata"))
		load(file.path(predpath,"Aligned","shift.Rdata"))
		load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
		load(file.path(predpath,"maxMZ.Rdata"))
		load(files[which.min(shift)])
		load(file.path(projectpath,"HMCR","MCR.Rdata"))
		
		windowlist		<-	as.numeric(substr(basename(list.files(file.path(projectpath,"HMCR",type,"win"))),4,6))
		EDGES_TIME		<-	SCAN_INFO[,2]
		
		load(file.path(projectpath,"SETTINGS.Rdata"))
		
		NL        <- SETTINGS$NL
		RT_LIMIT  <- SETTINGS$MPS
		DO_BL     <- SETTINGS$BC2
	 	#color     <- cbind("red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4")
    
    rm(Xbc,SETTINGS)
    
    dir.create(file.path(predpath,"Edges","dat"),recursive = TRUE,showWarnings = FALSE)
    temp <- list.files(file.path(projectpath,"Edges","dat"),full.names=TRUE)
    dir.create(file.path(predpath,"Edges","dat","Model samples_bg_corr"),showWarnings = FALSE,recursive=TRUE)
    dir.create(file.path(predpath,"HMCR",type,"win_png"),showWarnings = FALSE,recursive = TRUE)
		dir.create(file.path(predpath,"HMCR",type,"win"),showWarnings = FALSE,recursive = TRUE)
	  
    gc()
		
    Scores  <-  numeric()
		
		### attpempt to store the background correction data. bg correction storage in "find_spectrum.R" is done in arrays per window 
		### scan * m/z * sample. Here processing is sample wise (in "find_spectrum.R", the whole sampleset is in memory at once),
		### so the array has to be build step by step
		
		### First array variables for each window has to be defined. Then in the two main loops "sample" and "window" the variables
		### are filled sample by sample. After finishing the two loops, the arrays have to be stored in corresponding files 
		
		
		
		#for(win in windowlist){
		#	scans<-(edges[win+1]-edges[win]+1)
		#	mzs<-length(SCAN_RANGE)
		#	BL<-array(0,dim=c(scans,mzs,length(files)))
		#	#eval(parse(text=paste("BL_",win,"<-NULL",sep="")))
		#	save(BL,file=file.path(predpath,"Edges","dat","Model samples_bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
		#	rm(BL)
			
		#}
		
		if(length(windowlist)){
	 		for(i in 1:length(files)){
	    		load(files[i])
	    		cat("Resolving ",basename(files[i]),paste(" (",i,"/",length(files),")\n",sep=""))
		      	
		      	for(win in windowlist){
		      		
		      		if(edges[win]+shift[i] > 0 & edges[win+1]+shift[i] < nrow(Xbc)){
		      			x <-  Xbc[edges[win]:edges[win+1]+shift[i],]
						
						if(ncol(x) < maxMZ)
							x <-  cbind(x,matrix(0,nrow=nrow(x),ncol= maxMZ-ncol(x)))
	            			
	            		if(DO_BL){ # Removes baseline for prediction files
	      					BL<-matrix(apply(x,2,function(x) approx(c(1,length(x)),x[c(1,length(x))],1:length(x))$y),nrow=length(edges[win]:edges[win+1]),ncol=maxMZ)       # method = 'linear' by default
	        				x		<-	x-BL
	        				x[x<0]	<-	0
	        	
	        				###test code###
	        				#BL_new<-BL
	        				#load(file.path(predpath,"Edges","dat","Model\ samples","bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
	        				#BL[,,i]<-BL_new
	        				#save(BL,file=file.path(predpath,"Edges","dat","Model samples","bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
	        					
	        				#rm(BL_new)
	        				rm(BL)
	        			}
							
						load(file.path(projectpath,"HMCR", "REG","win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
						 # browser()
						if(!length(S))
							c2 <-  numeric()
						else
							c2 <-  do_AR_all_prediction(x,S,CP,RT_LIMIT)
						if(exists("PCApara"))
							scores  <-  (colSums(x)  <-  PCApara$mean)%*%PCApara$loadings[,1:2]
						else
							scores  <-  numeric()
							ok <- 1
						
					}else{
						load(file.path(projectpath,"HMCR", "REG","win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
						c2 <-  matrix(0,length(edges[win]:edges[win+1]),ncol(as.matrix(C)))
						ok  <-  0
						if(exists("PCApara"))
							scores <-  c(0,0)
						else
							scores  <-  numeric()
					}
						
						
					if(i == 1){
						if(ok){
							C<-matrix(colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))),nrow=1)
						    	CC  <-  c2
						}else{
							C<-matrix(colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))),nrow=1)*NA
							CC<-c2
						}
						
						Scores  <-  rbind(Scores,scores)
	    					save(C,CC,S,TIME,Scores,file=file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
					}else{
						load(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
						Scores  <-  rbind(Scores,scores)
						if(ok){
							C<-rbind(C,colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))))
							CC<-rbind(CC,c2)
						}else{
							C<-rbind(C,colSums(as.matrix(c2,ncol=ncol(as.matrix(C))))*NA)
							CC<-rbind(CC,c2)
						}
	    					save(C,CC,S,TIME,Scores,file=file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
	      			}
	      			gc()
	      		}
				gc()
			}
			
			# Plot windows
			# for(win in windowlist){
				
				# if(win %in% as.numeric(substr(basename(list.files(file.path(projectpath,"HMCR",type,"win"))),4,6))){
					# rm(Scores)
					# load(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
			
					# if(min(dim(CC)) > 0){
						# if(length(Scores)){
							   # #Plot functions for Cross validation method
						# }
			      		
			      		# Hz		<-	1/median(diff(SCAN_INFO[,2]))
						
						# for(i in 1:ncol(CC)){
							# c2 <-  matrix(CC[,i],length(edges[win]:edges[win+1]),length(files))
							# if(length(EDGES_TIME) < edges[win+1])
								# EDGES_TIME <-  c(EDGES_TIME,approxExtrap(1:length(EDGES_TIME),EDGES_TIME,xout=(length(EDGES_TIME)+1):edges[win+1])$y)
				          	# for(j in 1:ncol(c2)){
				        			# plot(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz,c2[,j],col=color[i],xlab="Time",ylab="Intensity",type="l",ylim=cbind(0,max(CC)),xlim=range(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz),main=paste("Window:",win),lwd=1)
				   				# par(new=TRUE)
							# }
				     	# }
						
						# par(new=FALSE)
			    		# #savePlot(file.path(predpath,"HMCR", "REG","win_png",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),sep="")),type="png")
					
					# }else
						# cat("No data found in window ", win,"!\n")
				# }
			# }
		
		}else
			cat("No windows processed.\n")
		
		return(list(type=type,windowlist=windowlist))
	}

}

