separate_windows <-
function(projectpath,windowlist,files,edges,model,prediction){
	
	if(length(windowlist)){
		load(file.path(projectpath,"maxMZ.Rdata"))
		load(file.path(projectpath,"Aligned","shift.Rdata"))
		dir.create(file.path(projectpath,"Edges","dat","Model samples"),showWarnings=FALSE,recursive=TRUE)
		dir.create(file.path(projectpath,"Edges","dat","Prediction samples"),showWarnings=FALSE,recursive=TRUE)
		dir.create(file.path(projectpath,"Edges","dat","Metadata"),showWarnings=FALSE,recursive=TRUE)
		cat("Checking integrity of existing window data files..\n")

		modelfilespath			<-  list.files(file.path(projectpath,"Edges","dat","Model samples"),full.names=TRUE)
	 	predictionfilespath <-  list.files(file.path(projectpath,"Edges","dat","Prediction samples"),full.names=TRUE)


		if(length(files) == length(c(model,prediction))){
			excludedfiles <-	numeric()
		
		}else{
			excludedfiles	<-	(1:length(files))[-c(model,prediction)]
		}

		## Check which windows that needs to be separated
		newlist 	<-  logical(max(windowlist))
		metadata  <-  list.files(file.path(projectpath,"Edges","dat","Metadata"),full.names=TRUE) # Alla f<U+00F6>nster med metadata
		metadata  <-  metadata[as.numeric(substr(basename(metadata),4,6)) %in% windowlist] # Alla f<U+00F6>nster med metadata som ska kollas
		newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(metadata),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns metadata till
		newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(modelfilespath),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns modellfiler till
		newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(predictionfilespath),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns prediktionsfiler till
		
		if(length(metadata)){
			n <-  0
	 		for(k in as.numeric(substr(basename(metadata),4,6))){
	 			n <-  n + 1
	 			load(metadata[n])
	 			if(!(all.equal(modelfiles,model) == TRUE & all.equal(predictionfiles,prediction) == TRUE))
	 			newlist[k]  <-  TRUE # Om metadata inte st<U+00E4>mmer <U+00F6>verens med nuvarande inst<U+00E4>llningar
			}
		}
		
		windowlist  <-  which(newlist)
		
		# Remove data files if they exists

		if(any(as.numeric(substr(basename(modelfilespath),4,6)) %in% windowlist))
 	    file.remove(modelfilespath[which(as.numeric(substr(basename(modelfilespath),4,6)) %in% windowlist)])

		if(any(as.numeric(substr(basename(predictionfilespath),4,6)) %in% windowlist))
 	  	file.remove(predictionfilespath[which(as.numeric(substr(basename(predictionfilespath),4,6)) %in% windowlist)])
		##

	 	cat("Windows to separate: [ ", windowlist," ]\n")
		cat("In total", length(windowlist),"windows\n")
	 	filenum 				<-  0
		modelfiles  		<-  model
		predictionfiles <-  prediction
	 	
	 	
	 	for(i in sort(c(model,prediction))){
	 		filenum <-  filenum +1
	 		cat("Loading ",basename(files[i]),paste("(",filenum,"/",length(c(model,prediction)),") ",sep=""))
			load(files[i])
			cat("... Exporting windows... \n")
		
			for(j in windowlist){
				if(nrow(Xbc)<(edges[j]+shift[i])){
					a <- 0
				}else if(nrow(Xbc)<(edges[j+1]+shift[i])){
					a <- (edges[j]+shift[i]):nrow(Xbc)
				}else
					a <- edges[j]:edges[j+1]+shift[i]  # Check so we are not exceeding matrix index
					
				Xbc <-  cbind(Xbc,matrix(0,nrow = nrow(Xbc),ncol = ifelse(maxMZ > ncol(Xbc),maxMZ-ncol(Xbc),0)))                                                           # Add columns if necessary so rows in outputfile will have same length
				x	<-	rbind(Xbc[a,],matrix(0,nrow=sum((0:(edges[j+1]+shift[i]-nrow(Xbc)))>0),ncol=ncol(Xbc)))
	   	  		write.table(t(x),file=file.path(projectpath,"Edges","dat",ifelse(i %in% model,"Model samples","Prediction samples"),paste("win",ifelse(j<=99,ifelse(j<=9,paste("00",j,sep=""),paste("0",j,sep="")),j),".dat",sep="")),append=TRUE,row.names=FALSE,col.names=FALSE)
			}
		}
		
		#All windows exported successfully; save metadata
		for(j in windowlist)
			save(modelfiles,predictionfiles,excludedfiles,file=file.path(projectpath,"Edges","dat","Metadata",paste("win",ifelse(j<=99,ifelse(j<=9,paste("00",j,sep=""),paste("0",j,sep="")),j),"_metadata.Rdata",sep="")))
	  
	}else
		cat("No windows selected!\n")
}

