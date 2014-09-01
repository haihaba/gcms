##' Function gui_MCR
##' 
##' Function gui_MCR
##' @export
##' @param projectpath
gui_MCR<-function(projectpath){
	
	if(file.exists(file.path(projectpath,"Aligned","files.Rdata")) | file.exists(file.path(projectpath,"Aligned","shift.Rdata"))){
		load(file.path(projectpath,"Aligned","files.Rdata"))
		load(file.path(projectpath,"Aligned","shift.Rdata"))

	}else{
		cat("Error! You must align data first!\n")
		Sys.sleep(1)
		return(NULL)
	}

	load(files[which.min(shift)])
	EDGES_TIME	<-	SCAN_INFO[,2]
	save(EDGES_TIME,file=file.path(projectpath,"EDGES_TIME.Rdata"))
	out	<-	2
	temp	<-	5
  	while(out == 2 | temp == 5){
  		
  		if(file.exists(file.path(projectpath,"HMCR","MCR.Rdata"))){
  			load(file.path(projectpath,"HMCR","MCR.Rdata"))
  			cat("Current settings:\n\n")
			cat("Number of samples: ",length(files),"\n")
			cat("Windows to be processed: ",grep("[P]",MCR$windowlist),"\n\n")
			cat("Regular method settings:\n")
			cat("Included files: ",MCR$reg$incl,"\n")
			cat("Prediction files: ",MCR$reg$pred,"\n")
			cat("Excluded files: ",MCR$reg$excl,"\n\n")

			
			cat("Similarity criterion:\n")
			cat("Chrom. profile: ",MCR$cp,"\n")
			cat("Spec. profile: ",MCR$sp,"\n\n")
			
			out	<-	1 #menu(choices=c("Process","Change settings","Quit"),title="\nMCR")
			
			if(out == 3)
				temp  <-  6
			if(out == 1)
				temp	<-	1  #menu(c("Regular","Back"),title="\nMCR processing")
		}
		
		if(out == 2){
			MCR	<-	MCR_settings(projectpath)
			if(is.null(MCR))
			  return(NULL)
		}
	}
	
	if(temp == 1){ #Regular
		cat("======== Regular method ========\n\n")
		dir.create(file.path(projectpath,"HMCR","REG"),recursive = TRUE,showWarnings = FALSE)
		
		
		samples				<-	MCR$reg$samples
		incl				<-	MCR$reg$incl
		excl				<-	MCR$reg$excl
		pred				<-	MCR$reg$pred
		MCR$reg$excluded	<-	excl
  	
		
	 	save(MCR=MCR,file=file.path(projectpath,"HMCR","MCR.Rdata"))
	 	
	 	
	 	win				<-	find_spectrum(projectpath,incl,pred,grep("[P]",MCR$windowlist))
	 	window_data		<-	read_win(projectpath,"REG",win)
	 	warnings()
	 	
	 
	
	}
	
	
	if(temp == 1 ){ 
		#cat("Exporting REG-data NIST text file.\n")
		#spec2NIST(projectpath,type = "REG")
	}
		
	cat("Ending Method 2.\n\n")
}

