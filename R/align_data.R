##' Chromatogram alignment by maximizing the covariance between the
##' ion-channels or total ion counts (TIC) of the individual chromatograms.
##'
##' creates a number of files:
##' \itemize{
##' \item{A_DATA}{list with 'TIC','ic','bp'}
##' \item{ADATA}{list with 'TIC','ic','bp'}
##' }
##' @export
##' @title Alignment of GC/MS chromatograms in \code{GCMS} package
##' @param projectpath base path from \code{GCMS} processing script
##' @return a directory \code{alignement} is created that includes the
##' related to chromatogram alignment. 
##' @author Lorenz Gerber
align_data <-function(projectpath){
  
  
  
  ## Initial Menu for Alignment
  k <- menu(c("Load files to align","Load previous settings","Quit"),title="\nAlign data")
  
  
  
  ## A_DATA is a list that will contain TIC, Basepeak and IC chromatogram
  A_DATA		<-	list()
	A_DATA$ok	<-	0
 	
  
  ## Load files to align, k == 1 
  if(k == 1){
     
     ## ask for Ion Chromatogram to extract (IC), check for numeric 
     while(!is.numeric(mz	<-	round(as.numeric(readline("Input IC (default = 298): ")))))
       cat("\n Ion must be numeric!\n")
     
     ## if no value entered, use default 298
     if(is.na(mz))
       mz  <-  298
     
     ## select files to be aligned and sort, usees TclTk
     alignfiles <- tk_choose.files(caption="Select .Rdata files to align.", multi = TRUE,filters = matrix(c("Rdata files (*.Rdata", "*.Rdata"),1,2,byrow=TRUE))
     alignfiles <- alignfiles[order(alignfiles)]
     
     
     ## select more files for processing
     more <- 1
     while(more){
       od  <-  menu(c("Yes","No"),title="Do you want to select more files from another directory?")
       
       if(od == 1)
         ## alignfiles contains the files with path
         alignfiles  <-  unique(c(alignfiles,tk_choose.files(caption="Select .Rdata file to align.", multi = TRUE,filters = matrix(c("Rdata files (*.Rdata", "*.Rdata"),1,2,byrow=TRUE))))
       else
         more <- 0
     }
		
		if(!length(alignfiles)){
			cat("No files selected!\n")
		
		}else{
      
      ## how many files
      n	<- length(alignfiles)
      
      ## create directory for aligned files
			dir.create(file.path(projectpath,"Aligned"),showWarnings=FALSE)
 			
      
			## construct filenames without extensions fort COLUMNID1
      COLUMNID1 		<-	cbind(sub("[.][^.]*$", "", basename(alignfiles)))
			A_DATA$files	<-	files	<-	alignfiles
			A_DATA$fsize	<-	matrix(0,n,2)
      
      ## Loop through the number of files to align
      for(i in 1:n){
        cat("Loading ", paste(basename(alignfiles[i])," (",i,"/",n,")\n",sep=""))
        
        ###load CDF imported RData files
        load(alignfiles[i])
        A_DATA$fsize[i,]		<-	dim(Xbc)
        num									<-	which(SCAN_RANGE == mz)       #Om length(num) == 0?
			
				if(i == 1){
					NUM_MZ			<-	ncol(Xbc)
 					NUM_scans		<-	nrow(Xbc)
 					A_DATA$tic	<-	matrix(0,n,NUM_scans)
 					A_DATA$bp		<-	matrix(0,n,NUM_scans)
 					A_DATA$ic		<-	matrix(0,n,NUM_scans)
 					SUM_MZ      <-  matrix(0,n,NUM_MZ)
 				}
   			
   				if(nrow(Xbc) > NUM_scans){
   					A_DATA$tic	<-	cbind(A_DATA$tic,matrix(0,n,nrow(Xbc)-NUM_scans))
  					A_DATA$bp		<-	cbind(A_DATA$bp,matrix(0,n,nrow(Xbc)-NUM_scans))
  					A_DATA$ic		<-  cbind(A_DATA$ic,matrix(0,n,nrow(Xbc)-NUM_scans))
  					NUM_scans 	<-  nrow(Xbc)
   				}
   			
   				if(ncol(Xbc) > NUM_MZ){
   					SUM_MZ  <-  cbind(SUM_MZ,matrix(0,n,ncol(Xbc)-NUM_MZ))
					NUM_MZ  <-  ncol(Xbc)
				}
   			
        
   				A_DATA$tic[i,1:nrow(Xbc)]		<-	rowSums(Xbc)
   				A_DATA$bp[i,1:nrow(Xbc)]		<-	apply(Xbc,1,max)
   				A_DATA$ic[i,1:nrow(Xbc)]		<-	t(Xbc[,num])
   				SUM_MZ[i,1:ncol(Xbc)]				<-	colSums(Xbc)
			}
			
			A_DATA$ion	<-	mz
			A_DATA$ok		<-  1
			save(A_DATA,SUM_MZ,COLUMNID1,SCAN_RANGE,NUM_scans,file=file.path(projectpath,"Aligned","A_DATA.Rdata"))
		}
	
	}else if(k == 2){
		
		if(file.exists(file.path(projectpath,"Aligned","A_DATA.Rdata"))){
			load(file.path(projectpath,"Aligned","A_DATA.Rdata"))
			
			if(!length(A_DATA$files) | !is.character(A_DATA$files) | !A_DATA$ok){
				cat("Invalid data!\n")
				
			}else
			files <-  A_DATA$files
		
		}else{
			cat("Settings and alignment data not found in this project folder, please specify location.\n")
			aligndata	<-	tk_choose.files(default=file.path(projectpath,"Aligned","A_DATA.Rdata"),caption="Select align data settings file (A_DATA.Rdata)",multi = FALSE)
			
			if(!length(grep(pattern="^.+\\.([rR][dD][aA][tT][aA])$",aligndata))){
				cat("File is not .Rdata!\n\n")
			
			}else{
				load(aligndata)
				
				if(!length(A_DATA$files) | !is.character(A_DATA$files) | !A_DATA$ok){
					cat("Invalid data!\n")
				
				}else
					files <-  A_DATA$files
			}
		}
	}
	
	if(A_DATA$ok & k != 3 & k != 0){
		##  GUI ALIGN
		A_DATA  <-  gui_align(projectpath,A_DATA)
		##
		runagain  <-  1
		
		while(runagain){
			
			if(A_DATA$ok){
				shift			<-	A_DATA$shift
				start 		<-  1+shift[1:nrow(A_DATA$tic)]
				stop			<-	ncol(A_DATA$tic)-max(shift)+shift[1:nrow(A_DATA$tic)]
				BASEPEAK	<-	matrix(0,nrow(A_DATA$bp),ncol(A_DATA$bp))
				TIC				<-	matrix(0,nrow(A_DATA$tic),ncol(A_DATA$tic))
				for(i in 1:nrow(A_DATA$bp)){
					TIC[i,1:(stop[i]-start[i]+1)]			<-	A_DATA$tic[i,start[i]:stop[i]]
		  			BASEPEAK[i,1:(stop[i]-start[i]+1)]	<-	A_DATA$bp[i,start[i]:stop[i]]
				}
				
				TIC_RAW	<-	A_DATA$tic

				cat("Saving variables..\n")
				save(TIC_RAW,file=file.path(projectpath,"Aligned","TIC_RAW.Rdata"))
				save(shift,file=file.path(projectpath,"Aligned","shift.Rdata"))
				save(TIC,file=file.path(projectpath,"Aligned","TIC.Rdata"))
				save(BASEPEAK,file=file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
				save(SUM_MZ,file=file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
				save(COLUMNID1,file=file.path(projectpath,"Aligned","COLUMNID1.Rdata"))
				save(files,file=file.path(projectpath,"Aligned","files.Rdata"))
				save(SCAN_RANGE,file=file.path(projectpath,"Aligned","SCAN_RANGE.Rdata"))
				save(NUM_scans,file=file.path(projectpath,"Aligned","NUM_scans.Rdata"))

				textstr <-  paste(files," Shift: ", shift," Scans: ", A_DATA$fsize[,1]," MZ: ",A_DATA$fsize[,2],sep="\n")
				#do_log(projectpath,textstr)

				target 		<-	A_DATA$target
				minshift	<-	A_DATA$minshift
				ion			<-	A_DATA$ion
				datasource	<-	A_DATA$datasource
				targetfile	<-	A_DATA$targetfile

				save(target, minshift, datasource, ion,targetfile,file=file.path(projectpath,"Aligned","target.Rdata"))
				#do_log(projectpath,paste("Target file: ", targetfile))
				#do_log(projectpath,switch(A_DATA$datasource, TIC = "Data used: TIC", Basepeak  = "Data used: BASEPEAK", IC  = paste("Data used: Ion Chromatogram, m/z: ", ion)))
				#do_log(projectpath,"----------------------")
				runagain  <-  0
			
			}else{
				runagain  <-  menu(c("Yes","No"),title="No alignment was done! Are you sure you want to quit?")
				
				if(runagain == 1)
					runagain  <-  0
				else
					A_DATA  <-  gui_align(projectpath,A_DATA)
		 	}
		}
	}
graphics.off()
}

