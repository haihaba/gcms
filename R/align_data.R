##' Chromatogram alignment by maximizing the covariance between the
##' ion-channels or total ion counts (TIC) of the individual chromatograms.
##'
##' creates a number of files:
##' \itemize{
##' \item{A_DATA}{list with 'TIC','ic','bp', 'ok', 'ion'}
##' \item{SUM_MZ}{sum of ion counts per mz}
##' \item{COLUMNID1}{filenames without extension}
##' \item{SCAN_RANGE}{sequence of all mz values, is obtained from rawimport files}
##' \item{NUM_scans}{nrow Xbc, how many data points}
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
  
  
  
  ## A_DATA is a list to collect most of the data related to the alignment.
  ## $ok, $files, $fsize, $tic, $bp, $ic, 
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
     } else {
       ## how many files
       n <- length(alignfiles)
       
       ## create directory for aligned files
       dir.create(file.path(projectpath,"Aligned"),showWarnings=FALSE)
       
       ## construct filenames without extensions for COLUMNID1
       COLUMNID1    <-  cbind(sub("[.][^.]*$", "", basename(alignfiles)))
       
       ## write filenames into alignfiles
       A_DATA$files <-  files <- alignfiles
       
       ## construct the matrix for $fsize, n rows for number of files
       ## and two columns, fsize will hold the dimensions of Xbc
       A_DATA$fsize <-  matrix(0,n,2)
       
       
       
       
       
       ## Loop through the number of files to align
       for(i in 1:n){
         cat("Loading ", paste(basename(alignfiles[i])," (",i,"/",n,")\n",sep=""))
         
         ## load CDF imported RData files
         load(alignfiles[i])
         
         ## fsize is the dimension of the main data table Xbc
         A_DATA$fsize[i,] <- dim(Xbc)
         
         ## num is the col id to choose the ion count channel
         num <- which(SCAN_RANGE == mz)
         
         ## On the first file to align, setting up variables
         if(i == 1){
           NUM_MZ     <-  ncol(Xbc)
           NUM_scans  <-  nrow(Xbc)
           A_DATA$tic <-  matrix(0,n,NUM_scans)
           A_DATA$bp  <-  matrix(0,n,NUM_scans)
           A_DATA$ic  <-  matrix(0,n,NUM_scans)
           SUM_MZ     <-  matrix(0,n,NUM_MZ)
         }
         
         ## Check in the following files (i>1) whether nrow(Xbc)>NUM_scans
         ## if so, extend the $tic, etc. matrices
         if(nrow(Xbc) > NUM_scans){
           A_DATA$tic <- cbind(A_DATA$tic,matrix(0,n,nrow(Xbc)-NUM_scans))
           A_DATA$bp <- cbind(A_DATA$bp,matrix(0,n,nrow(Xbc)-NUM_scans))
           A_DATA$ic <- cbind(A_DATA$ic,matrix(0,n,nrow(Xbc)-NUM_scans))
           NUM_scans <- nrow(Xbc)  
         }
         
         ## check in the files i>1 whether ncol(Xbc) > NUM_MZ
         ## in this case, correct SUM_MZ and NUM_MZ
         if(ncol(Xbc) > NUM_MZ){
           SUM_MZ  <-  cbind(SUM_MZ,matrix(0,n,ncol(Xbc)-NUM_MZ))
           NUM_MZ  <-  ncol(Xbc)  
         }
       
         
         ## calculate and write tic data
         A_DATA$tic[i,1:nrow(Xbc)] <- rowSums(Xbc)
         ## calculate and write basepeak data
         A_DATA$bp[i,1:nrow(Xbc)] <- apply(Xbc,1,max)
         ## calculate and write ion channel data
         A_DATA$ic[i,1:nrow(Xbc)] <- t(Xbc[,num])
         ## calculate and write sum per mz
         SUM_MZ[i,1:ncol(Xbc)] <- colSums(Xbc)
       }
       
       
       
       ## store the ion channel that has been
       ## used into $ion 
       A_DATA$ion <- mz
       
       
       ## flag to activate the gui_align menu
       A_DATA$ok <- 1
       
       
       ## store all data used for alignment
       save(A_DATA,SUM_MZ,COLUMNID1,SCAN_RANGE,NUM_scans,file=file.path(projectpath,"Aligned","A_DATA.Rdata"))
     }
     
     ## loading old settings (tic, basepeak, ic)
     }else if(k == 2){
       
       
       
       ## check if the file A_Data.Rdata exists, if yes, load it
       if(file.exists(file.path(projectpath,"Aligned","A_DATA.Rdata"))){
         
         load(file.path(projectpath,"Aligned","A_DATA.Rdata"))
         
         
         
         ## check for $files (length, character) and $ok, if ok, load $files   
         if(!length(A_DATA$files) | !is.character(A_DATA$files) | !A_DATA$ok){
           cat("Invalid data!\n")
         }else
           files <-  A_DATA$files 
       }else{
         
         
         ## when A_Data.Rdata doesn't exist, ask for the location where to look
         cat("Settings and alignment data not found in this project folder, please specify location.\n")
         aligndata <- tk_choose.files(default=file.path(projectpath,"Aligned","A_DATA.Rdata"),caption="Select align data settings file (A_DATA.Rdata)",multi = FALSE)
         
         
         
         ## check the type of file chose. has to be Rdata
         if(!length(grep(pattern="^.+\\.([rR][dD][aA][tT][aA])$",aligndata))){
           cat("File is not .Rdata!\n\n")
         }else{
           ## load file chosen by the user
           load(aligndata)
         
           
           
           ## check for $files (length, character) and $ok, if ok, load $files
           if(!length(A_DATA$files) | !is.character(A_DATA$files) | !A_DATA$ok){
             cat("Invalid data!\n")
           }else
             files <-  A_DATA$files
         }  
       }
     }
  
  ## Check condition to load gui_align
  if(A_DATA$ok & k != 3 & k != 0){
    ##  call gui_alig
		A_DATA  <-  gui_align(projectpath,A_DATA)
    
    ## flag for exiting on value 0
		runagain  <-  1
		
    ## loop around 
    while(runagain){
			
      ## if the data is ready for further processing
      ## A_DATA$ok should be set to 1 in gui_align
			if(A_DATA$ok){
				
        ## preparing for TIC and Basepeak correction
        shift			<-	A_DATA$shift
				start 		<-  1+shift[1:nrow(A_DATA$tic)]
				stop			<-	ncol(A_DATA$tic)-max(shift)+shift[1:nrow(A_DATA$tic)]
				BASEPEAK	<-	matrix(0,nrow(A_DATA$bp),ncol(A_DATA$bp))
				TIC				<-	matrix(0,nrow(A_DATA$tic),ncol(A_DATA$tic))
				
        ## Correcting TIC and Basepeak for the shift
        for(i in 1:nrow(A_DATA$bp)){
					TIC[i,1:(stop[i]-start[i]+1)]			<-	A_DATA$tic[i,start[i]:stop[i]]
		  			BASEPEAK[i,1:(stop[i]-start[i]+1)]	<-	A_DATA$bp[i,start[i]:stop[i]]
				}
				
				TIC_RAW	<-	A_DATA$tic
        
        ## writing all the data to file
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
				
				target 		<-	A_DATA$target
				minshift	<-	A_DATA$minshift
				ion			<-	A_DATA$ion
				datasource	<-	A_DATA$datasource
				targetfile	<-	A_DATA$targetfile

				save(target, minshift, datasource, ion,targetfile,file=file.path(projectpath,"Aligned","target.Rdata"))
				
        ## runagain set to 0, processing was succesful
				runagain  <-  0
			
			}else{
				## No alignment done yet, should it be tried again?
        runagain  <-  menu(c("Yes","No"),title="No alignment was done! Are you sure you want to quit?")
				if(runagain == 1)
					runagain  <-  0
				else
					A_DATA  <-  gui_align(projectpath,A_DATA)
		 	}
		}
	}
#graphics.off()
}

