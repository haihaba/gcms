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
  
  require(tcltk)
  
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

##' Function gui_align
##' 
##' Function gui_align
##' @param projectpath
##' @param A_DATA
##' @return A_DATA
gui_align <- function(projectpath,A_DATA){
  
  ## inital menuvector to switch on/off menupoints of 'menulist' that is 
  ## defined further down
  menuvector			<-	c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)
  p					<-	c(0,0)
  
  
  ## adding new list entries to A_DATA: $targetfile
  A_DATA$targetfile	<-	"Median sample"
  
  ## set status of A_DATA to 0
  A_DATA$ok <- 0
  
  ## set flag for main loop below
  k <- 1
  
  ## open a new plot
  #dev.new()
  plot.new()
  #par(bg = "white")
  split.screen(c(2,1))
  
  ## main loop for 
  while(k){
    
    ## menulist[menuvector] is used to adjust the available menupoints
    ## to the current situation in the script
    menulist <-  c("Data source","List files","Choose target","View current settings",
                   "New IC","Align data","List fileshift","Zoom controls","Save and quit")
    
    ## initial menu settings: Data source, List Files, View current settings, Save and quit
    k	<-	menu(menulist[menuvector],title="\nAlign data")
    
    
    ## sum(vector is 4 and the last menupoint, Save and Quit), this 'if' is followed
    ## by a number of 'else if' for each menu point.
    if(k == sum(menuvector) | !k){
      
      ## k = 0 will exit the while loop above.
      k <- 0
      cat("\nSaving data...\n")
      save(A_DATA,file=file.path(projectpath,"Aligned","ADATA.Rdata"))
      
      
      ## Because the menu is variable, matching the choice has to be
      ## done by string matching here. On true, a new menu for selecting
      ## the data source comes up. 1=TIC, 2=Basepeak, 3=IC
    }else if(menulist[menuvector][k] == "Data source"){
      datasource <- menu(choices=c("TIC","Basepeak","IC"),title="\nChoose data source")
      
      ## when the data source is chosen, the main menu is adapted to:
      ## Data source, List files, Choose target, View current Settings, New IC, Align data, Save and quit.
      if(datasource){
        menuvector[c(3,5,6,7)] <- c(TRUE,TRUE,TRUE,FALSE)
        
        ## loading the data into DATA and adding the type
        ## of data into A_DATA$datasource
        if(datasource == 1){
          DATA <- A_DATA$tic
          A_DATA$datasource <- "TIC"
          maintext <- "Total ion current"
          
        }else if(datasource == 2){
          DATA <- A_DATA$bp
          A_DATA$datasource	<- "Basepeak"
          maintext <- "Basepeak raw"
          
        }else if(datasource == 3){
          DATA <- A_DATA$ic
          A_DATA$datasource <- "IC"
          maintext  <-  "Ion Chromatogram raw"  
        }
        
        
        ## on first call "Median sample" is default
        ## otherwise, targetfile should be a fileindex
        if(A_DATA$targetfile != "Median sample"){
          A_DATA$target	<-	DATA[fileindex,]
          ## calculate median sample  
        }else
          A_DATA$target	<-	apply(DATA,2,median)
        
        ## plot the data
        screen(1)
        plot(1:length(A_DATA$target),A_DATA$target,main=c(maintext,paste("Target file: ",A_DATA$targetfile)),type="l",ylab=A_DATA$datasource,xlab="Time (scans)")	
      }
      
      ## Menu selection: show all files included in the alignment	
    }else if(menulist[menuvector][k] == "List files?"){
      cat("\n",paste(basename(A_DATA$files),"\n"))
      
      
      ## Menu selection: choose target 
    }else if(menulist[menuvector][k] == "Choose target"){
      menuvector[7] <- FALSE
      
      ## choose either a file or median samples as target to align
      A_DATA$targetfile	<-	select.list(c(basename(A_DATA$files),"Median sample"),multiple = FALSE)
      
      ## check user entry
      if(!nchar(A_DATA$targetfile)){
        cat("\nNo target file selected! Using median sample...\n")
        A_DATA$targetfile	<-	"Median sample"
        
        ## find the file index of the chosen file
      }else if(A_DATA$targetfile != "Median sample"){
        fileindex			<-  which(basename(A_DATA$files) == A_DATA$targetfile)
        
        ## Don't understand this one. Isn't this overkill?
      }else if(A_DATA$targetfile == "Median sample")
        A_DATA$targetfile	<-	"Median sample"
      
      ## write the data for the chosen sample inteo A_DATA$target
      if(A_DATA$targetfile != "Median sample"){
        A_DATA$target	<-	DATA[fileindex,]
        
        ## calculate and write the median into A_DATA$target
      }else
        A_DATA$target	<-	apply(DATA,2,median)
      
      
      ## Select upper screen and plot 'target' data 
      screen(1)	
      plot(1:length(A_DATA$target),A_DATA$target,main=c(maintext,paste("Target file: ",A_DATA$targetfile)),
           type="l",ylab=A_DATA$datasource,xlab="Time (scans)")
      
      
      
    }else if(menulist[menuvector][k] == "View current settings"){
      cat("\n\nData source:\t",A_DATA$datasource,"\n")
      cat("Target file:\t",A_DATA$targetfile,"\n")
      cat("Ion:\t\t",A_DATA$ion,"\n\n")
      
      ## choose new ion channel 
    }else if(menulist[menuvector][k] == "New IC"){
      menuvector[7] <- FALSE
      
      ## user input, check for numeric
      while(!is.numeric(A_DATA$ion	<-	round(as.numeric(readline("Input new IC (default = 298): ")))))
        cat("\n Ion must be numeric!\n")
      
      ## load single ion channel data
      for(i in 1:length(A_DATA$files)){
        cat("Loading ", paste(basename(A_DATA$files[i])," (",i,"/",length(A_DATA$files),")\n",sep=""))
        load(A_DATA$files[i])
        ionindex <- which(SCAN_RANGE == A_DATA$ion)
        A_DATA$ic[i,1:nrow(Xbc)]		<-	t(Xbc[,ionindex])  
      }
      
      rm(Xbc,SCAN_RANGE,SCAN_INFO)
      
      
      ## 'else if' sub loop for alignment
    }else if(menulist[menuvector][k] == "Align data"){
      
      ## loading settings for max shift, else set default 500
      if(file.exists(file.path(projectpath,"SETTINGS.RDATA"))){
        load(file.path(projectpath,"SETTINGS.Rdata"))
        max_shift	<-	SETTINGS$MS
      }else{
        cat("No settings for maximum shift found! Using default value (500)")
        max_shift <-  500  
      }
      
      ## prepare matrix for the shifted target
      TT	<-	matrix(0,max_shift*2+1,ncol(DATA)+  max_shift*2+1)
      
      ## put target into TT matrix with shifts from
      ## -maxShift until + maxShift applied
      for(i in 1:nrow(TT))
        TT[i,i:(ncol(DATA)-1+i)]	<- A_DATA$target
      
      ## matrix multiplication of TT with the original data
      ## then apply to find indices of max values. Max value is where
      ## there is the best overlap between target and chromatogram to shift
      I <- apply(TT[,(max_shift+1):(ncol(TT)-max_shift-1)]%*%t(DATA),2,which.max)
      
      ## store shift values
      A_DATA$shift <- I-(max_shift+1)
      A_DATA$minshift <- min(A_DATA$shift)
      A_DATA$shift <- A_DATA$shift - A_DATA$minshift
      
      ## calculate start/stop values to shift the chromatograms
      ADATA <- matrix(0,nrow(DATA),ncol(DATA))
      start <- 1+A_DATA$shift[1:nrow(DATA)]
      stop <- ncol(DATA)-max(A_DATA$shift)+A_DATA$shift[1:nrow(DATA)]
      
      ## application of shift values to the actual data
      for(i in 1:nrow(DATA))
        ADATA[i,1:(stop[i]-start[i]+1)]	<-	DATA[i,start[i]:stop[i]]
      
      ## plot the alignment
      replot_align(DATA,ADATA,start = start,stop = stop,zoomwidth = c(0,0),A_DATA$datasource,maintext,A_DATA$targetfile)
      
      ## store fileshifts for 'list fileshift' menu point 
      filesshift  <-  paste(basename(A_DATA$files)," ==> ",A_DATA$shift,sep="")
      
      ## activate 'list fileshift' and 'zoom control'
      menuvector[c(7,8)] <- TRUE
      
      ## data is ready for further processing
      A_DATA$ok <- 1
      
      ## list fileshifts
    }else if(menulist[menuvector][k] == "List fileshift"){
      select.list(as.character(filesshift),multiple=FALSE,title="Shift with current settings")
      
      ## Zoom control for the aligned chromatograms
    }else if(menulist[menuvector][k] == "Zoom controls"){
      zoomvector  <- c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
      zoommenu    <- c("Zoom in","Zoom default","Zoom out","Pan left","Pan right","Back")
      zoom        <- 1
      zoomwidth   <- c(0,0)
      
      while(zoom){
        if(sum(zoomwidth) != 0)
          cat("\n\nCurrent zoom width:\t",zoomwidth[2] - zoomwidth[1]," scans.\n")
        
        cat("============================\n")
        
        zoom  <-  menu(choices=zoommenu[zoomvector],title = "Zoom Controls")
        
        if(zoom == sum(zoomvector) | !zoom){
          zoom  		<-	0
          zoomwidth	<-	c(0,0)
          
          # Zoom back to default after quit
          replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
          
        }else if(zoom == 1){
          cat("\n\nSelect two points on the x-axis to zoom in.\n\n")
          zoomwidth <- round(sort(locator(2,"p")$x))
          replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
        }else if(zoom == 2){
          zoomwidth	<-	c(0,0)
          replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
        }else if(zoom == 3){
          if(sum(zoomwidth) == 0){
            cat("You must zoom in before you zoom out!\n")
          }else{
            zoomwidth[1] <- zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
            zoomwidth[2] <- zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
            replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
          }
          
        }else if(zoom == 4){
          if(sum(zoomwidth) == 0){
            cat("You must zoom in before you can pan!\n")
            
          }else if(zoomwidth[1] <= 0){
            cat("Beginning of chromatogram.\n")
            
          }else{
            tmp1      <- zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
            tmp2      <- zoomwidth[2] - (zoomwidth[2]-zoomwidth[1])/3
            zoomwidth <- c(tmp1,tmp2)
            
            replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
          }
          
        }else if(zoom == 5){
          
          if(sum(zoomwidth) == 0){
            cat("You must zoom in before you can pan!\n")
            
          }else if(zoomwidth[2] >= ncol(DATA)){
            cat("End of chromatogram.\n")
            
          }else{
            tmp1      <-  zoomwidth[1] + (zoomwidth[2]-zoomwidth[1])/3
            tmp2  		<-  zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
            zoomwidth <-  c(tmp1,tmp2)
            replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
          }
        }
      }
    }
  }
  close.screen(all = TRUE)
  return(A_DATA)
}

##' Function replot_align
##' 
##' Function replot_align
##' @param DATA  the target data
##' @param ADATA the shifted data
##' @param start start position for the chromatograms
##' @param stop stop position for the chromatograms
##' @param zoomwidth zoomwidth 
##' @param datasource
##' @param maintext
##' @param targetfile
replot_align<-function(DATA,ADATA,start,stop,zoomwidth = c(0,0),datasource,maintext,targetfile){
  
  zoomwidth[1]	<-	max(0,zoomwidth[1])
  zoomwidth[2]	<-	min(zoomwidth[2],ncol(DATA))
  if(sum(zoomwidth) == 0){
    ## Plotting without zoom
    screen(1)
    for(i in 1:nrow(DATA)){
      plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim=c(0,ncol(DATA)*1.01),ylim=c(0,max(DATA[,min(start):max(stop)])*1.01),main=c(maintext,paste("Target file: ", targetfile)),xlab="",ylab=datasource)
      par(new=TRUE)
    }
    
    if(!missing(ADATA)){
      
      screen(2)
      for(i in 1:nrow(ADATA)){
        plot(1:length(ADATA[i,]),ADATA[i,],type="l",col=i,xlim=c(0,ncol(ADATA)*1.01),ylim=c(0,max(DATA[,min(start):max(stop)])*1.01),main="Aligned data",xlab="Time (scans)",ylab=datasource)
        par(new=TRUE)
      }
    }
    
    ## Plot with zoom  
  }else{
    
    ## Plot unaligned data in upper half of graph:
    xlimit  <-  c(max(0,zoomwidth[1]),min(zoomwidth[2],ncol(DATA)))
    #ylimit  <-  c(0,max(ADATA[,zoomwidth[1]:zoomwidth[2]],DATA[,min(start):max(stop)])*1.01)  ### did not adjust the ylimit
    ylimit  <-  c(0,max(ADATA[,zoomwidth[1]:zoomwidth[2]])*1.01)
    screen(1)
    
    for(i in 1:nrow(DATA)){
      plot(1:length(DATA[i,]),DATA[i,],type="l",col=i,xlim = xlimit,ylim = ylimit,main=c(maintext,paste("Target file: ", targetfile)),xlab="",ylab=datasource)
      par(new=TRUE)
    }
    
    if(!missing(ADATA)){
      ## Plot aligned data in lower half of graph:
      screen(2)
      for(i in 1:nrow(ADATA)){
        plot(1:length(ADATA[i,]),ADATA[i,],type="l",col=i,xlim = xlimit,ylim = ylimit,main="Aligned data",xlab="Time (scans)",ylab=datasource)
        par(new=TRUE)
      }
    }
  }
}


