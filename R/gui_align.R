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
  dev.new()
  par(bg = "white")
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

