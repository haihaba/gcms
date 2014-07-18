gui_align <-
function(projectpath,A_DATA){

	menuvector			<-	c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)
 	p					<-	c(0,0)
	A_DATA$targetfile	<-	"Median sample"
	A_DATA$ok			<-	0
	k					<-	1
        X11.options(type="cairo")
        X11()
        par(bg = "white")
	split.screen(c(2,1))
	while(k){
		
		menulist <-  c("Data source","List files","Choose target","View current settings","New IC","Align data","List fileshift","Zoom controls","Save and quit")
		k	<-	menu(menulist[menuvector],title="\nAlign data")
		if(k == sum(menuvector) | !k){
			k <- 0
		  	cat("\nSaving data...\n")
		  	save(A_DATA,file=file.path(projectpath,"Aligned","ADATA.Rdata"))   # Save and quit
		  
			}else if(menulist[menuvector][k] == "Data source"){
				datasource	<-	menu(choices=c("TIC","Basepeak","IC"),title="\nChoose data source")
				if(datasource){
					menuvector[c(3,5,6,7)]  <-  c(TRUE,TRUE,TRUE,FALSE)
					if(datasource == 1){
						DATA				<-	A_DATA$tic
						A_DATA$datasource 	<-	"TIC"
		  				maintext			<-	"Total ion current"
		  			}else if(datasource == 2){
		  				DATA				<-	A_DATA$bp
		  				A_DATA$datasource	<-	"Basepeak"
		  				maintext			<-	"Basepeak raw"
		  			}else if(datasource == 3){
		  				DATA      <-  A_DATA$ic
						A_DATA$datasource <-  "IC"
		  				maintext  <-  "Ion Cromatogram raw"
					}

					if(A_DATA$targetfile != "Median sample"){
						A_DATA$target	<-	DATA[fileindex,]
					}else
						A_DATA$target	<-	apply(DATA,2,median)
					
					screen(1)
					plot(1:length(A_DATA$target),A_DATA$target,main=c(maintext,paste("Target file: ",A_DATA$targetfile)),type="l",ylab=A_DATA$datasource,xlab="Time (scans)")
				}
				
			}else if(menulist[menuvector][k] == "List files"){
				cat("\n",paste(basename(A_DATA$files),"\n"))
			
			}else if(menulist[menuvector][k] == "Choose target"){
				menuvector[7]		<-	FALSE
	  			A_DATA$targetfile	<-	select.list(c(basename(A_DATA$files),"Median sample"),multiple = FALSE)
	  			
	  			if(!nchar(A_DATA$targetfile)){
	  				cat("\nNo target file selected! Using median sample...\n")
					A_DATA$targetfile	<-	"Median sample"
				
				}else if(A_DATA$targetfile != "Median sample"){
					fileindex			<-  which(basename(A_DATA$files) == A_DATA$targetfile)
    			
    			}else if(A_DATA$targetfile == "Median sample")
					A_DATA$targetfile	<-	"Median sample"

				if(A_DATA$targetfile != "Median sample"){
					A_DATA$target	<-	DATA[fileindex,]
				}else
					A_DATA$target	<-	apply(DATA,2,median)
				
				screen(1)
				
				plot(1:length(A_DATA$target),A_DATA$target,main=c(maintext,paste("Target file: ",A_DATA$targetfile)),type="l",ylab=A_DATA$datasource,xlab="Time (scans)")
			}else if(menulist[menuvector][k] == "View current settings"){
				cat("\n\nData source:\t",A_DATA$datasource,"\n")
				cat("Target file:\t",A_DATA$targetfile,"\n")
				cat("Ion:\t\t",A_DATA$ion,"\n\n")
			}else if(menulist[menuvector][k] == "New IC"){
				menuvector[7]  			<-  FALSE
			
			while(!is.numeric(A_DATA$ion	<-	round(as.numeric(readline("Input new IC (default = 298): ")))))
				cat("\n Ion must be numeric!\n")
			 for(i in 1:length(A_DATA$files)){
			 	cat("Loading ", paste(basename(A_DATA$files[i])," (",i,"/",length(A_DATA$files),")\n",sep=""))
				load(A_DATA$files[i])
    			ionindex	<-	which(SCAN_RANGE == A_DATA$ion)
		    	A_DATA$ic[i,1:nrow(Xbc)]		<-	t(Xbc[,ionindex])
			}
			
			rm(Xbc,SCAN_RANGE,SCAN_INFO)

		}else if(menulist[menuvector][k] == "Align data"){
			if(file.exists(file.path(projectpath,"SETTINGS.RDATA"))){
				load(file.path(projectpath,"SETTINGS.Rdata"))
		    	max_shift	<-	SETTINGS$MS
			}else{
				cat("No settings for maximum shift found! Using default value (500)")
			  	max_shift <-  500
			}
			
			TT	<-	matrix(0,max_shift*2+1,ncol(DATA)+  max_shift*2+1)
			for(i in 1:nrow(TT))
    			TT[i,i:(ncol(DATA)-1+i)]	<- A_DATA$target
    		
      		I				<-	apply(TT[,(max_shift+1):(ncol(TT)-max_shift-1)]%*%t(DATA),2,which.max)
			A_DATA$shift	<-	I-(max_shift+1)
			A_DATA$minshift	<-	min(A_DATA$shift)
			A_DATA$shift	<-	A_DATA$shift - A_DATA$minshift

			ADATA			<-	matrix(0,nrow(DATA),ncol(DATA))
    		start			<-	1+A_DATA$shift[1:nrow(DATA)]
     		stop			<-	ncol(DATA)-max(A_DATA$shift)+A_DATA$shift[1:nrow(DATA)]

			for(i in 1:nrow(DATA))
        		ADATA[i,1:(stop[i]-start[i]+1)]	<-	DATA[i,start[i]:stop[i]]
			
			replot_align(DATA,ADATA,start = start,stop = stop,zoomwidth = c(0,0),A_DATA$datasource,maintext,A_DATA$targetfile)
   			filesshift  <-  paste(basename(A_DATA$files)," ==> ",A_DATA$shift,sep="")
   			menuvector[c(7,8)]  			<-  TRUE
   			A_DATA$ok <-  1
		}else if(menulist[menuvector][k] == "List fileshift"){
			select.list(as.character(filesshift),multiple=FALSE,title="Shift with current settings")
		}else if(menulist[menuvector][k] == "Zoom controls"){
			zoomvector  <-	c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
	    	zoommenu  	<-	c("Zoom in","Zoom default","Zoom out","Pan left","Pan right","Back")
			zoom		<-	1
   			zoomwidth	<-	c(0,0)
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
					zoomwidth			<-	round(sort(locator(2,"p")$x))
          			replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
				}else if(zoom == 2){
					zoomwidth	<-	c(0,0)
					replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
				}else if(zoom == 3){
					if(sum(zoomwidth) == 0){
						cat("You must zoom in before you zoom out!\n")
					}else{
						zoomwidth[1]  <-  zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
						zoomwidth[2]  <-  zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
            			replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
					}

				}else if(zoom == 4){
					if(sum(zoomwidth) == 0){
						cat("You must zoom in before you can pan!\n")
					
					}else if(zoomwidth[1] <= 0){
						cat("Beginning of chromatogram.\n")
					
					}else{
						tmp1  		<-	zoomwidth[1] - (zoomwidth[2]-zoomwidth[1])/3
						tmp2  		<-	zoomwidth[2] - (zoomwidth[2]-zoomwidth[1])/3
 						zoomwidth	<-	c(tmp1,tmp2)
 						
 						replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
					}

				}else if(zoom == 5){
					
					if(sum(zoomwidth) == 0){
						cat("You must zoom in before you can pan!\n")
					
					}else if(zoomwidth[2] >= ncol(DATA)){
						cat("End of chromatogram.\n")
						
					}else{
						tmp1  		<-  zoomwidth[1] + (zoomwidth[2]-zoomwidth[1])/3
						tmp2  		<-  zoomwidth[2] + (zoomwidth[2]-zoomwidth[1])/3
						zoomwidth <-  c(tmp1,tmp2)
						replot_align(DATA,ADATA,start,stop,zoomwidth,A_DATA$datasource,maintext,A_DATA$targetfile)
					}
				}
			}
		}
	}
	close.screen(all = TRUE)
	A_DATA
}

