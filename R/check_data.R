##' Function check_data
##' 
##' Function check_data
##' @export
##' @param projectpath
check_data <- function(projectpath){
  
  
  ## Loading the list of files and setting variables
	load(file.path(projectpath,"Aligned","files.Rdata"))
  OBS			<-	1:length(files)
	s     		<-	character()
	COMP  		<-	1:5
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
				s[i]	<-	paste("[0", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
      		else
				s[i]	<-	paste("[00", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
		}
		else
			s[i]	<-	paste("[", i,"] ",sub("(.+)[.][^.]+$", "\\1", basename(files[i])),sep="")
  	}
  		
 		stemp <-  s

	outliers 	<-  character()
	k			<- 1
        X11.options(type="cairo")
        X11()

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

			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s)
                        replot<-TRUE
                                
		}else if(k == "Scaling"){
			scalenum 		<-  menu(c("Center","UV-scaled","Pareto-scaled","None"),title="Choose scaling method")
			scale 			<-  scaling(Xin,scalenum)
			X				<-	scale$Xout
			scale_text		<-  scale$scale_text
			vec				<-  pca(X,2)$vec
 			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s)
                        replot<-TRUE
  	
                }else if(k == "Export to excel"){
                        write.table(Xin,file=file.path(projectpath,"Check Data",paste(maintext,".txt",sep="")),sep="\t",row.names = stemp)
                        cat("Data exported to: ",file.path(projectpath,"Check Data",paste(maintext,".txt\n",sep=""),fsep="\\"))
	
		}else if(k == "List files"){
			select.list(s,title="Outliers")
		
		}else if(k == "Mark outliers"){
      
      ## the data points in mark are excluded
			mark  <-  select.list(s,preselect=s[1],multiple=TRUE,title="Mark outliers to remove.")
			s[s %in% mark] <-  paste(stemp[s %in% mark]," [OUT]")
			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s)
                        replot<-TRUE
                        
		}else if(k == "Unmark outliers"){
			unmark  <-  select.list(s,preselect=s[1],multiple = TRUE,title = "Unmark outliers.")
			s[s %in% unmark] <-  stemp[s %in% unmark]
			replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s)
                        replot<-TRUE
                }else if(k == "Zoom in"){
                        coords<-locator(n=2,type="p")
                        replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s,coords)
		}else if(k == "Standard zoom"){
                        replot_checkdata(vec[,1],vec[,2],scale_text,maintext,OBS,ss(X),s)
                }
	}
	graphics.off()
	remove_outliers(projectpath, grep("[OUT]",s,fixed=TRUE))
}

