##' Function new_samples
##' 
##' Function new_samples
##' @param projectpath
new_samples<-function(projectpath){
	
	## Skapa prediktionskatalog
	startup <-  1
	
	while(startup){
		cat("\n\n===============================\n")
		cat("===== Process new samples =====\n")
		cat("===============================\n\n")
		startup	<-	menu(c("Set path","Back"),title="")
		
		if(startup == 1){
			temp	<-	tk_choose.dir(caption="Select directory.")
			
			if(!is.na(temp)){
				predpath <-  temp
				if(length(list.files(projectpath))){
					cat("¤¤¤¤\n¤¤¤¤\n")
					cat("Warning! ",predpath," is not empty! Some files might be overwritten!\n")
					cat("¤¤¤¤\n¤¤¤¤\n")
				}
				
				a		<-	1
				startup	<-  0
								
				file.copy(file.path(projectpath,"SETTINGS.Rdata"),file.path(predpath))
				#file.copy(file.path(projectpath,"sampleinfo.Rdata"),file.path(predpath))
				file.copy(file.path(projectpath,"maxMZ.Rdata"),file.path(predpath))
				dir.create(file.path(predpath,"Edges"),showWarnings=FALSE)
				dir.create(file.path(predpath,"ALKANE_SERIE"),showWarnings=FALSE)
				dir.create(file.path(predpath,"HMCR"),showWarnings=FALSE)
				
				if(file.exists(file.path(projectpath,"HMCR","MCR.Rdata")))
					file.copy(file.path(projectpath,"HMCR","MCR.Rdata"),file.path(predpath,"HMCR"))
				
				if(file.exists(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata")))
			  		file.copy(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata"),file.path(predpath,"ALKANE_SERIE"))
				
				file.copy(file.path(projectpath,"Edges","edges.Rdata"),file.path(predpath,"Edges"))
			
				if(file.exists(file.path(predpath,"newsamples.Rdata")))
 						load(file.path(predpath,"newsamples.Rdata"))
 					
 				}else{
				cat("\n\n Error! You must specify where to store prediction files to continue!\n")
		  		a <-  0
			}
			
		}else if(startup == 2 | !startup){
			a		<-	0
			startup	<-	0
		}
	}
		
		
	### MENY
	while(a){
		cat("\n\n===============================\n")
		cat("===== Process new samples =====\n")
		cat("===============================\n\n")
		cat("Path:\t",predpath,"\n\n")

		predmenu  	<-  c("Import and smooth data","Select samples","Align samples","Check for outliers","Compress data","Resolve data","Export to NIST","Quit")

		if(exists("newsamples")){
			
			if(length(newsamples)){
				predvector  <-  c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE)
				cat("Number of prediction files selected: ",length(newsamples),"\n")
			}else
				predvector	<-	c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
		}else
			predvector	<-	c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
		a  <-  menu(predmenu[predvector],title="")
				
		if(!a | a == sum(predvector)){
			a <-  0
			
		}else if(a == 1){
			#Import data
			importmenu  <-  1
		
			while(importmenu){
				importmenu  <-  menu(c("Leco CSV","Leco CSV (special)","Andi NetCDF","Done"),title = "Import files")
					
				if(importmenu == 4 | !importmenu)
					importmenu  <-  0
					
				if(importmenu == 1)
					cat("\nSorry, not implemented yet\n")
					
				else if(importmenu == 2)
					cat("\nSorry, not implemented yet\n")
					
				else if(importmenu == 3)
					read_cdf(predpath)
			}
		}else if(a == 2){
			b <-  1
				
			while(b){
				b  <-  menu(c("Select samples", "Remove files","List files","Done"),title="Select samples")
					
				if(b == 1){
					temp		<-	tk_choose.files(caption="Select .Rdata files.", multi = TRUE,filters = matrix(c("Rdata files (*.Rdata", "*.Rdata"),1,2,byrow=T))
				  		
			  		if(length(temp)){
				  			
			  			if(exists("newsamples"))
			  				newsamples  <-  unique(c(newsamples,temp))
						else
					  		newsamples  <-  temp
						predvector  <-  c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
						more  <-  1
							
						while(more == 1){
							more  <-  menu(c("Yes","No"),title="Choose files from another directory?")
							if(more == 1){
								temp		<-	tk_choose.files(,caption="Select .Rdata files.", multi = TRUE,filters = matrix(c("Rdata files (*.Rdata", "*.Rdata"),1,2,byrow=T))
								newsamples  <-  unique(c(newsamples,temp))
				
						 	}
						}
					}
					
				}else if(b == 2){
					remfiles  	<-  select.list(basename(newsamples),multiple=TRUE,title="Remove samples from prediction set")
					newsamples  <-  newsamples[!(basename(newsamples) %in% remfiles)]
					
				}else if(b == 3){
					if(exists("newsamples"))
						select.list(basename(newsamples))
					else
						cat("No samples selected!\n")
				}else if(b == 4 | !b){
					b <-  0
					save(newsamples,file=file.path(predpath,"newsamples.Rdata"))
				}
			}
			
		}else if(a == 3){
				
			if(file.exists(file.path(projectpath,"Aligned","target.Rdata"))){
				load(file.path(projectpath,"Aligned","target.Rdata"))
				cat("Aligning samples..\n")
				aligndatapred(projectpath,predpath,newsamples,target,minshift,datasource,ion)
				
			}else
				cat("Error! No original target file found!")
				
		}else if(a == 4){
			check_data(predpath)
							
		}else if(a == 5){
			HDcompdata	<-	sigma_unfold(predpath)
			scoresAR(predpath,HDcompdata)
				
		}else if(a == 6){
			cat("Resolving data.. \n")
			results	<-	find_spectrum2(predpath,projectpath)
					
			if(length(results)){
				read_win2(predpath,results$type,results$windowlist)
				cat("Exporting spectrum to NIST...\n")
				spec2NIST2(predpath,results$type,all=TRUE)
			}
				
		}else if(a == 7){
			#Export to NIST
			spec2NIST2(predpath,results$type)
				
		}else if(a == 8 | !a)
			a <-  0
	}
}

