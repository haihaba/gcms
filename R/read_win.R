##' Function read_win
##' 
##' Function read_win
##' @export
##' @param projectpath
##' @param type
##' @param win
##' @return DATA, SPECTRUM, VARID1, OBSID1, VARID2
read_win <- function(projectpath,type,win){
	
	SPECTRUM	<-	DATA	<-	VARID1	<-	VARID2	<-	numeric()
	load(file.path(projectpath,"Edges","edges.Rdata"))
	num_windows	<-	length(edges)-1
	load(file.path(projectpath,"Aligned","files.Rdata"))
	load(file.path(projectpath,"HMCR","MCR.Rdata"))
	
	if(!missing(win)){
		
		# if(length(list.files(file.path(projectpath,"HMCR",type,"win"))) > length(win)){
			# cat("\nRecently processed windows: ",win,"\n")
			# cat("Previously processed windows: ",	setdiff(as.numeric(substr(basename(list.files(file.path(projectpath,"HMCR",type,"win"))),4,6)),win),"\n")
	  	
	  		# if(1 == menu(c("Yes","No"),title="There are some processed windows that are not selected, do you wish to include some of them as well?"))
	    		# while(!length(win <-  sort(as.numeric(substr(basename(tk_choose.files(caption="Select windows to process")),4,6)))))
	      			# cat("You must select at least one window!\n")
		# }
	
	}else{
		cat("Select windows.")
   		while(!length(win <-  sort(as.numeric(substr(basename(tk_choose.files(caption="Select windows to process")),4,6)))))
			cat("You must select at least one window!\n")
	}
	
	cat("\n")
	
	for(i in win){
				if(file.exists(file.path(projectpath,"HMCR",type,"win",paste("win",ifelse(i<=99,ifelse(i<=9,paste("00",i,sep=""),paste("0",i,sep="")),i),".Rdata",sep="")))){
					load(file.path(projectpath,"HMCR",type,"win",paste("win",ifelse(i<=99,ifelse(i<=9,paste("00",i,sep=""),paste("0",i,sep="")),i),".Rdata",sep="")))
   			if(length(S)){
   				cat("Number of processed files: ",nrow(as.matrix(C)),"\n")
   				excludedfiles <-  eval(parse(text=paste("MCR$",tolower(type),"$excl",sep="")))
      				
   				if(nrow(as.matrix(C)) == (length(files) - length(excludedfiles))){
      				cat("--------------------------------------------------------\n")
      				SPECTRUM	<-	cbind(SPECTRUM,S)
      				DATA		<-	cbind(DATA,C)
      				colnames(DATA)<-NULL  
      				cat("Window num: ",i,"\nNumber of comps: ",ncol(S),"\nTotal number of comps: ", ncol(SPECTRUM),"\n")
      				VARID2	<-	c(VARID2,TIME)
      			
    				for(j in 1:ncol(S))
      					VARID1	<-	rbind(VARID1,ifelse(i>99,ifelse(j<10,paste("Win",i,"_C0",j,sep=""),paste("Win ",i,"_C",j,sep="")),ifelse(i<10,ifelse(j<10,paste("Win00",i,"_C0",j,sep=""),paste("Win00",i,"_C",j,sep="")),ifelse(j<10,paste("Win0",i,"_C0",j,sep=""),paste("Win0",i,"_C",j,sep="")))))
    
    			}else
		 			
	 			cat("Number of files and rows in data matrix does not match!\n")
			}
		}
  	}
		
	load(file.path(projectpath,"Aligned","COLUMNID1.Rdata"))
		
	if(length(excludedfiles))
		OBSID1		<-	COLUMNID1[-excludedfiles]
	else
		OBSID1		<-	COLUMNID1
	save(DATA,VARID1,VARID2,OBSID1,SPECTRUM,file=file.path(projectpath,"HMCR",type,"MVA_DATA.Rdata"))
	save(SPECTRUM = SPECTRUM,file=file.path(projectpath,"HMCR",type,"SPECTRUM.Rdata"))
	save(VARID1	= VARID1,file=file.path(projectpath,"HMCR",type,"VARID1.Rdata"))
	save(VARID2 = VARID2,file=file.path(projectpath,"HMCR",type,"VARID2.Rdata"))
	write.table(data.frame(ID = OBSID1,DATA),file=file.path(projectpath,"HMCR",type,"MVA_DATA.txt"),sep='\t',row.names = FALSE,col.names = c("Primary ID",VARID1),quote=FALSE)
	out <-  list(DATA = DATA,SPECTRUM = SPECTRUM,VARID1 = VARID1, OBSID1 = OBSID1,VARID2 = VARID2)
}

