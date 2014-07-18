spec2NIST <-
function(projectpath,type,all=FALSE){
	
	load(file.path(projectpath,"Aligned","SCAN_RANGE.Rdata"))
	load(file.path(projectpath,"HMCR",type,"MVA_DATA.Rdata"))
	load(file.path(projectpath,"HMCR","MCR.Rdata"))
	
 	excludedfiles <-  eval(parse(text=paste("MCR$",tolower(type),"$excluded",sep="")))
	
	load(file.path(projectpath,"Aligned","files.Rdata"))
	
	filenames <-  sub("(.+)[.][^.]+$", "\\1", basename(files))
	
	if(!all){
		cat("Select compounds to export.\n")
		choice	<-	which(VARID1 %in% select.list(VARID1,multiple=TRUE))
	}else
		choice  <-  1:length(VARID1)
	
	if(!length(choice)){
		cat("Error! No compounds selected!\n")
	
	}else{
		S		<-	t(SPECTRUM[,choice])
		RT		<-	VARID2[choice]
		NAME	<-	as.matrix(VARID1[choice])
		ID		<-	cbind(round(RT,2),NAME)
		
		#load(file.path(projectpath,"sampleinfo.Rdata"))
	 	
	 	if(file.exists(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata"))){
	 		load(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata"))
	 		RI	<-	approxExtrap(RT_INFO[,2],RT_INFO[,1],RT)$y
		
		}else
	 		RI	<-	rep(-99,length(VARID1[choice]))
	 		
	 	if(length(excludedfiles))
	 		row.names(DATA) <-  filenames[-excludedfiles]
	 	
		else
			row.names(DATA) <-  filenames
		
		out <-	data.frame(Window = VARID1[choice],RT_s = round(VARID2[choice],2),RI,Annotated = character(length(VARID1[choice])),Classified = character(length(VARID1[choice])),Comments = character(length(VARID1[choice])),t(DATA)[choice,])
				
		write.table(out,file=file.path(projectpath,"HMCR",type,"Spectra.txt"),quote=FALSE,sep="\t")
		

		if(file.exists(file.path(projectpath,"HMCR",type,"Spectrum.txt")))
			file.remove(file.path(projectpath,"HMCR",type,"Spectrum.txt"))
		
		fid	<-	file(file.path(projectpath,"HMCR",type,"Spectrum.txt"),"a+")
		
		for(j in 1:nrow(ID)){
			s				<-	S[j,]
			s				<-	s/max(s)*999
			s[round(s)<1.5]	<-	0
			N				<-	which(s == 0)
			
			if(length(N)){
				s	<-	s[-N]
    			MZ	<-	SCAN_RANGE[-N]
			
			}else
				MZ		<-	SCAN_RANGE
			DATA	<-	cbind(MZ,round(s))
			
			cat("Name:", paste(ID[j,],collapse=" "),"\n",file=fid,sep="")
			cat("Synon:Annotated:\n",file=fid,sep="")
			cat("Synon:Classified:\n",file=fid,sep="")
			cat("Synon:RI:", round(RI[j],1),"\n",file=fid,sep="")
			cat("Synon:RT(s):", round(RT[j],1),"\n",file=fid,sep="")
			cat("Synon:MCR-Reference:", VARID1[choice][j],"\n",file=fid,sep="")
			cat("Synon:Date:",as.character(Sys.Date()),"\n",file=fid,sep="")
			cat("Comments:","\n",file=fid,sep="")
			cat("DB#:",j,"\n",file=fid,sep="")
			cat("Num peaks:", length(MZ),"\n",file=fid,sep="")
			cat(file=fid,paste(apply(DATA,1,function(DATA) paste(DATA,collapse= " ")),";",rep(c(rep("",3),"\n"),len=length(MZ)),sep="",collapse=" "),"\n")
		}
		
		close(fid)
		cat("\nSpectrum exported!\n\n")
	}
}

