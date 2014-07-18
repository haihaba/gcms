##' Function remove_outliers
##' 
##' Function remove_outliers
##' @param projectpath
##' @param outliers
remove_outliers<-function(projectpath,outliers){
	
	OK	<-	menu(c("Yes","No"),title="Remove outliers from further analysis?")
	if(OK == 1 & length(outliers) & is.numeric(outliers)){
		#do_log(projectpath,"Removing outliers: ")
		cat("Removing outliers...\n\n")
    	
    	load(file.path(projectpath,"Aligned","TIC.Rdata"))
		load(file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
		load(file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
		load(file.path(projectpath,"Aligned","shift.Rdata"))
		load(file.path(projectpath,"Aligned","files.Rdata"))
		load(file.path(projectpath,"Aligned","COLUMNID1.Rdata"))

		#do_log(projectpath,paste(files[OUTLIER],"\n"))
		#do_log(projectpath,"----------------------")

		files			<-	files[-outliers]
		shift			<-	shift[-outliers]
		BASEPEAK	<-	BASEPEAK[-outliers,]
		TIC				<-	TIC[-outliers,]
		SUM_MZ		<-	SUM_MZ[-outliers,]
		COLUMNID1	<-	COLUMNID1[-outliers]

		save(TIC,file=file.path(projectpath,"Aligned","TIC.Rdata"))
		save(BASEPEAK,file=file.path(projectpath,"Aligned","BASEPEAK.Rdata"))
		save(SUM_MZ,file=file.path(projectpath,"Aligned","SUM_MZ.Rdata"))
		save(shift,file=file.path(projectpath,"Aligned","shift.Rdata"))
		save(files,file=file.path(projectpath,"Aligned","files.Rdata"))
		save(COLUMNID1,file=file.path(projectpath,"Aligned","COLUMNID1.Rdata"))
	
	}else
		cat("Aborted by user or no outliers selected!\n")

}

