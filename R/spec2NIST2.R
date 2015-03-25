##' Function spec2NIST2
##' 
##' Function spec2NIST2
##' @param predpath
##' @param type
##' @param all
spec2NIST2<-function(predpath,type,all=FALSE){
	
	load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
	load(file.path(predpath,"HMCR",type,"MVA_DATA.Rdata"))
	load(file.path(predpath,"Aligned","files.Rdata"))
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
		
		#load(file.path(predpath,"sampleinfo.Rdata"))
		
	 	if(file.exists(file.path(predpath,"ALKANE_SERIE","RT_INFO.Rdata"))){
	 		load(file.path(predpath,"ALKANE_SERIE","RT_INFO.Rdata"))
	 		RI	<-	approxExtrap(RT_INFO[,2],RT_INFO[,1],RT)$y
		
		}else
			RI	<-	rep(-99,length(VARID1[choice]))
		
		row.names(DATA) <-  filenames
		out <-	data.frame(Window = VARID1[choice],RT_s = round(VARID2[choice],2),RI,Annotated = character(length(VARID1[choice])),Classified = character(length(VARID1[choice])),Comments = character(length(VARID1[choice])),t(DATA)[choice,])
			
		xlsORtxt  <-  menu(c("xls","txt"),title="Write output to excel or text?")
			
			if(xlsORtxt == 1){
				require(xlsReadWrite)
 			
 			write.xls(out,file=file.path(predpath,"HMCR",type,"Spectra.xls"))
		
		}else{
			write.table(out,file=file.path(predpath,"HMCR",type,"Spectra.txt"),quote=FALSE,sep="\t")
		}

		if(file.exists(file.path(predpath,"HMCR",type,"Spectrum.txt")))
			file.remove(file.path(predpath,"HMCR",type,"Spectrum.txt"))
	
		fid	<-	file(file.path(predpath,"HMCR",type,"Spectrum.txt"),"a+")
	
		for(j in 1:nrow(ID)){
			s					<-	S[j,]
	   		s					<-	s/max(s)*999
	    	s[round(s)<1.5]		<-	0
	    
	    	N					<-	which(s == 0)
	    
	    	if(length(N)){
	    		s	<-	s[-N]
	    		MZ	<-	SCAN_RANGE[-N]
	    
	    	}else
				MZ  <-  SCAN_RANGE
	    
	    	DATA	<-	cbind(MZ,round(s))
	    
	    	cat("Name:", paste(ID[j,],collapse=" "),"\n",file=fid,sep="")
	    	cat("Synon:Annotated:\n",file=fid,sep="")
	    	cat("Synon:Classified:\n",file=fid,sep="")
		    cat("Synon:RI:", round(RI[j],1),"\n",file=fid,sep="")
			cat("Synon:RT(s):", round(RT[j],1),"\n",file=fid,sep="")
			cat("Synon:MCR-Reference:", VARID1[choice][j],"\n",file=fid,sep="")
		    #cat("Synon:Protocol:", sampleinfo$PR,"\n",file=fid,sep="")
		    #cat("Synon:Project:", sampleinfo$PN,"\n",file=fid,sep="")
	    	#cat("Synon:PI:", sampleinfo$PI,"\n",file=fid,sep="")
		    #cat("Synon:Species:",sampleinfo$SP,"\n",file=fid,sep="")
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

##' Function approxExtrap
##' 
##' Function approxExtrap
##' This function is called from find_spectrum.R, find_spectrum2.R, spec2NIST.R and spec2NIST2.R
##' @param x
##' @param y
##' @param xout
##' @param method
##' @param n
##' @param rule
##' @param f
##' @param ties
##' @param na.rm
##' @return list(x,y)
approxExtrap <- function (x, y, xout, method = "linear", n = 50, rule = 2, f = 0,ties = "ordered", na.rm = FALSE){
  
  if(is.list(x)){
    y <- x[[2]]
    x <- x[[1]]
  }
  
  if(na.rm){
    d <- !is.na(x + y)
    x <- x[d]
    y <- y[d]
  }
  
  d <- !duplicated(x)
  x <- x[d]
  y <- y[d]
  d <- order(x)
  x <- x[d]
  y <- y[d]
  w <- approx(x, y, xout = xout, method = method, n = n, rule = 2,f = f, ties = ties)$y
  r <- range(x)
  d <- xout < r[1]
  
  if (any(is.na(d)))
    stop("NAs not allowed in xout")
  
  if (any(d))
    w[d] <- (y[2] - y[1])/(x[2] - x[1]) * (xout[d] - x[1]) + y[1]
  
  d <- xout > r[2]
  n <- length(y)
  
  if (any(d))
    w[d] <- (y[n] - y[n - 1])/(x[n] - x[n - 1]) * (xout[d] - x[n - 1]) + y[n - 1]
  
  list(x = xout, y = w)
}

