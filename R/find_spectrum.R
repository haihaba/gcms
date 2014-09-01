##' Funcion find_spectrum
##' 
##' Function find_spectrum
##' @export
##' @param projectpath
##' @param incl
##' @param extra
##' @param windowlist
##' @return windowlist
find_spectrum<-function(projectpath,incl,extra,windowlist){
	
	require(MASS)
	X11.options(type="cairo")
	load(file.path(projectpath,"Aligned","files.Rdata"))
	load(file.path(projectpath,"Edges","edges.Rdata"))
	load(file.path(projectpath,"Aligned","shift.Rdata"))
	load(file.path(projectpath,"Aligned","SCAN_RANGE.Rdata"))
	load(file.path(projectpath,"maxMZ.Rdata"))
	load(files[which.min(shift)])
	EDGES_TIME	<-	SCAN_INFO[,2]
	
	dir.create(file.path(projectpath,"HMCR","REG","win_png"),showWarnings = FALSE,recursive = TRUE)
	dir.create(file.path(projectpath,"HMCR","REG","win"),showWarnings = FALSE,recursive = TRUE)
	load(file.path(projectpath,"SETTINGS.Rdata"))
	
	NL				<-	SETTINGS$NL
	RP				<-	SETTINGS$RP
	RT_LIMIT		<-	SETTINGS$MPS
	DO_BL			<-	SETTINGS$BC2
	color			<-	rep(c("red","green","blue","black","purple","grey","yellow4"),20)

	rm(Xbc,SETTINGS)
	gc()
  
	dir.create(file.path(projectpath,"Edges","dat","Model samples"),showWarnings = FALSE,recursive=TRUE)
	dir.create(file.path(projectpath,"Edges","dat","Model samples_bg_corr"),showWarnings = FALSE,recursive=TRUE)
	dir.create(file.path(projectpath,"Edges","dat","Prediction samples"),showWarnings = FALSE,recursive=TRUE)
 	
	temp1  					<-  list.files(file.path(projectpath,"Edges","dat","Model samples"),full.names=TRUE)
	windowdatapath1	<-  temp1[as.numeric(substr(basename(temp1),4,6)) %in% windowlist]
	temp2  					<-  list.files(file.path(projectpath,"Edges","dat","Prediction samples"),full.names=TRUE)
	windowdatapath2 <-  temp2[as.numeric(substr(basename(temp2),4,6)) %in% windowlist]
	metadatapath    <-  list.files(file.path(projectpath,"Edges","dat","Metadata"),full.names=TRUE)
	metadata				<-	as.numeric(substr(basename(metadatapath),4,6))
	ok  <-  0
 	
 	if(length(windowlist)){
		counter <-  0
 		
 		for(win in windowlist){
 			counter <-  counter + 1
			
			if(win %in% metadata){ #Check integrity of files. Unless ok, re-separate windows
				
				load(metadatapath[which(metadata %in% win)])
				ok	<-	ifelse(all.equal(modelfiles,incl) == TRUE & all.equal(predictionfiles,extra) == TRUE,1,0)
				ok	<-	ifelse(file.exists(windowdatapath1[counter]),1,0)
				ok	<-	ifelse(length(extra),ifelse(file.exists(windowdatapath2[counter]),1,0),ok)
			
			}else
				ok  <-  0
			
			if(!ok){
				cat("Window data file not found.\n")
				separate_windows(projectpath,windowlist,files,edges,incl,extra)
				temp1			<-	list.files(file.path(projectpath,"Edges","dat","Model samples"),full.names=TRUE)
				windowdatapath1	<-	temp1[as.numeric(substr(basename(temp1),4,6)) %in% windowlist]
				temp2			<-	list.files(file.path(projectpath,"Edges","dat","Prediction samples"),full.names=TRUE)
				windowdatapath2	<-	temp2[as.numeric(substr(basename(temp2),4,6)) %in% windowlist]
				metadatapath	<-	list.files(file.path(projectpath,"Edges","dat","Metadata"),full.names=TRUE)
				metadata		<-	as.numeric(substr(basename(metadatapath),4,6))
			}
			
			cat("Window ",win," of ",length(edges)-1," (processing ",length(windowlist)," window(s) in total)\n")
			cat("Loading model files... ")
			
			xIncl <-  array(scan(file=windowdatapath1[counter],what=numeric(0),skip=0,quiet=TRUE),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(incl)))
			
			if(DO_BL){ # Removes baseline for model samples
				cat("Removing baseline.. ")
				BL				<-	array(apply(xIncl,3,function(xIncl) apply(xIncl,2,function(xIncl) approx(c(1,length(xIncl)),xIncl[c(1,length(xIncl))],1:length(xIncl),ties = "ordered")$y)),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(incl)))       # method = 'linear' by default
        		xIncl			<-	xIncl-BL
        		xIncl[xIncl<0] 	<-	0
        		save(BL,file=file.path(projectpath,"Edges","dat","Model samples_bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
				rm(BL)
			}
			
			cat("Done\n")
			invisible(gc())
			
			out 	<-  do_AR_all(xIncl,projectpath,NL,RP,RT_LIMIT,SCAN_RANGE)
			rm(xIncl)
			invisible(gc())
			C		<-	out$C
			S		<-	out$S
			INDEX	<-	out$INDEX
			R2		<-	out$R2
			noise	<-	out$noise
			CP		<-	out$CP
			
			if(length(C)){
				if(length(extra)){
					cat("\nLoading prediction files..\n")
					xExtra <-  try(array(scan(file=windowdatapath2[counter],what=numeric(0),skip=0,quiet=TRUE),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(extra))),silent=TRUE)
					errorcounter  <-  0
					
					while(mode(xExtra) == "character"){
						errorcounter  <-  errorcounter + 1
						cat(xExtra[1])
						cat("Garbage collection... retrying..\n")
						gc()
						xExtra <-  try(array(scan(file=windowdatapath2[counter],what=numeric(0),skip=0,quiet=TRUE),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(extra))),silent=TRUE)
						if(errorcounter > 4 & mode(xExtra) == "character")
							stop("Error!")
					}
					
					if(DO_BL){ # Removes baseline for prediction files
						
						cat("Removing baseline.. \n")
						BL				<-	try(array(apply(xExtra,3,function(xExtra) apply(xExtra,2,function(xExtra) approx(c(1,length(xExtra)),xExtra[c(1,length(xExtra))],1:length(xExtra))$y)),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(extra))),silent=TRUE)       # method = 'linear' by default
						errorcounter	<-	0
						
						while(mode(BL) == "character"){
							errorcounter  <-  errorcounter + 1
							cat(BL[1])
							cat("Garbage collection... retrying..\n")
							gc()
							BL	<-	try(array(apply(xExtra,3,function(xExtra) apply(xExtra,2,function(xExtra) approx(c(1,length(xExtra)),xExtra[c(1,length(xExtra))],1:length(xExtra))$y)),dim=c(length(edges[win]:edges[win+1]),maxMZ,length(extra))),silent=TRUE)       # method = 'linear' by default
							
							if(errorcounter > 4 & mode(BL) == "character")
								stop("Error!")
						}
						
						xExtra				<-	xExtra-BL
						xExtra[xExtra<0] 	<-	0
						rm(BL)
					}
					
					C2  <-  numeric()
					apply(xExtra,3,function(xExtra){ C2 	<<-	rbind(C2,do_AR_all_prediction(xExtra,S,CP,RT_LIMIT))})
					rm(xExtra)
					C	<-	CC	<-	rbind(C,C2)

				}else{
					cat("No prediction samples found..\n")
					CC  <-  C
				}
				
				Hz	<-	1/median(diff(SCAN_INFO[,2]))
				
				for(i in 1:ncol(CC)){
					c <-  matrix(CC[,i],length(edges[win]:edges[win+1]),length(incl)+length(extra))
					
					if(length(EDGES_TIME) < edges[win+1])
						EDGES_TIME <-  c(EDGES_TIME,approxExtrap(1:length(EDGES_TIME),EDGES_TIME,xout=(length(EDGES_TIME)+1):edges[win+1])$y)
					
					#for(j in 1:ncol(c)){
					#	plot(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz,c[,j],col=color[i],xlab="Time",ylab="Intensity",type="l",ylim=cbind(0,max(CC)),xlim=range(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz),main=paste("Window:",win),lwd=1)
					#	par(new=TRUE)
					#}
					
					if(i == 1)
						C <- as.matrix(colSums(c))
					else
						C <- cbind(C,colSums(c))
     			}
				
				par(new=FALSE)
				TIME	<-	(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz)[INDEX]
				cat("TIMES = ",TIME,"\n")
				cat("===============================\n")
				CC	<-	ceiling(CC)
				C	<-	ceiling(C)
				C	<-	C[order(c(incl,extra)),]
				
			}else
				TIME	<- 	CC		<-	numeric()
			
			MZ		<-	SCAN_RANGE
    	
    			save(C,CC,S,MZ,INDEX,TIME,noise,CP,file=file.path(projectpath,"HMCR","REG","win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
    		invisible(gc())
   	  	}
   	
   	}else
   		cat("You must select at least one window to process!\n")
	return(windowlist)
}

