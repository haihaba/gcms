##' Function new_samples
##' 
##' Function new_samples
##' @export
##' @param projectpath
new_samples<-function(projectpath){
  require(tcltk)
	
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
				#dir.create(file.path(predpath,"ALKANE_SERIE"),showWarnings=FALSE)
				dir.create(file.path(predpath,"HMCR"),showWarnings=FALSE)
				
				if(file.exists(file.path(projectpath,"HMCR","MCR.Rdata")))
					file.copy(file.path(projectpath,"HMCR","MCR.Rdata"),file.path(predpath,"HMCR"))
				
				#if(file.exists(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata")))
			  #		file.copy(file.path(projectpath,"ALKANE_SERIE","RT_INFO.Rdata"),file.path(predpath,"ALKANE_SERIE"))
				
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
		
	
  
  
	### Menu
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
			DATA	<-	sigma_unfold(predpath)
			scoresAR(predpath,DATA)
				
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


##' Function scoresAR
##' 
##' Function scoresAR
##' @param predpath
##' @param DATA
##' @return DATA, SPECTRUM, VARID1, OBSID1, VARID2
scoresAR  <-	function(predpath,DATA){
  
  if(missing(DATA))
    load(file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))
  
  load(file.path(predpath,"info.Rdata"))
  load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
  D					<-	numeric()
  VARID1		<-	character()
  spectrum	<- numeric()
  for(i in 1:max(ROWID1[,1]))
  {
    textstr	<-	paste("Checking window number:", i)
    cat(textstr,"\n")
    #do_log(predpath,textstr)
    int   <-   which(ROWID1[,1]==i)
    XX		<-	DATA[,int]    # DATA = X
    SSX		<-	sum(XX^2)
    
    if(min(dim(XX))>2)
    {
      CS					<-	ARmethod1(XX,predpath)
      comp				<-	ncol(CS$C)
      D						<-	cbind(D,CS$C)
      S_temp			<-	matrix(0,nrow=comp,ncol=max(SCAN_RANGE))
      MZ					<-	ROWID1[int,2]
      S_temp[,MZ]	<-	t(CS$S)
      spectrum		<-	rbind(spectrum,S_temp)
      for(j in 1:comp)
      {
        if(i > 99)
        {
          if(j < 10)
          {
            VARID1	<-	rbind(VARID1,paste("W",i,"_C0",j))
          }else
            VARID1	<-	rbind(VARID1,paste("W",i,"_C",j))
        }else
        {
          if(i < 10)
          {
            if(j < 10)
            {
              VARID1	<-	rbind(VARID1,paste("W00",i,"_C0",j))
            }else
              VARID1	<-	rbind(VARID1,paste("W00",i,"_C",j))
          }else
            if(j < 10)
            {
              VARID1	<-	rbind(VARID1,paste("W0",i,"_C0",j))
            }else
              VARID1	<-	rbind(VARID1,paste("W0",i,"_C",j))
        }
      }
    }
  }
  #output	<-	list(D=D,VARID1=VARID1,spectrum=spectrum)
  DATA  <-  D
  load(file.path(predpath,"Aligned","COLUMNID1.Rdata"))
  OBSID1	<-	COLUMNID1
  save(DATA,OBSID1,VARID1,file=file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))
  write.table(data.frame(ID=OBSID1,DATA),file=file.path(predpath,"HMCR","REG","MVA_DATA.txt"),row.names=FALSE,col.names=c("PrimaryID",VARID1),sep="\t",quote=FALSE)
  save(spectrum,VARID1,file=file.path(predpath,"HMCR","REG","SPECTRUM.Rdata"))
  
  
  load(files[which.min(shift)])
  EDGES_TIME	<-	SCAN_INFO[,2]
  save(EDGES_TIME,file=file.path(predpath,"EDGES_TIME.Rdata"))
  
}



##' Function ARmethod1
##' 
##' Function ARmethod1
##' @param X
##' @param projectpath
##' @return C, S
ARmethod1 <- function(X,projectpath){
  require(MASS)
  X[X<0]	<-	0
  XX			<-	sweep(X,2,colMeans(X))
  vec    <-  pca(XX,min(c((dim(X)-1),15)))$vec
  SSX			<-	sum(XX^2)
  load(file.path(projectpath,"SETTINGS.Rdata"))
  limit		<-	SETTINGS$R2Xc
  R				<-	sum(cumsum(diag(t(vec)%*%vec)/SSX)<limit)+1
  # R f?r ej vara 0
  
  p		<-	pca(X,R)$p
  
  S				<-	abs(p)
  Sstart	<-	S
  diff		<-	1
  R2			<-	0
  z				<-	0
  zz			<-	100
  C				<-	X%*%S%*%ginv(t(S)%*%S)
  C[C<0]	<-	0
  Cstart	<-	C
  while(diff > 1e-6)
  {
    z				<-	z+1
    S				<-	t(X)%*%C%*%ginv(t(C)%*%C)
    S[S<=0]	<-	0
    if(R == 1) # Applyfunktionen nedan fungerar endast f?r R>1
    {
      S[,1]/sqrt(sum(S[,1]^2))
    }else
      S[,1:R] 	<-  apply(S[,1:R],2,function(S) S/sqrt(sum(S^2)))
    
    CC	<-	t(S)%*%S
    if(sum(CC>0.95) > R)
    {
      zz	<-	z
      z		<-	100
    }
    C				<-	X%*%S%*%ginv(t(S)%*%S)
    C[C<0]	<-	0
    R2new		<-	1-ss(X-C%*%t(S))/ss(X)
    diff		<-	abs(R2-R2new)
    R2			<-	R2new
    if(any(is.na(R2)))
      z	<-	100
    if(z == 100)
    {
      cat("Rank reduced after: ",zz, " iterations.\n")
      z			<-	0
      R			<-	R-1
      C			<-	Cstart[,1:R]
      diff	<-	1
      R2		<-	0
      zz  	<- 100
    }
  }
  cat("Rank:       ",R,"\n")
  cat("R2X:        ",R2,"\n")
  cat("Iterations: ",z,"\n")
  output	<-	list(C=C,S=S)
}


##' Function ss
##' 
##' Function ss
##' @param X
##' @return SSX
ss<-function(X)
  SSX	<-	sum(X^2)


##' Function pca
##' 
##' Function pca
##' @export
##' @param X
##' @param comp
##' @return vec, p
pca<-function(X,comp){
  
  #   Calculates Principal componets of X by calc eigvectors of X'*X or X*X'
  #   Depending on whats easiest to calculate....
  
  if(nrow(as.matrix(X)) < ncol(as.matrix(X))){
    tX  <-  t(X)
    p 	<-  eigen(X%*%tX)$vectors[,1:comp]
    
    if(comp>1){
      p <-  apply(p,2,function(p) tX%*%p)
      p <-  apply(p,2,function(p) p/sqrt(sum(p^2)))
      
    }else{
      p	<-	tX%*%p
      p 	<-	p/sqrt(sum(p^2))
    }
    
  }else
    p		<-	eigen(t(X)%*%X)$vectors[,1:comp]
  
  vec		<-	X%*%p
  vecp	<-	list(vec=vec,p=p)
}


##' function aligndatapred
##' 
##' Function aligndatapred
##' 
##' @param projectpath
##' @param predpath
##' @param newsamples
##' @param target
##' @param minshift
##' @param datasource
##' @param ion
aligndatapred <- function(projectpath,predpath,newsamples,target,minshift,datasource,ion){
  
  load(file.path(projectpath,"maxMZ.Rdata"))
  alignfiles  <-  newsamples[order(newsamples)]
  
  A_DATA  			<-  list()
  #do_log(predpath,"Aligning data")
  n 						<-  length(alignfiles)
  dir.create(file.path(predpath,"Aligned"),showWarnings=FALSE)
  COLUMNID1 		<-  cbind(sub("[.][^.]*$", "", basename(alignfiles)))   # Filenames without extensions
  A_DATA$files  <-  files	<-	alignfiles
  A_DATA$fsize  <-  matrix(0,n,2)
  
  for(i in 1:n){
    cat("Loading ", paste(basename(alignfiles[i])," (",i,"/",n,")\n",sep=""))
    load(alignfiles[i])
    A_DATA$fsize[i,]		<-	dim(Xbc)
    num									<-	which(SCAN_RANGE == ion)
    
    if(i == 1){
      NUM_MZ			<-	ncol(Xbc)
      NUM_scans		<-	nrow(Xbc)
      A_DATA$tic	<-	matrix(0,n,NUM_scans)
      A_DATA$bp		<-	matrix(0,n,NUM_scans)
      A_DATA$ic		<-	matrix(0,n,NUM_scans)
      SUM_MZ      <-  matrix(0,n,NUM_MZ)
    }
    
    if(nrow(Xbc) > NUM_scans){
      A_DATA$tic	<-	cbind(A_DATA$tic,matrix(0,n,nrow(Xbc)-NUM_scans))
      A_DATA$bp		<-	cbind(A_DATA$bp,matrix(0,n,nrow(Xbc)-NUM_scans))
      A_DATA$ic		<-  cbind(A_DATA$ic,matrix(0,n,nrow(Xbc)-NUM_scans))
      NUM_scans		<-	nrow(Xbc)
    }
    
    if(ncol(Xbc) > NUM_MZ)
      SUM_MZ      <-  cbind(SUM_MZ,matrix(0,n,ncol(Xbc)-NUM_MZ))
    
    A_DATA$tic[i,1:nrow(Xbc)]		<-	rowSums(Xbc)
    A_DATA$bp[i,1:nrow(Xbc)]		<-	apply(Xbc,1,max)
    A_DATA$ic[i,1:nrow(Xbc)]		<-	t(Xbc[,num])
    SUM_MZ[i,1:ncol(Xbc)]				<-	colSums(Xbc)
    
    if(ncol(Xbc) < maxMZ){
      Xbc					<-	cbind(Xbc,matrix(0,nrow = nrow(Xbc),ncol = maxMZ-ncol(Xbc)))
      SCAN_RANGE  <-  c(SCAN_RANGE,(max(SCAN_RANGE)+1):maxMZ)
      save(Xbc,SCAN_INFO,SCAN_RANGE,file = alignfiles[i])
      
    }else if(ncol(Xbc) >= maxMZ){
      Xbc					<-	Xbc[,1:maxMZ]
      SCAN_RANGE  <-  SCAN_RANGE[1:maxMZ]
      save(Xbc,SCAN_INFO,SCAN_RANGE,file = alignfiles[i])
    }
  }
  
  A_DATA$ion	<- ion
  A_DATA$ok		<-  1
  
  if(file.exists(file.path(predpath,"SETTINGS.Rdata"))){
    load(file.path(predpath,"SETTINGS.Rdata"))
    max_shift	<-	SETTINGS$MS
    
  }else{
    cat("No settings for maximum shift found! Using default value (500)")
    max_shift <-  500
  }
  
  if(datasource == "TIC"){
    DATA  <-  A_DATA$tic
    
  }else if(datasource == "Basepeak"){
    DATA	<-	A_DATA$bp
    
  }else if(datasource == "IC")
    DATA	<-  A_DATA$ic
  
  if(length(target) < ncol(DATA))
    target  <-  c(target,rep(0,ncol(DATA)-length(target)))
  
  shift <-  numeric()
  
  for(i in 1:length(alignfiles)){
    CO_VAR	<-	numeric()
    A				<-	target[(max_shift+1):(ncol(DATA)-max_shift)]
    
    for(j in -max_shift:max_shift){
      B				<-	DATA[i,(max_shift+1+j):(ncol(DATA)-max_shift+j)]
      CO_VAR	<-	c(CO_VAR,sum(A*B))
    }
    
    shift	<-	c(shift,(which.max(CO_VAR)-max_shift-1))
    
  }
  
  shift	<-	shift-minshift
  cat("Min shift: ", min(shift),"\n")
  cat("Max shift: ", max(shift),"\n")
  
  if(min(shift)<0)
    cat("Adjusting negative shifts...\n")
  
  shift <-  shift-min(shift)
  
  
  
  start 		<-  1+shift[1:nrow(A_DATA$tic)]
  stop			<-	ncol(A_DATA$tic)-max(shift)+shift[1:nrow(A_DATA$tic)]
  BASEPEAK	<-	matrix(0,nrow(A_DATA$bp),ncol(A_DATA$bp))
  TIC				<-	matrix(0,nrow(A_DATA$tic),ncol(A_DATA$tic))
  
  for(i in 1:nrow(A_DATA$bp)){
    TIC[i,1:(stop[i]-start[i]+1)]				<-	A_DATA$tic[i,start[i]:stop[i]]
    BASEPEAK[i,1:(stop[i]-start[i]+1)]	<-	A_DATA$bp[i,start[i]:stop[i]]
  }
  
  TIC_RAW	<-	A_DATA$tic
  
  cat("Saving variables..\n")
  save(TIC_RAW,file=file.path(predpath,"Aligned","TIC_RAW.Rdata"))
  save(shift,file=file.path(predpath,"Aligned","shift.Rdata"))
  save(TIC,file=file.path(predpath,"Aligned","TIC.Rdata"))
  save(BASEPEAK,file=file.path(predpath,"Aligned","BASEPEAK.Rdata"))
  save(SUM_MZ,file=file.path(predpath,"Aligned","SUM_MZ.Rdata"))
  save(COLUMNID1,file=file.path(predpath,"Aligned","COLUMNID1.Rdata"))
  save(files,file=file.path(predpath,"Aligned","files.Rdata"))
  save(SCAN_RANGE,file=file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
  save(NUM_scans,file=file.path(predpath,"Aligned","NUM_scans.Rdata"))
  
  textstr <-  paste(files," Shift: ", shift," Scans: ", A_DATA$fsize[,1]," MZ: ",A_DATA$fsize[,2],sep="\n")
  #do_log(predpath,textstr)
}

##' Function read_win2
##' 
##' Function read_win2
##' @param predpath
##' @param type
##' @param win
##' @return DATA, SPECTRUM, VARID1, OBSID1, VARID2
read_win2<-function(predpath,type,win){
  
  #require(xlsReadWrite)
  SPECTRUM	<-	DATA	<-	VARID1	<-	VARID2	<-	numeric()
  load(file.path(predpath,"Edges","edges.Rdata"))
  num_windows	<-	length(edges)-1
  load(file.path(predpath,"Aligned","files.Rdata"))
  
  if(!missing(win)){
    
    if(length(list.files(file.path(predpath,"HMCR",type,"win"))) > length(win)){
      cat("\nRecently processed windows: ",win,"\n")
      cat("Previously processed windows: ",	setdiff(as.numeric(substr(basename(list.files(file.path(predpath,"HMCR",type,"win"))),4,6)),win),"\n")
      if(1 == menu(c("Yes","No"),title="There are some processed windows that are not selected, do you wish to include some of them as well?"))
        while(!length(win <-  sort(as.numeric(substr(basename(tk_choose.files(caption="Select windows to process")),4,6)))))
          cat("You must select at least one window!\n")
    }
    
  }else{
    cat("Select windows.")
    while(!length(win <-  sort(as.numeric(substr(basename(tk_choose.files(caption="Select windows to process")),4,6)))))
      cat("You must select at least one window!\n")
  }
  
  cat("\n")
  for(i in win){
    if(file.exists(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(i<=99,ifelse(i<=9,paste("00",i,sep=""),paste("0",i,sep="")),i),".Rdata",sep="")))){
      
      load(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(i<=99,ifelse(i<=9,paste("00",i,sep=""),paste("0",i,sep="")),i),".Rdata",sep="")))
      if(length(S)){
        cat(nrow(as.matrix(C)),"\n")
        if(nrow(as.matrix(C)) == length(files)){
          cat("--------------------------------------------------------\n")
          SPECTRUM	<-	cbind(SPECTRUM,S)
          DATA			<-	cbind(DATA,C)
          cat("Window num: ",i,"\nNumber of comps: ",ncol(S),"\nTotal number of comps: ", ncol(SPECTRUM),"\n")
          VARID2	<-	c(VARID2,TIME)
          
          for(j in 1:ncol(S))
            VARID1	<-	rbind(VARID1,ifelse(i>99,ifelse(j<10,paste("Win",i,"_C0",j,sep=""),paste("Win ",i,"_C",j,sep="")),ifelse(i<10,ifelse(j<10,paste("Win00",i,"_C0",j,sep=""),paste("Win00",i,"_C",j,sep="")),ifelse(j<10,paste("Win0",i,"_C0",j,sep=""),paste("Win0",i,"_C",j,sep="")))))
        }else
          
          cat("Number of files and rows in data matrix does not match!\n")
      }
    }
  }
  
  load(file.path(predpath,"Aligned","COLUMNID1.Rdata"))
  OBSID1		<-	COLUMNID1
  
  save(DATA,VARID1,VARID2,OBSID1,SPECTRUM,file=file.path(predpath,"HMCR",type,"MVA_DATA.Rdata"))
  save(SPECTRUM = SPECTRUM,file=file.path(predpath,"HMCR",type,"SPECTRUM.Rdata"))
  save(VARID1	= VARID1,file=file.path(predpath,"HMCR",type,"VARID1.Rdata"))
  save(VARID2 = VARID2,file=file.path(predpath,"HMCR",type,"VARID2.Rdata"))
  
  write.table(data.frame(ID = OBSID1,DATA),file=file.path(predpath,"HMCR",type,"MVA_DATA.txt"),sep='\t',row.names = FALSE,col.names = c("Primary ID",VARID1),quote=FALSE)
  out <-  list(DATA = DATA,SPECTRUM = SPECTRUM,VARID1 = VARID1, OBSID1 = OBSID1,VARID2 = VARID2)
}


##' Function sigma_unfold
##' 
##' Function sigma_unfold
##' @param predpath
##' @return DATA
sigma_unfold  <-  function(predpath){
  
  
  if(!file.exists(file.path(predpath,"Aligned","shift.Rdata"))){
    cat("Error! You need to align data first!\n")
    return(NULL)
    
  }else{
    load(file.path(predpath,"Aligned","files.Rdata"))
    load(file.path(predpath,"Aligned","COLUMNID1.Rdata"))
    load(file.path(predpath,"Aligned","shift.Rdata"))
    OBSID1	<-	COLUMNID1
    load(file.path(predpath,"SETTINGS.Rdata"))
    DOBL	<-	SETTINGS$BC2
    
    if(file.exists(file.path(predpath,"Edges","edges.Rdata"))){
      load(file.path(predpath,"Edges","edges.Rdata"))
      load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
      ROWID1	<-	varnamegen(length(edges)-1,length(SCAN_RANGE),min(SCAN_RANGE))
      load(file.path(predpath,"Aligned","NUM_scans.Rdata"))
      
      for(i in 1:length(files)){
        DATA_temp	<- numeric()
        cat("Loading ",basename(files[i]),paste("(",i,"/",length(files),")",sep=""),"\n")
        load(files[i])
        
        if(shift[i]+tail(edges,1)>nrow(Xbc))    # Avoids "subscript out of bounds" if the shift is too big (outliers)
          Xbc <-  rbind(Xbc,matrix(0,shift[i]+tail(edges,1)-nrow(Xbc),ncol(Xbc)))
        
        for(k in 1:(length(edges)-1)){
          X_temp <- Xbc[edges[k]:edges[k+1]+shift[i],]
          
          if(DOBL){
            BL               <- apply(X_temp,2,function(X_temp) approx(c(1,length(X_temp)),X_temp[c(1,length(X_temp))],1:length(X_temp))$y)       # method = 'linear' by default
            X_temp           <- X_temp-BL
            X_temp[X_temp<0] <- 0
          }
          
          DATA_temp <- cbind(DATA_temp, matrix(colSums(X_temp),nrow=1))
        }
        
        if(i == 1)
          DATA  <-  matrix(0,length(files),length(DATA_temp))
        
        if(ncol(DATA) < ncol(DATA_temp))
          DATA  <-  cbind(DATA,matrix(0,nrow=nrow(DATA),ncol=ncol(DATA_temp)-ncol(DATA)))
        
        DATA[i,1:ncol(DATA_temp)]	<-	DATA_temp
      }
      
      dir.create(file.path(predpath,"HMCR","REG"),showWarnings=FALSE)
      save(DATA,OBSID1,file=file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))
      save(edges,shift,ROWID1,files,file=file.path(predpath,"info.Rdata"))
      output <- DATA
    }else{
      cat("Error! You need to set edges first!\n")
      return(NULL)
    }
  }
}

##' Function varnamegen
##' 
##' Function varnamegen
##' @param intervals
##' @param peks
##' @param start
##' @return A
varnamegen <- function(intervals,peks,start)
{
  #A <- cbind(paste("int ",sort(rep(1:intervals,peks)),"  M/Z ", rep(start:(start+peks-1),intervals)," ",sep=""))
  A <-  matrix(c(sort(rep(1:intervals,peks)),rep(start:(start+peks-1),intervals)),peks*intervals,2)
}

##' Function find_spectrum2
##' 
##' Function find_spectrum2
##' @param predpath
##' @param projectpath
##' @return type, windwolist
find_spectrum2<-function(predpath,projectpath){
  
  require(MASS)
  #X11.options(type="cairo")
  datamenu  <-  character()
  
  if(file.exists(file.path(projectpath,"HMCR","REG","MVA_DATA.Rdata")))
    datamenu  <-  "REG H-MCR DATA"
  
  
  if(length(datamenu)){
    dataexport  <-  menu(c(datamenu,"Cancel"),title="Select data")
    if(dataexport == 0 | dataexport == length(c(datamenu,"Cancel")))
      datamenu <- character()
  }
  
  if(!length(datamenu)){
    cat("Error! (Aborted by user or no data found)\n\n")
    return(character())
    
  }else{
    type = ifelse(length(datamenu) == 1,strsplit(datamenu," ")[[1]][1],c("REG","CV")[dataexport])
    load(file.path(predpath,"Aligned","files.Rdata"))
    load(file.path(predpath,"Edges","edges.Rdata"))
    load(file.path(predpath,"Aligned","shift.Rdata"))
    load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
    load(file.path(predpath,"maxMZ.Rdata"))
    load(files[which.min(shift)])
    load(file.path(projectpath,"HMCR","MCR.Rdata"))
    
    windowlist		<-	as.numeric(substr(basename(list.files(file.path(projectpath,"HMCR",type,"win"))),4,6))
    EDGES_TIME		<-	SCAN_INFO[,2]
    
    load(file.path(projectpath,"SETTINGS.Rdata"))
    
    NL        <- SETTINGS$NL
    RT_LIMIT  <- SETTINGS$MPS
    DO_BL     <- SETTINGS$BC2
    #color     <- cbind("red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4","red","green","blue","black","purple","grey","yellow4")
    
    rm(Xbc,SETTINGS)
    
    dir.create(file.path(predpath,"Edges","dat"),recursive = TRUE,showWarnings = FALSE)
    temp <- list.files(file.path(projectpath,"Edges","dat"),full.names=TRUE)
    dir.create(file.path(predpath,"Edges","dat","Model samples_bg_corr"),showWarnings = FALSE,recursive=TRUE)
    dir.create(file.path(predpath,"HMCR",type,"win_png"),showWarnings = FALSE,recursive = TRUE)
    dir.create(file.path(predpath,"HMCR",type,"win"),showWarnings = FALSE,recursive = TRUE)
    
    gc()
    
    Scores  <-  numeric()
    
    ### attpempt to store the background correction data. bg correction storage in "find_spectrum.R" is done in arrays per window 
    ### scan * m/z * sample. Here processing is sample wise (in "find_spectrum.R", the whole sampleset is in memory at once),
    ### so the array has to be build step by step
    
    ### First array variables for each window has to be defined. Then in the two main loops "sample" and "window" the variables
    ### are filled sample by sample. After finishing the two loops, the arrays have to be stored in corresponding files 
    
    
    
    #for(win in windowlist){
    #	scans<-(edges[win+1]-edges[win]+1)
    #	mzs<-length(SCAN_RANGE)
    #	BL<-array(0,dim=c(scans,mzs,length(files)))
    #	#eval(parse(text=paste("BL_",win,"<-NULL",sep="")))
    #	save(BL,file=file.path(predpath,"Edges","dat","Model samples_bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
    #	rm(BL)
    
    #}
    
    if(length(windowlist)){
      for(i in 1:length(files)){
        load(files[i])
        cat("Resolving ",basename(files[i]),paste(" (",i,"/",length(files),")\n",sep=""))
        
        for(win in windowlist){
          
          if(edges[win]+shift[i] > 0 & edges[win+1]+shift[i] < nrow(Xbc)){
            x <-  Xbc[edges[win]:edges[win+1]+shift[i],]
            
            if(ncol(x) < maxMZ)
              x <-  cbind(x,matrix(0,nrow=nrow(x),ncol= maxMZ-ncol(x)))
            
            if(DO_BL){ # Removes baseline for prediction files
              BL<-matrix(apply(x,2,function(x) approx(c(1,length(x)),x[c(1,length(x))],1:length(x))$y),nrow=length(edges[win]:edges[win+1]),ncol=maxMZ)       # method = 'linear' by default
              x		<-	x-BL
              x[x<0]	<-	0
              
              ###test code###
              #BL_new<-BL
              #load(file.path(predpath,"Edges","dat","Model\ samples","bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
              #BL[,,i]<-BL_new
              #save(BL,file=file.path(predpath,"Edges","dat","Model samples","bg_corr",paste("bg_win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
              
              #rm(BL_new)
              rm(BL)
            }
            
            load(file.path(projectpath,"HMCR", "REG","win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
            # browser()
            if(!length(S))
              c2 <-  numeric()
            else
              c2 <-  do_AR_all_prediction(x,S,CP,RT_LIMIT)
            if(exists("PCApara"))
              scores  <-  (colSums(x)  <-  PCApara$mean)%*%PCApara$loadings[,1:2]
            else
              scores  <-  numeric()
            ok <- 1
            
          }else{
            load(file.path(projectpath,"HMCR", "REG","win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
            c2 <-  matrix(0,length(edges[win]:edges[win+1]),ncol(as.matrix(C)))
            ok  <-  0
            if(exists("PCApara"))
              scores <-  c(0,0)
            else
              scores  <-  numeric()
          }
          
          
          if(i == 1){
            if(ok){
              C<-matrix(colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))),nrow=1)
              CC  <-  c2
            }else{
              C<-matrix(colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))),nrow=1)*NA
              CC<-c2
            }
            
            Scores  <-  rbind(Scores,scores)
            save(C,CC,S,TIME,Scores,file=file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
          }else{
            load(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
            Scores  <-  rbind(Scores,scores)
            if(ok){
              C<-rbind(C,colSums(as.matrix(c2,ncol=ncol(as.matrix(C)))))
              CC<-rbind(CC,c2)
            }else{
              C<-rbind(C,colSums(as.matrix(c2,ncol=ncol(as.matrix(C))))*NA)
              CC<-rbind(CC,c2)
            }
            save(C,CC,S,TIME,Scores,file=file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
          }
          gc()
        }
        gc()
      }
      
      # Plot windows
      # for(win in windowlist){
      
      # if(win %in% as.numeric(substr(basename(list.files(file.path(projectpath,"HMCR",type,"win"))),4,6))){
      # rm(Scores)
      # load(file.path(predpath,"HMCR",type,"win",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),".Rdata",sep="")))
      
      # if(min(dim(CC)) > 0){
      # if(length(Scores)){
      # #Plot functions for Cross validation method
      # }
      
      # Hz		<-	1/median(diff(SCAN_INFO[,2]))
      
      # for(i in 1:ncol(CC)){
      # c2 <-  matrix(CC[,i],length(edges[win]:edges[win+1]),length(files))
      # if(length(EDGES_TIME) < edges[win+1])
      # EDGES_TIME <-  c(EDGES_TIME,approxExtrap(1:length(EDGES_TIME),EDGES_TIME,xout=(length(EDGES_TIME)+1):edges[win+1])$y)
      # for(j in 1:ncol(c2)){
      # plot(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz,c2[,j],col=color[i],xlab="Time",ylab="Intensity",type="l",ylim=cbind(0,max(CC)),xlim=range(EDGES_TIME[edges[win]:edges[win+1]]+median(shift)/Hz),main=paste("Window:",win),lwd=1)
      # par(new=TRUE)
      # }
      # }
      
      # par(new=FALSE)
      # #savePlot(file.path(predpath,"HMCR", "REG","win_png",paste("win",ifelse(win<=99,ifelse(win<=9,paste("00",win,sep=""),paste("0",win,sep="")),win),sep="")),type="png")
      
      # }else
      # cat("No data found in window ", win,"!\n")
      # }
      # }
      
    }else
      cat("No windows processed.\n")
    
    return(list(type=type,windowlist=windowlist))
  }
  
}

##' Function do_AR_all_prediction
##' 
##' Function do_AR_all_prediction
##' Funciton to do Alternate Regression on a prediction set
##' @param x
##' @param S
##' @param CP
##' @param RT_LIMIT
##' @return Cnew
do_AR_all_prediction <- function(x,S,CP,RT_LIMIT){
  
  var 		<-  nrow(x)
  mz  		<-  ncol(x)
  C			<-	x%*%S%*%ginv(t(S)%*%S)
  C[C<0]	<-	0
  Cnew		<-	C*0
  
  for(i in 1:ncol(C)){
    c		<-	as.matrix(C[,i],nrow=var)
    Cnew[,i]	<-	unimodal_3(c,1,CP[i],RT_LIMIT)
  }
  
  C	<-	Cnew
}

##' Function unimodal_3
##' 
##' Function unimodal_3
##' @param C
##' @param obs
##' @param cp
##' @param RT_LIMIT
##' @return C
unimodal_3<-function(C,obs,cp,RT_LIMIT){
  
  for(i in 1:obs){
    xpeak	<-	peak_pick(t(C[,i]))$xpeak
    mp		<-	which(as.logical(xpeak))
    
    if(!length(mp)){
      C[,i] <-  C[,i]*0
      
    }else{
      if(which.min(abs(mp-cp)))
        mp	<-	min(mp[which.min(abs(mp-cp))])
      else
        mp	<-	1
      
      if(abs(mp-cp)<RT_LIMIT){
        D	<-	diff(C[1:mp,i])
        
        if(length(which(D<0)))
          poss	<-	max(which(D<0))
        else
          poss	<-	1
        
        for(j in poss:1)
          C[j,i]	<-	min(C[j:mp,i])
        
        D		<-	diff(C[mp:nrow(C),i])
        poss	<-
          
          if(length(which(D>0)))
            poss	<-	min(which(D>0))
        else
          poss	<-	FALSE
        
        if(poss){
          for(j in (mp-1+poss):nrow(C))
            C[j,i]	<-	min(C[mp:j,i])
        }
        
        # Correction of strange peaks
        
        knew	<-	C[(diff(C[,i]) != 0),i]
        mpnew	<-	which.max(knew)
        k		<-	C[,i]*0
        k[1:length(knew)+mp-mpnew]	<-	knew
        C[,i]	<-	k
      }else
        
        C[,i]	<-	C[,i]*0
    }
  }
  C
}

##' Function peak_pick
##' 
##' Function peak_pick
##' Function to find peaks
##' @param x
##' @return xpeak, xout
peak_pick<-function(x){
  
  xout				<-	x	<-	as.matrix(x)
  
  xd1         <-  sgolayfilt(xout,3,11,1)
  NOISE_LIMIT	<-	median(x)
  N1					<-	which(x>NOISE_LIMIT)
  N2					<-	matrix(0,nrow(x),ncol(x))
  
  for (i in 3:(ncol(x)-2))
    if(xd1[i-2] > 0)
      if(xd1[i-1] > 0)
        if(xd1[i+1] < 0)
          if(xd1[i+2] < 0)
            if(sum(N1 == i) == 1)
              N2[i]	<-	TRUE
  
  N	<-	which(as.logical(N2))
  
  if(length(N) > 1){
    
    while(min(diff(N)) < 3){
      p1	<-	min(which(diff(N) < 3))
      p2	<-	p1+1
      
      if(xout[N[p1]] < xout[N[p2]])
        N	<-	N[-p1]
      else
        N	<-	N[-p2]
      if(length(N) == 1)
        break
    }
  }
  xpeak			<-	matrix(0,nrow(x),ncol(x))
  xpeak[N]	<-	xout[N]
  out <-  list(xpeak = xpeak, xout = xout)
}

##' Function sgolay
##' 
##' Function sgolay
##' @param p
##' @param n
##' @param ts
##' @return F 
sgolay<-function(p, n, m = 0, ts = 1){ 
  
  library(MASS)
  if (n %% 2 != 1)
    stop("sgolay needs an odd filter length n")
  
  if (p >= n)
    stop("sgolay needs filter length n larger than polynomial order p")
  
  F = matrix(0., n, n)
  k = floor(n/2)
  for (row  in  1:(k+1)) {
    C = ( ((1:n)-row) %*% matrix(1, 1, p+1) ) ^ ( matrix(1, n) %*% (0:p) )
    A = ginv(C)
    F[row,] = A[1+m,]
  }
  
  F[(k+2):n,] = (-1)^m * F[k:1,n:1]
  
  if (m > 0)
    F = F * prod(1:m) / (ts^m)
  
  class(F) = "sgolayFilter"
  F
}

##' Function sgolayfilt
##' 
##' Function sgolayfilt
##' @param x
##' @param p
##' @param n
##' @param m
##' @param ts
sgolayfilt<-function(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
  
  len = length(x)
  if (class(p) == "sgolayFilter" || (!is.null(dim(p)) && dim(p) > 1)){
    F = p
    n = nrow(F)
  }else
    F = sgolay(p, n, m, ts)
  k = floor(n/2)
  
  #z = filter(F[k+1,n:1], 1, x)
  
  filt  <-  F[k+1,n:1]
  z 		<-	na.omit(stats:::filter(c(rep(0,length(filt) - 1), x), filt, sides = 1))
  c(F[1:k,] %*% x[1:n], z[n:len], F[(k+2):n,] %*% x[(len-n+1):len])
}

