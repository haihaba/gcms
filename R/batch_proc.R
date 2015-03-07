##' Function batch_proc
##' 
##' Function batch_proc
##' @export
##' @param projectpath
##' @param cores
batch_proc <- function(projectpath,cores,parallel=TRUE){
  load(file.path(projectpath,"HMCR","MCR.Rdata"))
  winlist<-grep("[P]",MCR$windowlist)
  incl<-MCR$reg$incl
  pred<-MCR$reg$pred
    if(parallel){
      require(parallel)
      mclapply(winlist, function(i) find_spectrum(projectpath,incl,pred,i), mc.cores=cores)
    }else{
      lapply(winlist,function(i) find_spectrum(projectpath,incl,pred,i))
    }
  read_win(getwd(),'REG',winlist)
}


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

##' Function find_spectrum
##' 
##' Function find_spectrum
##' @param projectpath
##' @param incl
##' @param extra
##' @param windowlist
##' @return windowlist
##' @export
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

##' Function separate_windows
##' 
##' Function separate_windows
##' @param projectpath
##' @param windowlist
##' @param files
##' @param edges
##' @param model
##' @param prediction
separate_windows<-function(projectpath,windowlist,files,edges,model,prediction){
  
  if(length(windowlist)){
    load(file.path(projectpath,"maxMZ.Rdata"))
    load(file.path(projectpath,"Aligned","shift.Rdata"))
    dir.create(file.path(projectpath,"Edges","dat","Model samples"),showWarnings=FALSE,recursive=TRUE)
    dir.create(file.path(projectpath,"Edges","dat","Prediction samples"),showWarnings=FALSE,recursive=TRUE)
    dir.create(file.path(projectpath,"Edges","dat","Metadata"),showWarnings=FALSE,recursive=TRUE)
    cat("Checking integrity of existing window data files..\n")
    
    modelfilespath			<-  list.files(file.path(projectpath,"Edges","dat","Model samples"),full.names=TRUE)
    predictionfilespath <-  list.files(file.path(projectpath,"Edges","dat","Prediction samples"),full.names=TRUE)
    
    
    if(length(files) == length(c(model,prediction))){
      excludedfiles <-	numeric()
      
    }else{
      excludedfiles	<-	(1:length(files))[-c(model,prediction)]
    }
    
    ## Check which windows that needs to be separated
    newlist 	<-  logical(max(windowlist))
    metadata  <-  list.files(file.path(projectpath,"Edges","dat","Metadata"),full.names=TRUE) # Alla f<U+00F6>nster med metadata
    metadata  <-  metadata[as.numeric(substr(basename(metadata),4,6)) %in% windowlist] # Alla f<U+00F6>nster med metadata som ska kollas
    newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(metadata),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns metadata till
    newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(modelfilespath),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns modellfiler till
    newlist[windowlist[!(windowlist %in% as.numeric(substr(basename(predictionfilespath),4,6)))]]  <-  TRUE # F<U+00F6>nster som det ej finns prediktionsfiler till
    
    if(length(metadata)){
      n <-  0
      for(k in as.numeric(substr(basename(metadata),4,6))){
        n <-  n + 1
        load(metadata[n])
        if(!(all.equal(modelfiles,model) == TRUE & all.equal(predictionfiles,prediction) == TRUE))
          newlist[k]  <-  TRUE # Om metadata inte st<U+00E4>mmer <U+00F6>verens med nuvarande inst<U+00E4>llningar
      }
    }
    
    windowlist  <-  which(newlist)
    
    # Remove data files if they exists
    
    if(any(as.numeric(substr(basename(modelfilespath),4,6)) %in% windowlist))
      file.remove(modelfilespath[which(as.numeric(substr(basename(modelfilespath),4,6)) %in% windowlist)])
    
    if(any(as.numeric(substr(basename(predictionfilespath),4,6)) %in% windowlist))
      file.remove(predictionfilespath[which(as.numeric(substr(basename(predictionfilespath),4,6)) %in% windowlist)])
    ##
    
    cat("Windows to separate: [ ", windowlist," ]\n")
    cat("In total", length(windowlist),"windows\n")
    filenum 				<-  0
    modelfiles  		<-  model
    predictionfiles <-  prediction
    
    
    for(i in sort(c(model,prediction))){
      filenum <-  filenum +1
      cat("Loading ",basename(files[i]),paste("(",filenum,"/",length(c(model,prediction)),") ",sep=""))
      load(files[i])
      cat("... Exporting windows... \n")
      
      for(j in windowlist){
        if(nrow(Xbc)<(edges[j]+shift[i])){
          a <- 0
        }else if(nrow(Xbc)<(edges[j+1]+shift[i])){
          a <- (edges[j]+shift[i]):nrow(Xbc)
        }else
          a <- edges[j]:edges[j+1]+shift[i]  # Check so we are not exceeding matrix index
        
        Xbc <-  cbind(Xbc,matrix(0,nrow = nrow(Xbc),ncol = ifelse(maxMZ > ncol(Xbc),maxMZ-ncol(Xbc),0)))                                                           # Add columns if necessary so rows in outputfile will have same length
        x	<-	rbind(Xbc[a,],matrix(0,nrow=sum((0:(edges[j+1]+shift[i]-nrow(Xbc)))>0),ncol=ncol(Xbc)))
        write.table(t(x),file=file.path(projectpath,"Edges","dat",ifelse(i %in% model,"Model samples","Prediction samples"),paste("win",ifelse(j<=99,ifelse(j<=9,paste("00",j,sep=""),paste("0",j,sep="")),j),".dat",sep="")),append=TRUE,row.names=FALSE,col.names=FALSE)
      }
    }
    
    #All windows exported successfully; save metadata
    for(j in windowlist)
      save(modelfiles,predictionfiles,excludedfiles,file=file.path(projectpath,"Edges","dat","Metadata",paste("win",ifelse(j<=99,ifelse(j<=9,paste("00",j,sep=""),paste("0",j,sep="")),j),"_metadata.Rdata",sep="")))
    
  }else
    cat("No windows selected!\n")
}

##' Function do_AR_all
##' 
##' Function do_AR_all
##' Funciton to do Alternate Regression
##' @param x
##' @param projectpath
##' @param NL
##' @param RP
##' @param RT_LIMIT
##' @param SCAN_RANGE
##' @return C, S, INDEX, R2, noise, CP
do_AR_all <- function(x,projectpath,NL,RP,RT_LIMIT,SCAN_RANGE){
  
  var <-  dim(x)[1]
  mz  <-  dim(x)[2]
  obs <-  dim(x)[3]
  cat("\n[TIMEPOINTS,M/Z-Channels,Observations] = [",var,",",mz,",",obs,"]\n\n")
  R					<-	OK	<-	1
  R2					<-	0
  C					<-	S	<-	numeric()
  OK_out				<-	3
  X 					<-  t(matrix(aperm(x,c(2,1,3)),nrow=mz,ncol=var*obs))
  ind					<-	w_mor2(X)
  noise				<-	which(ind<NL)
  x[,noise,]			<-	0
  ssX					<-  sum(apply(x,3,ss))
  Xmean				<-	apply(x,2,function(x) apply(x,1,sum))
  rm(x)
  X[,noise]		<-	0
  max_R				<-	qr(X[,setdiff(1:length(ind),noise)])$rank
  S						<-	abs(pca(X,1)$p)
  vars				<-	which(as.logical(SCAN_RANGE))
  cat("-----------------------------\n")
  gc()
  
  while(OK_out>0){
    iter	<-	0
    C 		<-  numeric()
    if(max_R >= R){
      if(R > 1){
        p		<-	pure(Xmean[,vars],R,0.01)
        sp  <-  p$sp
        imp <-  p$imp
        S	<-	matrix(0,mz,R)
        
        for(i in 1:R)
          S[vars[imp[i]],i]	<-	1
      }
      
      R2			<-	0.001
      I			<-	numeric()
      dif			<-	1
      C			<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)
      
      
      
      C[C<0]		<-	0
      out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
      C	  		<-  out$C
      CP			<-  out$CP
      
      while(dif > 1e-6){
        r2			<-	R2
        iter		<-	iter+1
        tC  		<-  t(C)
        S			<-	t(ginv(tC%*%C,tol=.Machine$double.eps^2)%*%(tC%*%X))
        S[S<0]		<-	0
        S 			<-  apply(S,2,function(S) S/sum(S))
        
        if(any(is.na(S)))
          break
        
        C	<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)
        
        if(any(is.na(C)))
          break
        C[C<0]	<-	0
        out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
        C       <-  out$C
        CP      <-  out$CP
        R2		<-	1-ss(X-C%*%t(S))/ssX
        #				cat(iter,": ",R2,"| sum(C) = ",sum(C),"\n")
        
        if(iter == 50)
          dif <-  0
        
        else
          dif <-  abs(R2-r2)
      }
      
      if(any(is.na(S)) | any(is.na(C))){
        
        cat("Warning: NAs found.\n")
        C		<-	matrix(1,nrow(C),ncol(C))
        S		<-	matrix(0,nrow(S),ncol(S))
        R2		<-  0
      }
      
      INDEX	<-	matrix(0,obs,R)
      
      for(i in 1:R){
        c	<-	matrix(C[,i],var,obs)
        for(j in 1:obs){
          if(max(c[,j])>0)
            INDEX[j,i]  <-  min(which.max(c[,j]))
          else
            INDEX[j,i]  <-  -99
        }
      }
      
      if(sum(INDEX == -99) > 0){
        I <-  row(INDEX)[INDEX == -99]
        J <-  col(INDEX)[INDEX == -99]
        
        for(i in 1:length(J)){
          index				<-	INDEX[,J[i]]
          index				<-	index[-(index == -99)]
          INDEX[I[i],J[i]]	<-	median(index)
        }
      }
      
      I		<-  order(colMeans(INDEX))
      Y		<-  sort(colMeans(INDEX))
      C		<-	C[,I]
      S		<-	S[,I]
      CP		<-	CP[I]
      INDEX	<-	INDEX[,I]
      cat("Rank: ", R,"\n")
      
      if(R>1)
        PERMUTED_ORDER	<-	sum(diff(t(INDEX)) < -RP)
      else
        PERMUTED_ORDER	<-	0
      
      cat("PERMUTED_ORDER = ",PERMUTED_ORDER,"\n")
      cat("Number of iterations: ", iter,"\n")
      cat("R2X = ",round(R2,5),"\n")
      
    }else{
      R2				<-	0
      S				<- 	C 	<-  numeric()
      PERMUTED_ORDER	<-	99
    }
    
    s		<-	S
    c		<-	C
    r2		<-	R2
    
    if(r2>0){
      if(!PERMUTED_ORDER){
        C_out		<-	c
        S_out		<-	s
        CP_out		<-	CP	
        R2_out		<-	r2
        OK_out		<-	3
        OK2_out		<-	0
        cat("OK\n")
        
      }else{
        OK2_out	<-	1
        cat("NOT OK [1]\n")
      }
      
    }else{
      OK2_out	<-	1
      cat("NOT OK [2]\n")
    }
    
    OK_out	<-	OK_out-OK2_out
    R		<-	R+1
    
    cat("Timestamp: ",format(Sys.time(), "%X"),"\n")
    cat("-----------------------------\n")
  }
  
  if(!exists("C_out")){
    C_out		<- numeric()
    S_out		<-	numeric()
    R2_out	<-	numeric()
    CP_out	<-	numeric()
  }
  
  C		<-	as.matrix(C_out)
  S		<-	as.matrix(S_out)
  R2		<-	R2_out
  CP		<-	CP_out
  R		<-	ncol(C)
  
  cat("FINAL MODEL\n")
  cat("Rank:\t",R,"\n")
  cat("R2X:\t",R2,"\n")
  cat("-----------------------------")
  
  
  if(!is.null(R)){ #produces errors when R=1
    #if(R>1){
    INDEX	<-	matrix(0,obs,R)
    
    for(i in 1:R){
      c	<-	matrix(C[,i],var,obs)
      for(j in 1:obs)
        INDEX[j,i]	<-	min(which.max(c[,j]))
      
    }
    
    I		<-	order(round(colMeans(INDEX)))
    INDEX 	<-  sort(round(colMeans(INDEX)))
  }else
    INDEX	<-	numeric()
  
  gc()
  out <-  list(C=C,S=S,INDEX=INDEX,R2=R2,noise=noise,CP=CP)
}

##' Function pure
##' 
##' Function pure
##' [sp,imp]=pure(d,nr,f)
##' sp purest row/column profiles
##' imp indexes of purest variables
##' d data matrix; nr (rank) number of pure components to search
##' if d(nspectra,nwave) imp gives purest nwave => sp are conc. profiles (nr,nspectra)
##' if d(nwave,nspectra) imp gives purest nspectra => sp are spectra profiles (nr,nwave)
##' f percent of noise allowed respect maximum of the average spectrum given in % (i.e. 1% or 0.1%))
##' @param d
##' @param nr
##' @param f
##' @return sp, imp
pure <-function(d,nr,f){
  
  #[nrow,ncol]=size(d);
  # calculation of the purity spectrum
  w   		<-  p	<-	s	<-	matrix(0,nrow=max(nr,1),ncol=ncol(d))
  f			<-	f/100
  s[1,]		<-	apply(d,2,sd)
  m			<-	colMeans(d)
  ll			<-	s[1,]^2+m^2
  f			<-	max(m)*f
  p[1,]		<-	s[1,]/(m+f)
  
  imp 		<-  numeric()
  mp 			<-  max(p[1,])
  imp[1]  <-  which.max(p[1,])
  
  # calculation of the correlation matrix
  l	<-	sqrt(s[1,]^2+(m+f)^2)
  for(i in 1:nrow(d))
    d[i,]	<-	d[i,]/l
  
  c	<-	(t(d)%*%d)/nrow(d)
  
  # calculation of the weights
  # first weight
  w[1,]	<-	ll/l^2
  p[1,]	<-	w[1,]*p[1,]
  s[1,]	<-	w[1,]*s[1,]
  
  # next weights
  if(nr > 1){
    for(i in 2:nr){
      for(j in 1:ncol(d)){
        dm			<-	wmat(c,imp,i,j)
        w[i,j]	<-	det(dm)
        p[i,j]	<-	p[1,j]*w[i,j]
        s[i,j]	<-	s[1,j]*w[i,j]
      }
      
      mp[i] 	<-  max(p[i,])
      imp[i]  <-  which.max(p[i,])
    }
  }
  
  sp	<-	normv2(t(d[,imp[1:nr]]))
  purelist  <-  list(sp=sp,imp=imp)
}

##' Function normv2
##' 
##' Functino normv2
##' @param s
##' @return sn
normv2 <-function(s){
  
  if(nrow(s) == 1){
    sn  <-  s/sqrt(sum(s^2))
    
  }else
    sn  <-  apply(s,1,function(s) s/sqrt(sum(s^2)))
}

##' Function wmat
##' 
##' Function wmat
##' @param tempmat
##' @param imp
##' @param irank
##' @param jvar
##' @return dm 
wmat<-function(tempmat,imp,irank,jvar){
  dm								<-	matrix(0,irank,irank)
  dm[1,1]							<-	tempmat[jvar,jvar]
  dm[1,2:irank]					<-	tempmat[jvar,imp[1:(irank-1)]]
  dm[2:irank,1]					<-	tempmat[imp[1:(irank-1)],jvar]
  dm[2:irank,2:irank] <-  tempmat[imp[1:(irank-1)],imp[1:(irank-1)]]
  
  dm
  
  
}

##' Function unimodal2
##' 
##' Function unimodal2
##' @param C
##' @param R
##' @param var
##' @param obs
##' @param RT_LIMIT
##' @return C, CP
unimodal2<-function(C,R,var,obs,RT_LIMIT){
  
  N 	<-  nrow(C)
  CP  <-  numeric(R)
  for(i in 1:R){
    c       <-  matrix(C[,i],var,obs)
    out		<-	unimodal3(c,obs,RT_LIMIT)
    C[,i]	<-	matrix(out$C,N,1)
    CP[i]	<-	out$cp
  }
  out <-  list(C=C,CP=CP)
}

##' Function w_mor2
##' 
##' Function w_mor2
##' @param X
##' @return ind
w_mor2<-function(X){
  
  X								<-	sweep(X,2,colMeans(X))
  ind								<-  apply(X,2,function(X) sqrt(sum(X^2))/sqrt(sum(diff(X)^2)))
  ind[!is.finite(ind)]  <-  0
  ind
}

##' Function unimodal3
##' 
##' Function unimodal3
##' @param C
##' @param obs
##' @param RT_LIMIT
##' @return C, CP
unimodal3<-function(C,obs,RT_LIMIT){
  
  #cp	<-	excel_round(median(which.max(colMeans(t(C)))))
  cp  <-	round(median(which.max(colMeans(t(C)))))
  if(!length(cp)){
    C		<-	C*0
    cp		<-	0
    
  }else{
    for(i in 1:obs){
      xpeak	<-	peak_pick(t(C[,i]))$xpeak
      mp		<-	which(as.logical(xpeak))
      if(length(mp)){
        if(which.min(abs(mp-cp)))
          mp	<-	mp[which.min(abs(mp-cp))]
        if(abs(mp-cp)<RT_LIMIT){
          D			<-	diff(C[1:mp,i])
          if(any(D<0))
            poss	<-	max(which(D<0))
          
          else
            poss	<-	1
          
          for(j in poss:1)
            C[j,i]	<-	min(C[j:mp,i])
          
          D			<-	diff(C[mp:nrow(C),i])
          
          if(length(which(D>0)))
            poss	<-	min(which(D>0))
          
          else
            poss	<-	FALSE
          
          if(poss)
            for(j in (mp-1+poss):nrow(C))
              C[j,i]	<-	min(C[mp:j,i])
          
          # Correction of strange peaks
          knew						<-  C[-which(diff(C[,i]) == 0),i]
          mpnew						<-	which.max(knew)
          k							<-	C[,i]*0
          k[1:length(knew)+mp-mpnew]	<-	knew
          C[,i]	<-	k
          
        }else
          C[,i]	<-	C[,i]*0
      }else
        C[,i]	<-	C[,i]*0
    }
  }
  out <-  list(C=C,cp=cp)
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

##' Function ss
##' 
##' Function ss
##' @param X
##' @return SSX
ss<-function(X)
  SSX	<-	sum(X^2)


