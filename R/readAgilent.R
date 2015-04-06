##' Function readAgilent
##' 
##' Function readAgilent
##' @export
##' @param projectpath
##' @param filepath     can be used to import a single file and return the data to the caller
##' @return Xbc, SCAN_INFO, SCAN_RANGE, file  
readAgilent <-function(projectpath,filepath){
  #require(ncdf)
  #require(tcltk)
  
  
  ### file selection by TclTk
  if(missing(filepath)){
    #cdffiles <- tk_choose.files(caption="Select .D directories to import.",multi=TRUE,filters = matrix(c(".D files (*.D)","*.D","all","*.*"),2,2,byrow=TRUE))
    allDirectories <- list.dirs(path=projectpath)[grep('*.D$',list.dirs(path=projectpath))]
    cdffiles <- select.list(choices = allDirectories, preselect = allDirectories, 
                            multiple = TRUE, graphics = TRUE, 'Choose .D Files')
  }
  else if(file.exists(filepath))
    cdffiles <- filepath
  
  if(length(cdffiles)){
    cat("Selected files: \n",paste(cdffiles,"\n ",collapse=""))
  } else {
    cat("No .D directores selected.")
    return(NULL)
  }
  
  
  
  ## create directory for .D raw data import
  dir.create(file.path(projectpath,"Filtered","CDF"),showWarnings=FALSE,recursive=TRUE)
  
  
  
  
  
  ## Settings Handling
  SETTINGS  <-  NULL
  
  
  ## Check if a setting file exists, if yes read 
  ## it else create by setting default values
  if(file.exists(file.path(projectpath,"SETTINGS.Rdata"))){
    load(file.path(projectpath,"SETTINGS.Rdata"))
    FL   <- SETTINGS$FL
    MZP  <- SETTINGS$MZP
    MZR  <- SETTINGS$MZR
  }else{
    cat("No filterlength settings found, using FL = 0.\n")
    FL  <-  0
    cat("No mass channel precision settings found, using MZP = 0.\n")
    MZP <-  0
    cat("No mass channel range settings found, using MZR min = 50, max = 800.\n")
    MZR <-  c(50,800)
  }
  
  
  ## Check if maximum MZ Rdata file exist, then 
  ## load, otherwise create the variable maxMZ
  if(file.exists(file.path(projectpath,"maxMZ.Rdata")))
    load(file.path(projectpath,"maxMZ.Rdata"))
  else
    maxMZ <- numeric()
  
  
  
  ## construct the filepaths for the files to save
  filepaths <-  file.path(projectpath,"Filtered","CDF",paste(sub("(.+)[.][^.]+$", "\\1", basename(cdffiles)),".Rdata",sep=""))
  cat("================================\n")
  
  
  
  
  
  ## looping through all the CDF files for importing
  for(i in 1:length(cdffiles)){
    ## ncdf file import
    errorcounter  <- 0
    cat("File ",basename(cdffiles[i]), "opened. ",paste("(",i,"/",length(cdffiles),")",sep=""),"\n")
    DATA       <- readDFile(cdffiles[i])
    TIME  <- as.numeric(rownames(DATA))
    importMz <- as.numeric(colnames(DATA))
    MZmin <- min(importMz)
    MZmax <- max(importMz)
    
    
    
    cat("File ",basename(cdffiles[i]), "closed.\n")
    
    ## Smoothing of the data, Xbc is the smoothed (moving average) data 
    ## matrix with rows mz and columns for data points
    cat(paste("Smoothing (FL = ",FL,")\n",sep=""))
    Xbc <- baseline(DATA,projectpath,FL)
    rm(DATA)
    
    ## prepare and save data
    
    ## SCAN_RANGE is the sequence of MZ's in the current dataset
    mz        <- SCAN_RANGE <- seq(max(1,round(MZmin)),min(ncol(Xbc),round(MZmax),MZR[2]))
    
    ## Reduce Xbc to SCAN_RANGE
    Xbc       <- Xbc[,mz]
    
    ## SCAN_INFO, time values from CDF import
    SCAN_INFO <- matrix(cbind(1:length(TIME),TIME),nrow=length(TIME),ncol=2)
    
    ## maxMZ, number of mz values 
    maxMZ     <- max(maxMZ,ncol(Xbc))
    
    ## Save data
    save(Xbc,SCAN_INFO,SCAN_RANGE,file = filepaths[i])
    rm(TIME)
    
    ## when function argument 'filepath' missing, delete the variables
    ## otherwise they are still needed to return to the caller
    if(missing(filepath))
      rm(Xbc,SCAN_INFO,SCAN_RANGE)
    
    cat("================================\n")
    
  }
  
  ### save maxMZ data
  save(maxMZ,file=file.path(projectpath,"maxMZ.Rdata"))
  
  
  ## return data in case filepath was provided as function argument
  if(missing(filepath)){
    cat("Done! \nAll files imported and stored in ",file.path(projectpath,"Filtered","CDF"),"\n\n")
  }else
    import = list(Xbc=Xbc,SCAN_INFO=SCAN_INFO,SCAN_RANGE=SCAN_RANGE,file=basename(filepath))    # Return variables Xbc, SCAN_INFO and SCAN_RANGE if a valid file path was supplied
}


##' Function Baseline
##' 
##' Function Baseline
##' @param X              raw imported Data. Time as rows and mz as columns
##' @param projectpath    baseproject directory
##' @param FL             Filterlength
##' @return Xbc           moving average filtered data table, one file
baseline <- function(X,projectpath,FL = 0){
  
  N <-  nrow(X)
  K <-  ncol(X)
  #X	<-	sweep(X,2,apply(X,2,min))
  
  if(FL){
    Xbc								<-	X*0
    Xbc[1:FL,]						<-	X[1:FL,]
    Xbc[(nrow(X)-FL+1):nrow(X),] 	<-	X[(nrow(X)-FL+1):nrow(X),]
    
    
    start<-FL+1
    end<-dim(X)[1]-FL
    
    ## moving average filter along rows (scantimes)
    Xbc[start:end,]<-filter(X,filter=rep(1/(FL*2+1),(FL*2+1)))[start:end,]
    
    
    if(any(Xbc<0))
      Xbc[Xbc<0]	<-	0
  }else
    Xbc <-  X
  return(Xbc)
}


##' Function readDFile
##' 
##' Function readDFile
##' @param pathname       the pathname of the directory containing the data to import
##' @return outData       list containing importInt, importMz, importCount, importScanTime
##' @export
readDFile<-function(pathname){
  
  filename<-file.path(pathname,'DATA.MS')
  
  cat('Opening Agilent file...\n')
  to.read<-file(filename,'rb')
  agilent<-readBin(to.read,integer(),size=1,signed=FALSE,n=20000000, endian='little')
  close(to.read)
  
  ### preparing vector with counts per scan
  cat('...extracting counts per scan...\n')
  countsNumber<-agilent[5782]
  counts<-agilent[5782]
  currentPosition<-5782
  
  while(currentPosition<length(agilent)){
    jumpLength<-countsNumber*4+7*4
    currentPosition<-currentPosition+jumpLength
    countsNumber<-agilent[currentPosition]
    counts<-c(counts,agilent[currentPosition])
  }
  
  
  ### counts will be too long and the last entry is NA
  ### counts will be corrected when the number of scans
  ### is known.
  #counts<-countNumber(agilent)
  ### cut away NA. 
  counts<-counts[-which(is.na(counts))]
  
  ### the second period is extracted. This one is currently
  ### used to determine the number of scans. As the useful length
  ### is not known yet, it is extracted in overlength. Na's have
  ### to be removed then
  cat('...determine number of scans...\n')
  secondPeriod<-agilent[seq(5772,2000000,4)]
  ### remove Na's
  secondPeriod<-secondPeriod[-which(is.na(secondPeriod))]
  
  
  ### inline function to extract the data paddings
  betweenSequence<-function(period,counts){
    tempSequence<-period[1:4]
    currentPosition<-4
    for(ii in counts){
      currentPosition<-(currentPosition+ii)
      tempSequence<-c(tempSequence,period[currentPosition:(currentPosition+7)])
      currentPosition<-(currentPosition+7)
    }
    return(tempSequence)  
  }
  
  ### the betweenSequence removes the variable length ion count
  ### data. Left are 8 numbers of padding betweeen each scan
  betweenSecond<-betweenSequence(secondPeriod,counts)
  ### in the betweenSequence of the second Period, the 8th
  ### field is currently used for determination of scan numbers
  ### when it is 3 times in sequence zero,
  numberOfScans<-which.max(abs(diff(betweenSecond[seq(1,100000,8)]))>1)
  counts<-counts[1:numberOfScans]
  rawExtractLength<-sum(counts)*4+(numberOfScans*4*7)+5770
  
  ### extract periods with correct length
  cat('...separate rawdata in four sequences...\n')
  firstPeriod<-agilent[seq(5771,rawExtractLength,4)]
  secondPeriod<-agilent[seq(5772,rawExtractLength,4)]
  thirdPeriod<-agilent[seq(5773,rawExtractLength,4)]
  fourthPeriod<-c(agilent[seq(5770,rawExtractLength,4)][-1],0)
  
  ## extract second, third and fourth between for the SCAN TIME
  cat('...extracting scantimes...\n')
  betweenFirst<-betweenSequence(firstPeriod,counts)
  betweenSecond<-betweenSequence(secondPeriod,counts)
  betweenThird<-betweenSequence(thirdPeriod,counts)
  betweenFourth<-betweenSequence(fourthPeriod,counts)
  scanTime<-betweenFirst[seq(1,8*numberOfScans,8)]*16777216
                          +betweenSecond[seq(1,8*numberOfScans,8)]*65536
                          +betweenThird[seq(1,8*numberOfScans,8)]*256
                          +betweenFourth[seq(1,8*numberOfScans,8)]
  scanTime<-round(scanTime/1000/60,4)
  
  ## extract main sequence, reverse them scan wise
  mainSequence<-function(period,counts){
    tempSequence<-NULL
    currentPosition<-5
    for(ii in counts){
      tempSequence<-c(tempSequence,period[(currentPosition+ii-1):currentPosition])
      currentPosition<-(currentPosition+ii+7)    
    }
    return(tempSequence)
  }
  
  ### extract main data for INT and MZ
  cat('...extract intensity and Mz data...\n')
  mainFirst<-mainSequence(firstPeriod,counts)
  mainSecond<-mainSequence(secondPeriod,counts)
  mainThird<-mainSequence(thirdPeriod,counts)
  mainFourth<-mainSequence(fourthPeriod,counts)
  
  ### calculate MZs. This will result in the *real*, used MZs. For the case that the
  ### detector is switched off, zero values are included.
  importMz<-round(mainFirst*12.8+mainSecond*0.05)
  
  ### calculate intensity values
  cat('...calculate intensity values...\n')
  mainFourth[which(floor(mainThird/64)==3)]<-(mainFourth[which(floor(mainThird/64)==3)]*512)
  mainFourth[which(floor(mainThird/64)==2)]<-(mainFourth[which(floor(mainThird/64)==2)]*64)
  mainFourth[which(floor(mainThird/64)==1)]<-(mainFourth[which(floor(mainThird/64)==1)]*8)
  
  mainThird[which(floor(mainThird/64)==3)]<-((mainThird[which(floor(mainThird/64)==3)] %% 192)*131072)
  mainThird[which(floor(mainThird/64)==2)]<-((mainThird[which(floor(mainThird/64)==2)] %% 128)*16384)
  mainThird[which(floor(mainThird/64)==1)]<-((mainThird[which(floor(mainThird/64)==1)] %% 64)*2048)
  mainThird[which(floor(mainThird/64)==0)]<-(mainThird[which(floor(mainThird/64)==0)]*256)
  
  importInt<-(mainThird+mainFourth)
  
  cat('...assembling data matrix...\n')
  DATA<-matrix(rep(0,numberOfScans*max(importMz)),numberOfScans,max(importMz))
  position<-1
  for(ii in 1:numberOfScans){
    
    for(jj in counts[ii]){
      
      if(counts[ii]){
        DATA[ii,importMz[position:(position+jj-1)]]<-importInt[position:(position+jj-1)]
        position<-position+jj
      }else{
        position<-position+1
      }
    } 
  }
  
  rownames(DATA)<-scanTime
  colnames(DATA)<-seq(1:dim(DATA)[2])
  
  #### Output is a matrix of ion counts with rows as scantime and columns as mass, 
  #### and the respective values as labels 
  return(DATA) 
}
