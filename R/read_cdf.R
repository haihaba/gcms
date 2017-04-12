##' Function read_cdf
##' 
##' Function read_cdf
##' @export
##' @param projectpath
##' @param filepath     can be used to import a single file and return the data to the caller
##' @return Xbc, SCAN_INFO, SCAN_RANGE, file
read_cdf <-function(projectpath,filepath){
  require(ncdf4)
  require(tcltk)
  
  
  ### file selection by TclTk
  if(missing(filepath))
    cdffiles <- tk_choose.files(caption="Select CDF files to import.",multi=TRUE,filters = matrix(c("CDF files (*.CDF)","*.CDF","all","*.*"),2,2,byrow=TRUE))
  else if(file.exists(filepath))
    cdffiles <- filepath
  
  if(length(cdffiles)){
    cat("Selected files: \n",paste(cdffiles,"\n ",collapse=""))
  } else {
    cat("No CDF files selected.")
    return(NULL)
  }
  
  
  
  ## create directory for CDF raw data import
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
    errorcounter  <- 0
    
    
    ## ncdf file import
    cdffile       <- nc_open(cdffiles[i])
    cat("File ",basename(cdffiles[i]), "opened. ",paste("(",i,"/",length(cdffiles),")",sep=""),"\n")
    MZ    <- ncvar_get(cdffile,varid=cdffile$var[["mass_values"]])
    cat("Reading variables...\n")
    INT   <- ncvar_get(cdffile,varid=cdffile$var[["intensity_values"]])
    TIME  <- ncvar_get(cdffile,varid=cdffile$var[["scan_acquisition_time"]])
    count <- ncvar_get(cdffile,varid=cdffile$var[["point_count"]])
    MZmin <- min(ncvar_get(cdffile,varid=cdffile$var[["mass_range_min"]]))
    MZmax <- max(ncvar_get(cdffile,varid=cdffile$var[["mass_range_max"]]))

    nc_close(cdffile)
    cat("File ",basename(cdffiles[i]), "closed.\n")

    
    
    ## Data imoprt for Agilent Chemstation CDF files 
    ## For ther CDF's eventually some adaptation is needed.
    DATA <- do_DATA1(INT,MZ,TIME,count,MZmin,MZmax,MZP)
  
    
    
    
    ## the raw CDF import variables are discarded to free memory
    rm(INT,MZ,count)

    
    
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
    
    ## when funciton argument 'filepath' missing, delete the variables
    ## otherwise they are still needed to return to the caller
    if(missing(filepath))
      rm(Xbc,SCAN_INFO,SCAN_RANGE)

    cat("================================\n")

    }

  save(maxMZ,file=file.path(projectpath,"maxMZ.Rdata"))
  
  
  ## return data in case filepath was provided as function argument
  if(missing(filepath)){
    cat("Done! \nAll files imported and stored in ",file.path(projectpath,"Filtered","CDF"),"\n\n")
    }else
      import = list(Xbc=Xbc,SCAN_INFO=SCAN_INFO,SCAN_RANGE=SCAN_RANGE,file=basename(filepath))	  # Return variables Xbc, SCAN_INFO and SCAN_RANGE if a valid file path was supplied
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
    Xbc[start:end,]<-stats:::filter(X,filter=rep(1/(FL*2+1),(FL*2+1)))[start:end,]
    
    
    if(any(Xbc<0))
      Xbc[Xbc<0]	<-	0
  }else
    Xbc <-  X
  return(Xbc)
}

##' Function do_DATA1
##' 
##' Function do_DATA1 is the method used with Agilent Chemstation data 
##' and is set as default. It is called from read_cdf.R
##' This function is used to arrange the initial data table after import
##' @param INT    intensity values from CDF import
##' @param MZ     mass values from CDF import
##' @param TIME   scantime values from CDF import
##' @param count  point count from CDF import
##' @param MZmin  min value mass range from CDF import
##' @param MZmax  max value mass range from CDF import
##' @param MZP    mass precision from settings file
##' @return DATA  rawdata matrix with rows for timepoints
##' and columns for mz values. The mz values range from 
##' zero to maxMZ.
do_DATA1<-function(INT,MZ,TIME,count,MZmin,MZmax,MZP){
  
  ## prepare/process rawdata
  MZ        <- round(MZ,digits = MZP)*10^MZP
  INT       <- round(INT)
  lenINT    <- length(INT)
  lencount  <- length(count)
  lenTIME   <- length(TIME)
  
  ## check if all mass values from min to max are abundant. This is usually not 
  ## the case in Agilent datasets. Therefore the dataset has to be constructed
  ## by looping through the individual datavalues. Otherwise, a matrix assignment
  ## can be done directly
  if(length(abs(diff(MZ))[abs(diff(MZ)) != 1 & abs(diff(MZ)) != (MZmax-MZmin)])){ #, !all(count == (MZmax-MZmin+1))))
    cat("Using alternative method to arrange data matrix (slower)\n")
    
    DATA  <-  matrix(0,lenTIME,max(MZ))
    scans <-  min(which(count != 0))
    points <- count[min(which(count != 0))]
    
    for(i in 1:lenINT){
      points <- points-1
      DATA[scans,MZ[i]]  <- DATA[scans,MZ[i]] + INT[i]
      if(i != lenINT)
        if(points == 0){
          scans  <- scans + min(which(count[(scans+1):lencount] != 0))
          points <- count[scans]
        }
    }
  }else
    DATA <- cbind(matrix(0,nrow=lenTIME,ncol=max((min(MZ)-1),0)),matrix(INT,nrow=lenTIME,ncol=(max(MZ)-min(MZ)+1),byrow=TRUE))
  return(DATA)
}


