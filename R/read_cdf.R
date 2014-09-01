##' Function read_cdf
##' 
##' Function read_cdf
##' @export
##' @param projectpath
##' @param filepath
##' @param method
##' @return Xbc, SCAN_INFO, SCAN_RANGE, file
read_cdf <-function(projectpath,filepath, method=1){
  require(ncdf)
  require(tcltk)

  if(missing(filepath))
    cdffiles <- tk_choose.files(caption="Select CDF files to import.",multi=TRUE,filters = matrix(c("CDF files (*.CDF)","*.CDF","all","*.*"),2,2,byrow=TRUE))

  else if(file.exists(filepath))
    cdffiles <- filepath

  if(length(cdffiles)){
    cat("Selected files: \n",paste(cdffiles,"\n ",collapse=""))
    }else{
      cat("No CDF files selected.")
      return(NULL)
      }
  
  dir.create(file.path(projectpath,"Filtered","CDF"),showWarnings=FALSE,recursive=TRUE)

  SETTINGS  <-  NULL

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

  if(file.exists(file.path(projectpath,"maxMZ.Rdata")))
    load(file.path(projectpath,"maxMZ.Rdata"))
  else
    maxMZ <- numeric()

  filepaths <-  file.path(projectpath,"Filtered","CDF",paste(sub("(.+)[.][^.]+$", "\\1", basename(cdffiles)),".Rdata",sep=""))

  cat("================================\n")
  for(i in 1:length(cdffiles)){
    errorcounter  <- 0
    cdffile       <- open.ncdf(cdffiles[i])

    cat("File ",basename(cdffiles[i]), "opened. ",paste("(",i,"/",length(cdffiles),")",sep=""),"\n")
    MZ    <- get.var.ncdf(cdffile,varid=cdffile$var[["mass_values"]])
    cat("Reading variables...\n")
    INT   <- get.var.ncdf(cdffile,varid=cdffile$var[["intensity_values"]])
    TIME  <- get.var.ncdf(cdffile,varid=cdffile$var[["scan_acquisition_time"]])
    count <- get.var.ncdf(cdffile,varid=cdffile$var[["point_count"]])
    MZmin <- min(get.var.ncdf(cdffile,varid=cdffile$var[["mass_range_min"]]))
    MZmax <- max(get.var.ncdf(cdffile,varid=cdffile$var[["mass_range_max"]]))

    close.ncdf(cdffile)

    cat("File ",basename(cdffiles[i]), "closed.\n")
    if(method==1)
      DATA <- do_DATA1(INT,MZ,TIME,count,MZmin,MZmax,MZP)
    if(method==2)
      DATA <- do_DATA2(INT,MZ,TIME,count,MZmin,MZmax,MZP)
    if(method==3)
      DATA <- do_DATA3(INT,MZ,TIME,count,MZmin,MZmax,MZP)

    rm(INT,MZ,count)

    cat(paste("Smoothing (FL = ",FL,")\n",sep=""))

    Xbc <- try(baseline(DATA,projectpath,FL),silent=TRUE)

    while(mode(Xbc) == "character"){
      errorcounter  <-  errorcounter + 1
      cat(Xbc[1])
      cat("Garbage collection<U+2026> retrying..\n")
      gc()
      Xbc <- try(baseline(DATA,projectpath,FL),silent=TRUE)
      if(errorcounter > 4 & mode(Xbc) == "character")
        stop("Error!")
      }

    rm(DATA)

    mz        <- SCAN_RANGE <- seq(max(1,round(MZmin)),min(ncol(Xbc),round(MZmax),MZR[2]))
    Xbc       <- Xbc[,mz]
    SCAN_INFO <- matrix(cbind(1:length(TIME),TIME),nrow=length(TIME),ncol=2)
    maxMZ     <- max(maxMZ,ncol(Xbc))

    save(Xbc,SCAN_INFO,SCAN_RANGE,file = filepaths[i])
    rm(TIME)
    

    if(missing(filepath))
      rm(Xbc,SCAN_INFO,SCAN_RANGE)

    cat("================================\n")

    gc()

    }

  save(maxMZ,file=file.path(projectpath,"maxMZ.Rdata"))
  if(missing(filepath)){
    cat("Done! \nAll files imported and stored in ",file.path(projectpath,"Filtered","CDF"),"\n\n")
    gc()           # Final garbage collection
    }else
      import = list(Xbc=Xbc,SCAN_INFO=SCAN_INFO,SCAN_RANGE=SCAN_RANGE,file=basename(filepath))	  # Return variables Xbc, SCAN_INFO and SCAN_RANGE if a valid file path was supplied
}

