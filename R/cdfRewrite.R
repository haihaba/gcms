##' Function cdfRewrite
##' 
##' Function cdfRewrite
##' @export
##' @param projectpath
##' @param filepath can be used to rewrite files in batchmode
cdfRewrite <-function(projectpath,filepath){
  require(ncdf)
  require(tcltk)
  
  ### file selection by TclTk
  if(missing(filepath))
    cdffiles <- tk_choose.files(caption="Select CDF files to rewrite.",multi=TRUE,filters = matrix(c("CDF files (*.CDF)","*.CDF","all","*.*"),2,2,byrow=TRUE))
  else if(file.exists(filepath))
    cdffiles <- filepath
  
  if(length(cdffiles)){
    cat("Selected files: \n",paste(cdffiles,"\n ",collapse=""))
  } else {
    cat("No CDF files selected.")
    return(NULL)
  }
  
  
  ## create directory for CDF raw data import
  ##dir.create(file.path(projectpath,"Filtered","CDF"),showWarnings=FALSE,recursive=TRUE)
  
  
  
  
  
  ## Settings Handling
  #SETTINGS  <-  NULL
  
  
  ## Check if a setting file exists, if yes read 
  ## it else create by setting default values
  #if(file.exists(file.path(projectpath,"SETTINGS.Rdata"))){
  #  load(file.path(projectpath,"SETTINGS.Rdata"))
  #  FL   <- SETTINGS$FL
  #  MZP  <- SETTINGS$MZP
  #  MZR  <- SETTINGS$MZR
  #}else{
  #  cat("No filterlength settings found, using FL = 0.\n")
  #  FL  <-  0
  #  cat("No mass channel precision settings found, using MZP = 0.\n")
  #  MZP <-  0
  #  cat("No mass channel range settings found, using MZR min = 50, max = 800.\n")
  #  MZR <-  c(50,800)
  #}
  
  
  ## Check if maximum MZ Rdata file exist, then 
  ## load, otherwise create the variable maxMZ
  #if(file.exists(file.path(projectpath,"maxMZ.Rdata")))
  #  load(file.path(projectpath,"maxMZ.Rdata"))
  #else
  #  maxMZ <- numeric()
  
  
  
  ## construct the filepaths for the files to save
  #filepaths <-  file.path(projectpath,"Filtered","CDF",paste(sub("(.+)[.][^.]+$", "\\1", basename(cdffiles)),".Rdata",sep=""))
  #cat("================================\n")
  
  
  
  
  
  ## looping through all the CDF files for importing
  for(i in 1:length(cdffiles)){
    errorcounter  <- 0
    
    
    ## ncdf file import
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
  }
}