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
  
  

  
  
  
  ## construct the filepaths for the files to save
  #filepaths <-  file.path(projectpath,"Filtered","CDF",paste(sub("(.+)[.][^.]+$", "\\1", basename(cdffiles)),".Rdata",sep=""))
  #cat("================================\n")
  
  
  
  
  
  ## looping through all the CDF files for importing
  for(i in 1:length(cdffiles)){
    errorcounter  <- 0
    
    
    ## ncdf files from agilent
    ## variabels
    ## error_log: 63 chars, not useful
    ## a_d_sampling_rate, -9999, not useful
    ## a_d_coaddition_factor, -9999 not useful
    ## scan_acquisition_time, decimal seconds useful
    ## scan_duration, -9999, not useful
    ## inter_scan_time, 
    
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