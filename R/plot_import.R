##' Function plot_import
##' 
##' Function plot_import
##' @param projectpath
plot_import<-function(projectpath){
  require(tcltk)

  if(file.exists(file.path(projectpath,'Filtered')))
    setwd(file.path(projectpath,'Filtered'))

  rdatafiles <- tk_choose.files(caption="Select imported Rdata files", multi=TRUE)
  setwd(projectpath)
     
  load(rdatafiles[1])
  plot(apply(Xbc,1,sum),type='l', ylab='signal',xlab='scans')
  rm(list=c('Xbc','SCAN_INFO','SCAN_RANGE'))

  if(length(rdatafiles)>1){

    for(i in 2:length(rdatafiles)){
      load(rdatafiles[i])
      points(apply(Xbc,1,sum),type='l')
      rm(list=c('Xbc','SCAN_INFO','SCAN_RANGE'))
    }
  }
}
