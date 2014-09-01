##' Parallel by mclapply
##' 
##' @export
parallelMCR<-function(projectpath){
  
  load(file.path(projectpath,"HMCR","MCR.Rdata"))
  samples<-MCR$reg$samples
  incl<-MCR$reg$incl
  excl<-MCR$reg$excl
  pred<-MCR$reg$pred
  MCR$reg$excluded<-excl
  win<-find_spectrum(projectpath,incl,pred,wincores)


}