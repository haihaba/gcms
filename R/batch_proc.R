##' Function batch_proc
##' 
##' Function batch_proc
##' @export
##' @param projectpath
##' @param cores
batch_proc <- function( projectpath,cores){
  require(parallel)
  load(file.path(projectpath,"HMCR","MCR.Rdata"))
  winlist<-grep("[P]",MCR$windowlist)
  incl<-MCR$reg$incl
  pred<-MCR$reg$pred
  win<-unlist(mclapply(winlist, function(i) find_spectrum(projectpath,incl,pred,i)))
  read_win(getwd(),'REG',win)

}
