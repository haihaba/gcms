##' Function submitbatch
##' 
##' Function submitbatch
##' @export
##' @param projectpath
##' @param cores
submitbatch<-function(projectpath,cores){
  if(!file.exists(file.path(projectpath,"batch.R"))){
    cat("library(batch)\n",file="batch.R")
    cat("projectpath<-getwd()\n",file="batch.R",append=T)
    cat("library(GCMS)\n",file="batch.R",append=T)
    cat("load(file.path(projectpath,\"HMCR\",\"MCR.Rdata\"))\n",file="batch.R",append=T)
    cat("wincores<-0\n",file="batch.R",append=T)
    cat("parseCommandArgs()\n",file="batch.R",append=T)
    cat("samples<-MCR$reg$samples\n",file="batch.R",append=T)
    cat("incl<-MCR$reg$incl\n",file="batch.R",append=T)
    cat("excl<-MCR$reg$excl\n",file="batch.R",append=T)
    cat("pred<-MCR$reg$pred\n",file="batch.R",append=T)
    cat("MCR$reg$excluded<-excl\n",file="batch.R",append=T)
    cat("win<-find_spectrum(projectpath,incl,pred,wincores)\n",file="batch.R",append=T)
  }

  
  load(file.path(projectpath,"HMCR","MCR.Rdata"))
  winlist<-as.numeric(substr(MCR$windowlist[grep("[P]",MCR$windowlist)],8,10))
  if(length(winlist)<cores)
  	cores<-length(winlist)

  batches<-list()
  corecounter<-1

  for(i in winlist){
    eval(parse(text=paste("batches$core",corecounter,"<-c(batches$core",corecounter,",i)",sep="")))
    corecounter<-corecounter+1

    if(corecounter>cores)
      corecounter<-1
  }

  seed<-1000

  for (i in 1:cores){
    eval(parse(text=paste("k<-batches$core",i,sep="")))

    a<-as.character()

    for (j in 1:length(k)){
      a<-paste(a,",",k[j],sep="")
    }
    

    system(paste("R --vanilla --args wincores \"c(",substr(a,2,nchar(a)),")\" < \"batch.R\"",sep=""),wait=F,ignore.stdout=T)
  }
}

