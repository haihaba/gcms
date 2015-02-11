scoresAR	<-	function(predpath,DATA){
  
	if(missing(DATA))
 	  load(file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))

 	load(file.path(predpath,"info.Rdata"))
	load(file.path(predpath,"Aligned","SCAN_RANGE.Rdata"))
	D					<-	numeric()
	VARID1		<-	character()
 	spectrum	<- numeric()
	for(i in 1:max(ROWID1[,1]))
	{
    textstr	<-	paste("Checking window number:", i)
		cat(textstr,"\n")
    #do_log(predpath,textstr)
    int   <-   which(ROWID1[,1]==i)
    XX		<-	DATA[,int]    # DATA = X
    SSX		<-	sum(XX^2)

    if(min(dim(XX))>2)
    {
    	CS					<-	ARmethod1(XX,predpath)
    	comp				<-	ncol(CS$C)
    	D						<-	cbind(D,CS$C)
    	S_temp			<-	matrix(0,nrow=comp,ncol=max(SCAN_RANGE))
    	MZ					<-	ROWID1[int,2]
    	S_temp[,MZ]	<-	t(CS$S)
    	spectrum		<-	rbind(spectrum,S_temp)
    	for(j in 1:comp)
    	{
    		if(i > 99)
    		{
      		if(j < 10)
					{
        		VARID1	<-	rbind(VARID1,paste("W",i,"_C0",j))
        	}else
        		VARID1	<-	rbind(VARID1,paste("W",i,"_C",j))
      	}else
      	{
					if(i < 10)
      		{
      			if(j < 10)
      			{
        			VARID1	<-	rbind(VARID1,paste("W00",i,"_C0",j))
        		}else
        			VARID1	<-	rbind(VARID1,paste("W00",i,"_C",j))
      		}else
						if(j < 10)
      			{
      				VARID1	<-	rbind(VARID1,paste("W0",i,"_C0",j))
      			}else
      				VARID1	<-	rbind(VARID1,paste("W0",i,"_C",j))
				 }
		  }
		}
	}
  #output	<-	list(D=D,VARID1=VARID1,spectrum=spectrum)
  DATA  <-  D
  load(file.path(predpath,"Aligned","COLUMNID1.Rdata"))
  OBSID1	<-	COLUMNID1
  save(DATA,OBSID1,VARID1,file=file.path(predpath,"HMCR","REG","MVA_DATA.Rdata"))
  write.table(data.frame(ID=OBSID1,DATA),file=file.path(predpath,"HMCR","REG","MVA_DATA.txt"),row.names=FALSE,col.names=c("PrimaryID",VARID1),sep="\t",quote=FALSE)
	save(spectrum,VARID1,file=file.path(predpath,"HMCR","REG","SPECTRUM.Rdata"))


  load(files[which.min(shift)])
	EDGES_TIME	<-	SCAN_INFO[,2]
  save(EDGES_TIME,file=file.path(predpath,"EDGES_TIME.Rdata"))

  #do_log(predpath,"Ending method 1")
  #do_log(predpath,"=================")
  
  # Not done
  #if(file.exists.file.path(predpath,"Alkane serie","RT_INFO.Rdata"))
  #{
	#	load(file.path(predpath,"Alkane_serie","RT_INFO.Rdata"))
	#	show_spec2(predpath)
	#}
}
