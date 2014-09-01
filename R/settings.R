##' Basic Settings for GCMS Data Processing 
##'
##' \code{settings} creates the file \code{settings.Rdata} in the current
##' workdirectory with default settings for the \code{GCMS} data processing
##' script. If the file already exists a console based user interface to
##' modify the values. 
##' @export
##' @title GUI for general GCMS Data processing settings
##' @param projectpath Projectpath where the directory structure is built up in 
##' @param setdefault If \code{TRUE}, all values are set to default. 
##' @return none
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
settings <-function(projectpath,setdefault = FALSE){
  if(setdefault){
    SETTINGS     <-    list()
    SETTINGS$FL  <-    3
    SETTINGS$MS  <-    500
    SETTINGS$MZP <-    0
	  	SETTINGS$MZR  <- 	c(35,250)
	  	SETTINGS$BC1 	<-  TRUE
	  	SETTINGS$R2Xc	<-	0.95
	  	SETTINGS$BC2  <-	TRUE
	  	SETTINGS$MPS  <-	45
	  	SETTINGS$NL   <-	3
	  	SETTINGS$RP   <-	1
	  	save(SETTINGS,file=file.path(projectpath,"SETTINGS.Rdata"))

	}else if(file.exists(file.path(projectpath,"SETTINGS.Rdata"))){
		load(file.path(projectpath,"SETTINGS.Rdata"))
		cat("===================================\n")
		cat("== GCMS Data Processing Settings ==\n")
		cat("===================================\n")

		j <-  1
		
		while(j){
			j <-  menu(c("General","MCR","Done"),title="Settings")
			
			if(j == 1){
				k <-  1
				
				while(k){
					k	<-	menu(c(paste("Filterlength:\t",SETTINGS$FL),paste("Max shift:\t\t",SETTINGS$MS),paste("MZ precision:\t",SETTINGS$MZP),paste("MZ minimum:\t\t",SETTINGS$MZR[1]),paste("MZ maximum:\t\t",SETTINGS$MZR[2]),"Default","Done"),title="General settings:")
					if(k == 1){
						SETTINGS$FL  <-  max(as.numeric(menu(as.character(c(0,1,2,3,4,5)),graphics=TRUE))-1,0)
					
					}else if(k == 2){
						
						while(is.na(SETTINGS$MS  <- as.numeric(readline("Max shift: ")))){}
					
					}else if(k == 3){
						SETTINGS$MZP  <-  max(as.numeric(menu(as.character(c(0,1,2)),graphics=TRUE))-1,0)
						
					}else if(k == 4){
						
						while(is.na(SETTINGS$MZR[1]  <- as.numeric(readline("MZ minimum: ")))){}
						
					}else if(k == 5){
						
						while(is.na(SETTINGS$MZR[2]  <- as.numeric(readline("MZ maximum: ")))){}
						
					}else if(k == 6){
						SETTINGS$FL		<-	3
						SETTINGS$MS		<-	500
						SETTINGS$MCP	<-	0
	  	    
	  	    		}else if(k == 7)
	  	    			k <-  0
	  	    			
	  	    	}
			
			}else if(j == 2){
				k <-  1

				while(k){
					k	<-	menu(c(paste("Baseline correction:\t",SETTINGS$BC2),paste("Max peak shift:\t",SETTINGS$MPS),paste("Noise limit:\t\t",SETTINGS$NL),paste("Retention precision:\t",SETTINGS$RP),"Default settings","Done"),title="MCR settings:")

					if(k == 1){
						SETTINGS$BC2  <-  as.logical(max(menu(as.character(c(FALSE,TRUE)),graphics=TRUE)-1,0))
					
					}else if(k == 2){
						
						while(is.na(SETTINGS$MPS  <- as.numeric(readline("Max peak shift: ")))){}
					
					}else if(k == 3){
						
						while(is.na(SETTINGS$NL  <- as.numeric(readline("Noise limit: ")))){}
					
					}else if(k == 4){
						
						while(is.na(SETTINGS$RP  <- as.numeric(readline("Retention precision: ")))){}
					
					}else if(k == 5){
						SETTINGS$BC2	<-	TRUE
						SETTINGS$MPS	<-	45
						SETTINGS$NL		<-	3
						SETTINGS$RP		<-	1
						
					}else if(k == 6)
		  				k	<-  0

				}
			
			}else if(j == 3)
				j <-  0
		}
	  	
	  	save(SETTINGS,file=file.path(projectpath,"SETTINGS.Rdata"))
	}else
    
    settings(projectpath,TRUE)
}

