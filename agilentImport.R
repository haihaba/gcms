readDFile<-function(pathname){
  
  filename<-file.path(pathname,'DATA.MS')
  
  cat('Opening Agilent file...\n')
  to.read<-file(filename,'rb')
  agilent<-readBin(to.read,integer(),size=1,signed=FALSE,n=20000000)
  close(to.read)
  
  ### preparing vector with counts per scan
  cat('...extracting counts per scan...\n')
  countsNumber<-agilent[5782]
  counts<-agilent[5782]
  currentPosition<-5782
  
  while(currentPosition<length(agilent)){
    jumpLength<-countsNumber*4+7*4
    currentPosition<-currentPosition+jumpLength
    countsNumber<-agilent[currentPosition]
    counts<-c(counts,agilent[currentPosition])
  }
  
  
  ### counts will be too long and the last entry is NA
  ### counts will be corrected when the number of scans
  ### is known.
  #counts<-countNumber(agilent)
  ### cut away NA. 
  counts<-counts[-which(is.na(counts))]
  
  ### the second period is extracted. This one is currently
  ### used to determine the number of scans. As the useful length
  ### is not known yet, it is extracted in overlength. Na's have
  ### to be removed then
  cat('...determine number of scans...\n')
  secondPeriod<-agilent[seq(5772,2000000,4)]
  ### remove Na's
  secondPeriod<-secondPeriod[-which(is.na(secondPeriod))]
  
  
  ### inline function to extract the data paddings
  betweenSequence<-function(period,counts){
    tempSequence<-period[1:4]
    currentPosition<-4
    for(ii in counts){
      currentPosition<-(currentPosition+ii)
      tempSequence<-c(tempSequence,period[currentPosition:(currentPosition+7)])
      currentPosition<-(currentPosition+7)
    }
    return(tempSequence)  
  }
  
  ### the betweenSequence removes the variable length ion count
  ### data. Left are 8 numbers of padding betweeen each scan
  betweenSecond<-betweenSequence(secondPeriod,counts)
  ### in the betweenSequence of the second Period, the 8th
  ### field is currently used for determination of scan numbers
  ### when it is 3 times in sequence zero,
  numberOfScans<-which.max(abs(diff(betweenSecond[seq(1,100000,8)]))>1)
  counts<-counts[1:numberOfScans]
  rawExtractLength<-sum(counts)*4+(numberOfScans*4*7)+5770
  
  ### extract periods with correct length
  cat('...separate rawdata in four sequences...\n')
  firstPeriod<-agilent[seq(5771,rawExtractLength,4)]
  secondPeriod<-agilent[seq(5772,rawExtractLength,4)]
  thirdPeriod<-agilent[seq(5773,rawExtractLength,4)]
  fourthPeriod<-c(agilent[seq(5770,rawExtractLength,4)][-1],0)
  
  ## extract third and fourth between for the scan time
  cat('...extracting scantimes...\n')
  betweenSecond<-betweenSequence(secondPeriod,counts)
  betweenThird<-betweenSequence(thirdPeriod,counts)
  betweenFourth<-betweenSequence(fourthPeriod,counts)
  scanTime<-betweenSecond[seq(1,8*numberOfScans,8)]*65536+betweenThird[seq(1,8*numberOfScans,8)]*256+betweenFourth[seq(1,8*numberOfScans,8)]
  scanTime<-round(scanTime/1000/60,4)
  
  ## extract main sequence, reverse them scan wise
  mainSequence<-function(period,counts){
    tempSequence<-NULL
    currentPosition<-5
    for(ii in counts){
      tempSequence<-c(tempSequence,period[(currentPosition+ii-1):currentPosition])
      currentPosition<-(currentPosition+ii+7)    
    }
    return(tempSequence)
  }
  
  ### extract main data
  cat('...extract intensity and Mz data...\n')
  mainFirst<-mainSequence(firstPeriod,counts)
  mainSecond<-mainSequence(secondPeriod,counts)
  mainThird<-mainSequence(thirdPeriod,counts)
  mainFourth<-mainSequence(fourthPeriod,counts)
  
  ### calculate MZs
  importMz<-round(mainFirst*12.8+mainSecond*0.05)
  
  ### calculate intensity values
  cat('...calculate intensity values...\n')
  mainFourth[which(floor(mainThird/64)==3)]<-(mainFourth[which(floor(mainThird/64)==3)]*512)
  mainFourth[which(floor(mainThird/64)==2)]<-(mainFourth[which(floor(mainThird/64)==2)]*64)
  mainFourth[which(floor(mainThird/64)==1)]<-(mainFourth[which(floor(mainThird/64)==1)]*8)
  
  mainThird[which(floor(mainThird/64)==3)]<-((mainThird[which(floor(mainThird/64)==3)] %% 192)*131072)
  mainThird[which(floor(mainThird/64)==2)]<-((mainThird[which(floor(mainThird/64)==2)] %% 128)*16384)
  mainThird[which(floor(mainThird/64)==1)]<-((mainThird[which(floor(mainThird/64)==1)] %% 64)*2048)
  mainThird[which(floor(mainThird/64)==0)]<-(mainThird[which(floor(mainThird/64)==0)]*256)
  
  importInt<-(mainThird+mainFourth)
  ### assign data to a Scan by Mz Matrix
  cat('...assemble full data matrix...\n')
  fullData<-matrix(rep(0,numberOfScans*max(importMz)),numberOfScans,max(importMz))
  position<-1
  for(ii in 1:numberOfScans){
    
    for(jj in counts[ii]){
      
      if(counts[ii]){
        fullData[ii,importMz[position:(position+jj-1)]]<-importInt[position:(position+jj-1)]
        position<-position+jj
      }else{
        position<-position+1
      }
    } 
  }
  rownames(fullData)<-scanTime
  return(fullData)
  
}