## trying to figure out the binary file of the Agilent Chemstation
## using the original binary file and a CDF export from chemstation

### loading cdf file
library(ncdf)
cdffile<-open.ncdf('~/Google Drive/Programming/R//importAgilent//0107.CDF')
cdffileMz    <- get.var.ncdf(cdffile,varid=cdffile$var[["mass_values"]])
cdffileInt<- get.var.ncdf(cdffile,varid=cdffile$var[["intensity_values"]])
cdffileTime  <- get.var.ncdf(cdffile,varid=cdffile$var[["scan_acquisition_time"]])
cdffileCount <- get.var.ncdf(cdffile,varid=cdffile$var[["point_count"]])
cdffileMzMin <- min(get.var.ncdf(cdffile,varid=cdffile$var[["mass_range_min"]]))
cdffileMzMax <- max(get.var.ncdf(cdffile,varid=cdffile$var[["mass_range_max"]]))
close.ncdf(cdffile)


### loading agilent file
to.read<-file('~/Google Drive/Programming/R//importAgilent//0107.D//DATA.MS','rb')
agilent.big<-readBin(to.read,integer(),size=1,signed=FALSE,n=20000000, endian = 'big')
close(to.read)

agilent<-readBin(to.read,character(), n=7000)

### plot every thousandst point from the agilent file
### looks like the end could be time values
plot(agilent[seq(7,2000000,1000)],type='l')

### there seems to be a period of 4 in the data when read as single bytes
plot(agilent[seq(150003,169000,2)],type='l')

### 2 of the 4 series are clearly periodic
par(mfcol=c(2,2))
plot(agilent[seq(150001,169000,4)],type='l')
plot(agilent[seq(150002,169000,4)],type='l')
plot(agilent[seq(150003,169000,4)],type='l')
plot(agilent[seq(150004,169000,4)],type='l')

par(mfcol=c(4,1))
plot(agilent[seq(5773,25000,4)],type='l')
plot(agilent[seq(5774,25000,4)],type='l')
plot(agilent[seq(5775,25000,4)],type='l')
plot(agilent[seq(5776,25000,4)],type='l')

### stepping through the whole dataset
par(mfcol=c(5,1))
for(i in 0:72){
  plot(agilent[seq((5773+(i*25000)),(25000+(i*25000)),4)],type='l')
  par(ask=FALSE)
  plot(agilent[seq((5774+(i*25000)),(25000+(i*25000)),4)],type='l')
  plot(agilent[seq((5775+(i*25000)),(25000+(i*25000)),4)],type='l')
  plot(agilent[seq((5776+(i*25000)),(25000+(i*25000)),4)],type='l')
  par(ask=TRUE)
}

### the 5775 starting sequence could be some sort of counts or TIC but what for???



### this one's nice... looks like there is just one type 
### of data throughout the whole file
plot(agilent[seq(5001,2015000,4)],type='l')

plot(agilent[seq(5001,200000,4)],type='l')

### oops, plotting all series, seems like there are two different
### kinds of data. In the end shows up another type of data.
### Around 1´680´000 points.
### seems like in the series starting at 5771, the value 16 forms the begin of 
### a new scan or half scan. There are around 14´000 such peaks with the value 16,
### double the number of scans. 

### what is needed? Int values and mass values. There should be around 368´000
### Checked cdffileCount against the series of value 16. They clearly correlate
par(mfcol=c(2,1))
plot(diff(which(agilent==16))[which(diff(which(agilent==16))>4)][100:11000],type='l')
plot(cdffileCount,type='l')

### converter from decimal to binary
intToBin <- function(x, ndigits=0, b=2){
    xi <- as.integer(x)
    if(any(is.na(xi) | ((x-xi)!=0)))
      print(list(ERROR="x not integer", x=x))
    N <- length(x)
    xMax <- max(x)
    if(!ndigits)
      ndigits <- (floor(logb(xMax, base=2))+1)
    Base.b <- array(NA, dim=c(N, ndigits))
    for(i in 1:ndigits){#i <- 1
      Base.b[, ndigits-i+1] <- (x %% b)
      x <- (x %/% b)
    }
    if(N ==1) Base.b[1, ] else Base.b
  }

### convert binary to decimal
binToInt <- function(x){
  x<-rev(x)
  if(any(x>1))
    stop(cat('x is not binary'))
  multiplicator<-1
  for(i in 1:(length(x)-1))
    multiplicator<-c(multiplicator,2^i)
  decimalResult<-sum(decimal<-x*multiplicator)
  return(decimalResult)
}

## looking at the first scan
plot(agilent[seq(5773,5904,4)],type='l')
## alternatively these are two full scans with one '21' value inbetween
plot(agilent[seq(5773,6064,4)],type='l')

## first scan
par(mfcol=c(2,1))
plot(cdffileInt[1:34],type='l')
plot(agilent[seq(5770,5901,4)],type='l')
plot(agilent[seq(5771,5902,4)],type='l')
plot(
  ,type='l')
plot(rev(agilent[seq(5773,5904,4)]),type='l')

## Check in the 5775 series how many 16 tops there are
## result is 6075 and as such not too far away from the number of 
## scans in the CDF file.
length(which(diff(which(agilent[seq(5775,2000000,4)]==16))>1))

agilent[5773:5907]

## looks like the scan is from high to low mz!!!!!
## the 5773 sequence is to be multiplied by 256 and the value
## of the 5774 sequence is to be added for the INT value

rev(agilent[seq(5770,5901,4)])
rev(agilent[seq(5771,5902,4)])
rev(agilent[seq(5772,5903,4)])
rev(agilent[seq(5773,5904,4)])

## 5770 or 5774??? Seems to be 5774. There is a value (66) when starting
## on 5770. But that could also be a general initiater value. Otherwise,
## the series seems to be off by one compared to 5771 and 5772
rev(agilent[seq(5774,5905,4)])

## So, the 5773 and 5774 series are the intensity. What is ..71+72?
## Should be either M/Z or TIME
## seems to be time 3x+250=4x+75 => x=325

rev(agilent[seq(5771,5902,4)])*325+rev(agilent[seq(5772,5903,4)])

## 5771 seriescontains a multiplicator used for either time or M/Z (within one scan)
## but there seems to be another information, a peak at the end of each scan. Is somehow
## correlated with TIC maybe. Two zero values before this peak and three zero
## values after the peak

## general sequence is [4 digits][#counts digits][3 digits]. 


## in the 5774 sequence, the 3rd position signifies the number of
## count values per scan

## make an iterator that extracts the number of counts per scan

countNumber<-function(import){
  countsNumber<-import[5782]
  tempSequence<-import[5782]
  currentPosition<-5782
  
  while(currentPosition<length(import)){
    jumpLength<-countsNumber*4+7*4
    currentPosition<-currentPosition+jumpLength
    countsNumber<-import[currentPosition]
    tempSequence<-c(tempSequence,import[currentPosition])                  
  }
  return(tempSequence)
}

### now we have the counts, and we have the int. Missing is the detailed mz
### altough how about 

### funciton below helps to find a certain int value. This is now needed to 
### check how int values above 256*256 are handled.
summedIntPosition <- function(countNumber){
  
  tempSequence<-NULL
  
  for(ii in seq(along=countNumber)){
    tempSequence<-c(tempSequence,sum(countNumber[1:ii]))
  }
  
  return(tempSequence)
}


## 23 February Monday Morning
## make four variables for the four periods of data

## how to determine the starting point of the actual data from the header?
## ??? Will have to check several files

## first period contains a multiplier for the M/Z values. The data is 
## reversed for every scan. The data contains another value that somehow
## correlates with the 'signal amount'. The rest of the measurement is cleary
## the TIC with some multiplier to be applied
firstPeriod<-agilent[seq(5771,1645323,4)]

## extract the addiontal data sequence in the first period using the countNumber function
## hmm, this really could be the TIC
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


## overshoot
secondPeriod<-agilent[seq(5772,2000000,4)]

### until all real scans
secondPeriod<-agilent[seq(5772,1645324,4)]

### including the zero values. Includes seven values per scan
secondPeriod<-agilent[seq(5772,1675980,4)]

## in the first period, the number correlating to the TIC needs to be multiplied
## by 256 in the first few values to come to the closest of the real TIC value. 


## third Period
thirdPeriod<-agilent[seq(5773,1645325,4)]

## according to the the cdffileMZ, mainSecond and mainFirst, the MZ composition
## works as following: Mz 1 = 20 mainSecond digits (or 1 mainSecond digit is 0.05).
## Which makes one mainFirst digit to be 12.8
agilentMz<-mainFirst*12.8+mainSecond*0.05



mainFirst<-mainSequence(firstPeriod,counts)
mainSecond<-mainSequence(secondPeriod,counts)
mainThird<-mainSequence(thirdPeriod,counts)
mainFourth<-mainSequence(fourthPeriod,counts)


## multiplicators that work for calculation of the Int values:
## 1-64: multiplication of the mainThird by 256 and adding main Fourth
## and that can go to something above 16´640

## 65-128: multiplication of the modulus from division by 64 with 
## mainThird by 2048 and adding multiplication of mainFourth by 8

## 129-192: multiplication of the modulus from division by 128 with
## mainThird by 16'384 and adding multiplication of mainFourth by 64

## 193-255: multiplication of the modulus from division by 192 with
## mainThird by 131´072 and adding multiplication of mainFourth by 512


mainFourth[which(floor(mainThird/64)==3)]<-(mainFourth[which(floor(mainThird/64)==3)]*512)
mainFourth[which(floor(mainThird/64)==2)]<-(mainFourth[which(floor(mainThird/64)==2)]*64)
mainFourth[which(floor(mainThird/64)==1)]<-(mainFourth[which(floor(mainThird/64)==1)]*8)

mainThird[which(floor(mainThird/64)==3)]<-((mainThird[which(floor(mainThird/64)==3)] %% 192)*131072)
mainThird[which(floor(mainThird/64)==2)]<-((mainThird[which(floor(mainThird/64)==2)] %% 128)*16384)
mainThird[which(floor(mainThird/64)==1)]<-((mainThird[which(floor(mainThird/64)==1)] %% 64)*2048)
mainThird[which(floor(mainThird/64)==0)]<-(mainThird[which(floor(mainThird/64)==0)]*256)

importInt<-(mainThird+mainFourth)



readDFile<-function(filename){
  to.read<-file(filename,'rb')
  agilent<-readBin(to.read,integer(),size=1,signed=FALSE,n=20000000)
  close(to.read)
  
  ### preparing vector with counts per scan
  countsNumber<-agilent[5782]
  counts<-agilent[5782]
  currentPosition<-5782
    
  while(currentPosition<length(agilent)){
    jumpLength<-countsNumber*4+7*4
    currentPosition<-currentPosition+jumpLength
    countsNumber<-agilent[currentPosition]
    counts<-c(counts,import[currentPosition])
  }
    

  ### counts will be too long and the last entry is NA
  ### counts will be corrected when the number of scans
  ### is known.
  counts<-countNumber(agilent)
  ### cut away NA. 
  counts<-counts[-which(is.na(counts))]
  
  ### the second period is extracted. This one is currently
  ### used to determine the number of scans. As the useful length
  ### is not known yet, it is extracted in overlength. Na's have
  ### to be removed then
  secondPeriod<-agilent[seq(5772,2000000,4)]
  ### remove Na's
  secondPeriod<-secondPeriod[-which(is.na(secondPeriod))]
  
  ## extract the addiontal data sequence in the first period using the countNumber function
  ## hmm, this really could be the TIC
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
  firstPeriod<-agilent[seq(5771,rawExtractLength,4)]
  secondPeriod<-agilent[seq(5772,rawExtractLength,4)]
  thirdPeriod<-agilent[seq(5773,rawExtractLength,4)]
  fourthPeriod<-c(agilent[seq(5770,rawExtractLength,4)][-1],0)
  
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
  mainFirst<-mainSequence(firstPeriod,counts)
  mainSecond<-mainSequence(secondPeriod,counts)
  mainThird<-mainSequence(thirdPeriod,counts)
  mainFourth<-mainSequence(fourthPeriod,counts)
  
  ### calculate MZs
  importMz<-round(mainFirst*12.8+mainSecond*0.05)
  
  ### calculate intensity values
  mainFourth[which(floor(mainThird/64)==3)]<-(mainFourth[which(floor(mainThird/64)==3)]*512) ### 512 is 9 full bit
  mainFourth[which(floor(mainThird/64)==2)]<-(mainFourth[which(floor(mainThird/64)==2)]*64) ### 64 is 6 full bit
  mainFourth[which(floor(mainThird/64)==1)]<-(mainFourth[which(floor(mainThird/64)==1)]*8) ### 8 is 3 full bit
  
  mainThird[which(floor(mainThird/64)==3)]<-((mainThird[which(floor(mainThird/64)==3)] %% 192)*131072) ## 131072 is 17 full bit
  mainThird[which(floor(mainThird/64)==2)]<-((mainThird[which(floor(mainThird/64)==2)] %% 128)*16384) ### 16384 is 14 full bit
  mainThird[which(floor(mainThird/64)==1)]<-((mainThird[which(floor(mainThird/64)==1)] %% 64)*2048) ### 2048 is 11 full bit
  mainThird[which(floor(mainThird/64)==0)]<-(mainThird[which(floor(mainThird/64)==0)]*256)
  
  importInt<-(mainThird+mainFourth)
  
  ### assign data to a Scan by Mz Matrix
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
      
  return(fullData)

}


