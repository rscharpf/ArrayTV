CorrectM <-
function(gcFracBoth,useM,Ms,starts,narrays,nparts,chr,samplechr,remainingChr,increm,formaxs1,formaxs2,revmaxs1,revmaxs2,maxorigins,gmaxvals,gmaxvals2){

correctedM=foreach(i=1:narrays, .combine='cbind')%dopar% {

## print highest tv ##
print(paste('i is', i,'and high tv is',ifelse(maxorigins[i]==F,gmaxvals2[i],gmaxvals[i])))

ar=i

## Get GC Fraction for each location from sampled data
reversedirMax2=ifelse(maxorigins[i]==F,ifelse(revmaxs2[i,2] > formaxs2[i,2],T,F),F)
  
## Use First tvScore matrix to decide which direction to extend first
reversedirMax= revmaxs1[i,2] > formaxs1[i,2] & reversedirMax2 == F  #gmaxs[i]>(nrow(tvScore)/2)
maxuse1=ifelse(reversedirMax,revmaxs1[i,1],formaxs1[i,1])
vv1=ifelse(reversedirMax,2,1)
v1=c(3,1);v2=c(4,2)
priorFrac=priorFracs(gcFracBoth,maxuse1,vv1,nparts,v1,v2)

if(maxorigins[i]==F){
### Extend in a second direction if necessary  
maxuse2=ifelse(reversedirMax2,revmaxs2[i,1],formaxs2[i,1])
vv2=ifelse(reversedirMax2,2,1)
priorFrac2=priorFracs(gcFracBoth,maxuse2,vv2,nparts,v1,v2)
priorFrac=priorFrac*(maxuse1/(maxuse1+maxuse2))+priorFrac2*(maxuse2/(maxuse1+maxuse2))
}
priorFrac=round(priorFrac,2)

newsp=split(useM[,i],priorFrac)
#correctionVals=sapply(newsp,median)
correctionVals=sapply(newsp,mean)
names(correctionVals)=names(newsp)
## verify correct TV ##
fsampled=sum(useM[,i]);
 n=nrow(useM)
  lambda=fsampled/n
ngc=sapply(newsp,length)
tvscore=sum(ngc/n * (abs(correctionVals-lambda)))
print(paste('i is now',i,'tv score is',tvscore))

if(length(remainingChr)==0){
  priorFracWremaining=priorFrac
}else{
## Get GC for locations from non-sampled regions of genome
forwardExtend=ifelse(!reversedirMax,maxuse1*increm,ifelse(!reversedirMax2 && maxorigins[i]==F,maxuse2*increm,1)) ## must include at least the starting position
reverseExtend=ifelse(reversedirMax,maxuse1*increm,ifelse(reversedirMax2,maxuse2*increm,0))

schr=chr[chr %in% samplechr]
#priorFracWremaining=round(priorFracsRestOfGenome(forwardExtend,reverseExtend,samplechr,startToCalcRemainingGC[chr %in% samplechr],schr),2)#remainingChr,startToCalcRemainingGC,chr),2)
priorFracWremaining=round(priorFracsRestOfGenome(forwardExtend,reverseExtend,remainingChr,starts,chr),2)
print(paste('forward extend',forwardExtend,'reverse extend',reverseExtend))
#print(length(which(priorFracWremaining[chr %in% samplechr]>0)))


priorFracWremaining[chr %in% samplechr]=priorFrac
}

## calculate new TV, We only used a sample of the data to get the best window but we will use all the data to get correction values  ##
newsp=split(Ms[,i],paste(chr,priorFracWremaining,sep='.'))
#correctionVals=sapply(newsp,median)
correctionVals=sapply(newsp,mean)
names(correctionVals)=names(newsp)
fsampled=sum(Ms[,i]);
 n=nrow(Ms)
  lambda=fsampled/n
ngc=sapply(newsp,length)
tvscore=sum(ngc/n * (abs(correctionVals-lambda)))
print(paste('i is now',i,'NEW* tv score is',tvscore))

correctedM=Ms[,i]-correctionVals[match(paste(chr,priorFracWremaining,sep='.'),(names(correctionVals)))]

if(length(which(is.na(correctedM)))>0){
  print(priorFracWremaining[is.na(correctedM)][1:20])
}

## comparison metrics
print(paste('mad is',mad(correctedM)))

print(paste('old mad is',mad(Ms[,i])))
###

as(correctedM,"matrix")
}
return(correctedM)
}
