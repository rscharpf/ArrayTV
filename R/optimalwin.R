optimalwin <-
function(useM,gcFracBoth,tvScore,vv,ar,nparts){
optnum=4  
    v1=c(3,1);v2=c(4,2)
    
  tvScoreForwards=tvScore[1:(nrow(tvScore)/2),]
  tvScoreReverses=tvScore[((nrow(tvScore)/2)+1):nrow(tvScore),]
  formaxs=apply(tvScoreForwards,2,which.max)
  revmaxs=apply(tvScoreReverses,2,which.max)
  tvScorepart=vector()#matrix(0,nrow=nparts/nodes,ncol=8)
  maxuse=ifelse(vv==1,formaxs[ar],revmaxs[ar])
  ## start at first going forward and at first going backwards


  priorFrac=priorFracs(gcFracBoth,maxuse,vv,nparts,v1,v2)  
 priorWeight=maxuse              

 for(hh in 1:((nparts/optnum)*2)){
   if(hh > (nparts/optnum)){
     partitionedgcFrac=gcFracBoth[[v2[vv]]]
     }else{
       partitionedgcFrac=gcFracBoth[[v1[vv]]]
     }
   cat('.');
 
  if(hh==1){   
  gcFracA=(partitionedgcFrac[seq(hh,length(partitionedgcFrac),nparts/optnum)]+(priorFrac*priorWeight))/(priorWeight+1)
}else{
  gcFracA=(gcFracA*((priorWeight+hh-1)/(priorWeight+hh)))+partitionedgcFrac[seq(hh,length(partitionedgcFrac),nparts/optnum)]/(priorWeight+hh)
      }
  
  newsp=split(useM[,ar],round(gcFracA,2))
  toplot=sapply(newsp,mean)
   names(toplot)=names(newsp)
  fsampled=sum(useM[,ar]);
  n=nrow(useM)
  lambda=fsampled/n
  ngc=sapply(newsp,length)
  tvscore=sum(ngc/n * (abs(toplot-lambda)))

  tvScorepart=c(tvScorepart,tvscore)

 }
               
   tvScorepart
}
