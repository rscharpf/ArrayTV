calctvScore <-
function(gcFracBoth,samplechr,nparts,useM,narrays){
optnum=4
narrays=ncol(useM)
  ## Again, optimized for 4 nodes
### the order of row binding is forward first half, forward second half, reverse first half, reverse second half
tvScore=foreach(vv=1:4,.combine=rbind) %dopar% {

    tvScorepart=matrix(0,nrow=nparts/optnum,ncol=narrays)
  partitionedgcFrac=gcFracBoth[[vv]]


  if(vv%in% c(1,3)){
     priorFrac=0
     priorWeight=0
  }else if(vv %in% c(2,4)){    
    priorGC=gcFracBoth[[vv-1]]
    priorFrac=(cumsum(priorGC)[(nparts/optnum):length(priorGC)]-c(0,cumsum(priorGC)[1:(length(priorGC)-(nparts/optnum))]))[seq(1,(length(priorGC)-(nparts/optnum)+1) ,nparts/optnum)]/(nparts/optnum)
    priorWeight=nparts/optnum    
  }
### maxwin/increm * 2 must be divisible by 4 right now
 for(hh in 1:(nparts/optnum)){

  if(hh==1){
  gcFracA=(partitionedgcFrac[seq(hh,length(partitionedgcFrac),nparts/optnum)]+(priorFrac*priorWeight))/(priorWeight+1)
  }else{
  gcFracA=(gcFracA*((priorWeight+hh-1)/(priorWeight+hh)))+partitionedgcFrac[seq(hh,length(partitionedgcFrac),nparts/optnum)]/(priorWeight+hh)
  }

  
  newsp=split(useM,round(gcFracA,2))

  toplot=sapply(newsp,function(x){cc=matrix(x,ncol=narrays);y=colMeans(cc);y})


  rownames(toplot)=colnames(useM);
  colnames(toplot)=names(newsp)
  fsampled=colSums(useM)
  n=nrow(useM)
  lambda=fsampled/n
  ngc=sapply(newsp,length)/narrays
  tvscore=rowSums(rep((ngc/n),each=narrays)*(abs(toplot-lambda)))

  tvScorepart[hh,]=na.omit(tvscore)
}

  tvScorepart
  }
return(tvScore)
}
