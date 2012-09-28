correctionTVscore <-
function(useM,gcFrac,arrayn,windowl){
  newsp=split(useM,gcFrac)
  correctionVals=sapply(newsp,mean)
  fsampled=sum(useM);
  n=length(useM)
  lambda=fsampled/n
  ngc=sapply(newsp,length)
  tvscore=sum(ngc/n * (abs(correctionVals-lambda)))
  #print(paste('for array',arrayn,'new tv for window',windowl,'is',tvscore))
  return(tvscore)
}
