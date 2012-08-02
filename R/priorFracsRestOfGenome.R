priorFracsRestOfGenome <-
function(forwardExtend,reverseExtend,chrToInvestigate,locs,chrs){
gcFrac=rep(0,length(locs))
winlen=reverseExtend+forwardExtend
print(paste('calculating gc content for the rest of the genome for window',forwardExtend,'bp'))
  for(x in chrToInvestigate){
  li=which(chrs==x)
  ll=locs[li]
  schr=ifelse(grep('chr',x),x,paste('chr',x,sep=''))
  schr=gsub('23','X',schr)
  schr=gsub('24','Y',schr)  
  ## old strategy
  #preinDNAset=DNAStringSet(Views(Hsapiens[[genchr]],IRanges(start=ll-reverseExtend,end=ll+forwardExtend-1)))
  #inDNAset=unlist(preinDNAset)
  #gcIndivChr=letterFrequencyInSlidingView(inDNAset,view.width=winlen,'CG',as.prob=T)[seq(1,length(inDNAset),winlen)]
  ####
  
  ## new strategy
  startinds=ll-reverseExtend
  startinds[startinds<1]=1
  startinds[startinds>(length(Hsapiens[[schr]])-winlen)]=(length(Hsapiens[[schr]])-winlen)
  if(diff(range(startinds))>=winlen){
  gcIndivChr=letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=winlen,'CG',as.prob=T)[startinds]
  gcIndivChr[is.na(gcIndivChr)]=tail(gcIndivChr,1)
}else{
   gcIndivChr=rep(letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=length(Hsapiens[[schr]]),'CG',as.prob=T),length(startinds))
 }
  ###
  
  cat('.')
  if(length(gcFrac[li]) != length(ll)){
    print('something is wrong with gc dnastring set')
  }
  quants=sort(unique(quantile(gcIndivChr,probs=seq(0,1,.01))))
gcFrac[li]=quants[findInterval(gcIndivChr,quants)]
}  
gcFrac
}
