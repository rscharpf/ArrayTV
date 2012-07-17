priorFracsRestOfGenome <-
function(forwardExtend,reverseExtend,chrToInvestigate,locs,chrs){
gcFrac=rep(0,length(locs))
winlen=reverseExtend+forwardExtend
  for(x in chrToInvestigate){
  li=which(chrs==x)
  ll=locs[li]
  schr=ifelse(grep('chr',x),x,paste('chr',x,sep=''))
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
gcFrac[li]=gcIndivChr  
}  
gcFrac
}
