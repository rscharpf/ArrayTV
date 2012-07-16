gcFracAllWin <-
function(maxwin,increm,chr,allstarts,samplechr,uniqchrs,strategyuse){
  ### This is currently optimizied for 4 nodes, this should be made
  ### more generalizeable if someone has more than 4 nodes, for instance
  gcFracBoth=foreach(vv=1:4,.combine=list,.multicombine=T) %:%
  foreach(schr=uniqchrs[uniqchrs %in% samplechr],.combine=c) %dopar% {
  starts=allstarts[chr %in% schr]

### STRATEGY 1
if(strategyuse==1){
print(paste('vv is',vv)  )
  if(vv==1){
  ranges=IRanges(start=rep(starts,each=maxwin/increm/2)+rep(seq(0,(maxwin/2)-increm,increm),length(starts)),
         width=increm)
}else if(vv==2){
  ranges=IRanges(start=rep(starts,each=maxwin/increm/2)+rep(seq((maxwin/2),maxwin-increm,increm),length(starts)),
         width=increm)
  }else if(vv==3){
ranges=IRanges(start=rep(starts,each=maxwin/increm/2)-rep(seq(increm,maxwin/2,increm),length(starts)),
  width=increm)
}else if(vv==4){
ranges=IRanges(start=rep(starts,each=maxwin/increm/2)-rep(seq((maxwin/2)+increm,maxwin,increm),length(starts)),
  width=increm)
}
  
  ranges[which(start(ranges)<1)]=IRanges(start=1,end=increm)

  print('before All views')
  inDNAset=DNAStringSet()

    preinDNAset2=DNAStringSet(Views(Hsapiens[[schr]],ranges))
    preinDNAset2[which(width(preinDNAset2)<increm)]=DNAStringSet((paste(rep('N',increm),collapse='')))
   inDNAset=unlist(preinDNAset2)

    
  dnasetLen=length(inDNAset)

  gcFrac= letterFrequencyInSlidingView(inDNAset,view.width=increm,'CG',as.prob=T)[seq(1,dnasetLen,increm)]
}
### END Strategy 1  
  
### STRATEGY 2
if(strategyuse==2){  
pregcFrac=letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=increm,'CG',as.prob=T)
if(vv==1){
startinds=rep(starts,each=maxwin/increm/2)+rep(seq(0,(maxwin/2)-increm,increm),length(starts))  
}else if(vv==2){
startinds=rep(starts,each=maxwin/increm/2)+rep(seq((maxwin/2),maxwin-increm,increm),length(starts))  
}else if(vv==3){
startinds=rep(starts,each=maxwin/increm/2)-rep(seq(increm,maxwin/2,increm),length(starts))  
}else if(vv==4){
startinds=rep(starts,each=maxwin/increm/2)-rep(seq((maxwin/2)+increm,maxwin,increm),length(starts))  
}
print(paste('gc window chunk',vv))
startinds[which(startinds<1)]=1
startinds[which(startinds>length(pregcFrac))]=length(pregcFrac)
gcFrac=pregcFrac[startinds]
gcFrac[is.na(gcFrac)]=0

}
  gcFrac
}
return(gcFracBoth)  
}
