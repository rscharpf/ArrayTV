gcFracAllWin <-
function(maxwin,maxwin2,increm,increm2,chr,allstarts,samplechr,uniqchrs,strategyuse){
  ### This is currently optimizied for 4 nodes, this should be made
  ### more generalizeable if someone has more than 4 nodes, for instance
  #gcFracBoth=foreach(vv=1:2,.combine=list,.multicombine=T) %:%
  gcFracBoth=foreach(schr=uniqchrs[uniqchrs %in% samplechr],.combine=rbind) %dopar% {
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
schr=gsub('23','X',schr)
schr=gsub('24','X',schr)

print('Getting gc content From BS genome Object')
pregcFrac=letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=increm,'CG',as.prob=T)
#sgc=cumsum(pregcFrac)
pregcFrac2=letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=increm2,'CG',as.prob=T)
#pregcFrac2=sgc[(increm2/increm):length(sgc)]-c(0,sgc[1:(length(sgc)-(increm2/increm))])
print('gc content stored')


startinds=rep(starts,each=maxwin/increm)+rep(seq(0,maxwin-increm,increm),length(starts))
startinds2=rep(starts,each=maxwin2/increm2)+rep(seq(0,maxwin2-increm2,increm2),length(starts))  

startindbackwards=rep(starts,each=maxwin/increm)-rep(seq(increm,maxwin,increm),length(starts))
startindbackwards2=rep(starts,each=maxwin2/increm2)-rep(seq(increm2,maxwin2,increm2),length(starts))



startinds[which(startinds<1)]=1
startinds[which(startinds>length(pregcFrac))]=length(pregcFrac)
gcFrac1=pregcFrac[startinds]
gcFrac1[is.na(gcFrac1)]=0

startinds2[which(startinds2<1)]=1
startinds2[which(startinds2>length(pregcFrac2))]=length(pregcFrac2)
gcFrac2=pregcFrac2[startinds2]
gcFrac2[is.na(gcFrac2)]=0


startindbackwards[which(startindbackwards<1)]=1
startindbackwards[which(startindbackwards>length(pregcFrac))]=length(pregcFrac)
gcFracbackwards1=pregcFrac[startindbackwards]
gcFracbackwards1[is.na(gcFracbackwards1)]=0


startindbackwards2[which(startindbackwards2<1)]=1
startindbackwards2[which(startindbackwards2>length(pregcFrac2))]=length(pregcFrac2)
gcFracbackwards2=pregcFrac2[startindbackwards2]
gcFracbackwards2[is.na(gcFracbackwards2)]=0

gcFrac=cbind((gcFrac1+gcFracbackwards1)/2,(gcFrac2+gcFracbackwards2)/2)

}
  gcFrac
}
return(gcFracBoth)  
}
