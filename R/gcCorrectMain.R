gcCorrectMain <-
function(Ms,chr,starts,samplechr,nodes,increms,maxwins,outputTVs,outputChrReg,platform,build){
  do.call(library,list(paste("BSgenome.Hsapiens.UCSC.",build,sep='')))
  library('RColorBrewer')
  library('multicore')
  library('foreach')
  library('grDevices')
  library('doMC')

  registerDoMC(nodes)
if(exists('correctedM')){
  rm(correctedM)}
  


if(outputTVs & outputChrReg){
pdf(paste('chrViews',platform,'.pdf',sep=''))
}

  
  
#  samplechr=paste('chr',c(1,4,8,10),sep='')
#  samplechr=paste('chr',c(1,4,8,10),sep='')  
  uniqchrs=unique(chr)
  remainingChr=uniqchrs[is.na(match(uniqchrs,samplechr))]
  narrays=ncol(Ms)

for(correctind in 1:(length(increms))){
repeat{

#if(exists('correctedM')){
#  oldMs=Ms  ## this will hold onto the original
#  Ms=correctedM
#}
  increm=increms[correctind]
  maxwin=maxwins[correctind]
  nparts=(maxwin/increm)*2
  strategyuse=2

gcFracBoth=gcFracAllWin(maxwin,increm,chr,starts,samplechr,uniqchrs,strategyuse)


useM=as(Ms[chr  %in% samplechr,],"matrix")


## first TVscore  
tvScore=calctvScore(gcFracBoth,samplechr,nparts,useM,narrays)
rownames(tvScore)=c(seq(increm,maxwin,increm),-seq(increm,maxwin,increm))
assign("tvScore",tvScore,envir=.GlobalEnv)


## Second TVscore
tvScore2=calctvScore2(useM,gcFracBoth,tvScore,nparts,narrays)
assign("tvScore2",tvScore2,envir=.GlobalEnv)

#### CHECK TO SEE IF TVSCORES ARE THE SAME AS ANTICIPATED

maxTVscores=rep(0,narrays)
correctionVals=list()
optWins=rep(0,narrays)


## ProbeCorrect to get TV scores ####
formaxs1=tvScoreFirstHalfMax(tvScore)
revmaxs1=tvScoreSecondHalfMax(tvScore)

formaxs2=tvScoreSecondHalfMax(tvScore2)
revmaxs2=tvScoreFirstHalfMax(tvScore2)

### Extend in opposite Directions ####

gmaxvals=apply(tvScore,2,max)
gmaxvals2=apply(tvScore2,2,max)
maxorigins=gmaxvals >= gmaxvals2

maxoverall=max(c(gmaxvals,gmaxvals2))
  if(maxoverall < .01){break}

## Plot TV vals ##
if(outputTVs){
plotTV(tvScore,tvScore2,narrays,increm)
}

correctedM=CorrectM(gcFracBoth,useM,Ms,starts,narrays,nparts,chr,samplechr,remainingChr,increm,formaxs1,formaxs2,revmaxs1,revmaxs2,maxorigins,gmaxvals,gmaxvals2)


  oldMs=Ms  ## this will hold onto the original
  Ms=correctedM


  
##compare with plots###
if(outputChrReg){


#if(exists('oldMs')){
#  Ms=oldMs
#}
for(jj in 1:narrays){
  #dev.new()
  par(mfrow=c(2,1))
  minrange=85e6
  maxrange=89e6
  chrplot=samplechr[1]
  Msplot=oldMs[starts>minrange & starts < maxrange & chr==chrplot,jj]
  ylim1=quantile(Msplot,c(.002,.998))
    
plot(starts[starts>minrange & starts < maxrange & chr==chrplot],Msplot,ylim=ylim1)
plot(starts[starts>minrange & starts < maxrange & chr==chrplot],Ms[starts>minrange & starts < maxrange & chr==chrplot,jj],ylim=ylim1)
}
}
}
####
}
graphics.off()  
return(correctedM)
}
