gcCorrectMain <-
function(Ms,chr,starts,samplechr,nodes,increms,maxwins,jittercorrection=F,returnOnlyTV=F,build){
  do.call(library,list(paste("BSgenome.Hsapiens.UCSC.",build,sep='')))
  library('RColorBrewer')
  library('multicore')
  library('foreach')
  library('grDevices')
  library('doMC')
  library('DNAcopy')
  registerDoMC(nodes)

  uniqchrs=unique(chr)
  remainingChr=uniqchrs[is.na(match(uniqchrs,samplechr))]
  narrays=ncol(Ms)

  increm=increms[1]
  increm2=increms[2]
  maxwin=maxwins[1]
  maxwin2=maxwins[2]
  nparts=(maxwin/increm)
  if(nparts != maxwin2/increm2){readline('Bad Increment Values')}
  strategyuse=2

  gcFracBoth=gcFracAllWin(maxwin,maxwin2,increm,increm2,chr,starts,samplechr,uniqchrs,strategyuse)


useM=as(Ms[chr  %in% samplechr,],"matrix")


## first TVscore
tvScore=calctvScore(gcFracBoth,samplechr,nparts,useM,narrays,increm,increm2)

rownames(tvScore)=c(seq(increm,maxwin,increm),seq(increm2,maxwin2,increm2))

if(returnOnlyTV){
  return(tvScore)
}

#### CHECK TO SEE IF TVSCORES ARE THE SAME AS ANTICIPATED

maxTVscores=rep(0,narrays)
correctionVals=list()
optWins=rep(0,narrays)


gmaxvals=apply(tvScore,2,max)
gmaxvalsInd=apply(tvScore,2,which.max)

## Plot TV vals ##
#if(outputTVs){
#plotTV(tvScore,tvScore2,narrays,increm)
#}

correctedM=CorrectM(gcFracBoth,useM,Ms,starts,narrays,nparts,chr,samplechr,remainingChr,increm,increm2,gmaxvals,gmaxvalsInd,tvScore,jittercorrection)

return(correctedM)
}
