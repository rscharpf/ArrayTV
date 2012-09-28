gcFracAllWin <- function(maxwins, increms,
			 chr, allstarts,
			 samplechr, uniqchrs, strategyuse, annotation.pkg,
			 verbose=FALSE){
    incremuse <- NULL
    gcFracBoth <- foreach(schr=uniqchrs[uniqchrs %in% samplechr], .combine=rbind, .packages=c(annotation.pkg, "ArrayTV")) %:%
        foreach(incremuse=increms,.combine=cbind) %dopar% {
            starts <- allstarts[chr %in% schr]
            ## STRATEGY 1: This is not implemented
##             if(strategyuse==1){
##                 if(verbose) print(paste('vv is',vv)  )
##                 if(vv==1){
##                     ranges <- IRanges(start=rep(starts,each=maxwin/increm/2)+rep(seq(0,(maxwin/2)-increm,increm),length(starts)),
##                                       width=increm)
##                 }else if(vv==2){
##                     ranges <- IRanges(start=rep(starts,each=maxwin/increm/2)+rep(seq((maxwin/2),maxwin-increm,increm),length(starts)),
##                                       width=increm)
##                 }else if(vv==3){
##                     ranges <- IRanges(start=rep(starts,each=maxwin/increm/2)-rep(seq(increm,maxwin/2,increm),length(starts)),
##                                       width=increm)
##                 }else if(vv==4){
##                     ranges <- IRanges(start=rep(starts,each=maxwin/increm/2)-rep(seq((maxwin/2)+increm,maxwin,increm),length(starts)),
##                                       width=increm)
##                 }
##                 ranges[which(start(ranges) < 1)] <- IRanges(start=1,end=increm)
##                 if(verbose) print('before All views')
##                 inDNAset <- DNAStringSet()
##                 preinDNAset2 <- DNAStringSet(Views(Hsapiens[[schr]],ranges))
##                 preinDNAset2[which(width(preinDNAset2) < increm)] <- DNAStringSet((paste(rep('N',increm),collapse='')))
##                 inDNAset <- unlist(preinDNAset2)

##                 dnasetLen <- length(inDNAset)
##                 gcFrac <- letterFrequencyInSlidingView(inDNAset,view.width=increm,'CG',as.prob=T)[seq(1,dnasetLen,increm)]
##             }
            ## END Strategy 1

            ## STRATEGY 2
            if(strategyuse==2){
                schr <- gsub('23', 'X', schr)
                schr <- gsub('24', 'X', schr)

                if(verbose) print('Getting gc content From BS genome Object')
                maxwinuse <- maxwins[which(increms==incremuse)]
                gcFrac <- as.matrix(gcFracOneRange(schr,starts,incremuse,maxwinuse,verbose))
            } ## end strategy 2
            gcFrac
	} ## end foreach statement
    return(gcFracBoth)
}
