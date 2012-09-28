getGCinBestWindowGenome <-
    function(gcFracBoth,starts,nparts,chr,samplechr,remainingChr,increms,gmaxvals,gmaxvalsInd,tvScore,verbose=FALSE){
        i <- NULL
        priorFracWremaining=foreach(i=which(!duplicated(gmaxvalsInd)), .combine='list',.multicombine=T)%dopar% {

        if(verbose) print(paste('A maximum first pass  TV is',gmaxvals[i],'in window',rownames(tvScore)[gmaxvalsInd[i]]))

        ## maximum tv score window
        maxuse1=gmaxvalsInd[i]

        priorFrac=priorFracs(gcFracBoth,maxuse1,nparts,tvScore,increms)

        if(length(remainingChr)==0){
            priorFracWremaining=priorFrac
        }else{
            ## Get GC for locations from non-sampled regions of genome
            forwardExtend=as.numeric(rownames(tvScore)[maxuse1])
            reverseExtend=forwardExtend

            priorFracWremaining=priorFracsRestOfGenome(forwardExtend,reverseExtend,remainingChr,starts,chr)
            if(verbose) print(paste('forward extend',forwardExtend,'reverse extend',reverseExtend))

            priorFracWremaining[chr %in% samplechr]=priorFrac
        }

        priorFracWremaining
    }
    return(priorFracWremaining)
}
