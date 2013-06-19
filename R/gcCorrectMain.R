gcCorrectMain <- function(Ms, chr, starts, samplechr='',
			  increms=0,
			  maxwins,
			  jittercorrection=FALSE,
			  returnOnlyTV=FALSE,
			  onlyGC=FALSE,
			  providedGC=NULL,
			  build, verbose=FALSE){
	do.call(library, list(paste("BSgenome.Hsapiens.UCSC.", build, sep='')))
	pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
	## EITAN:  These should be declared in the DESCRIPTION
	## library('RColorBrewer')
	## library('multicore')
	## library('foreach')
	## library('grDevices')
	## library('doMC')
	## library('DNAcopy')
	## registerDoMC(nodes) ## do at command line
	narrays <- ncol(Ms)
	if(is.null(providedGC)){ ##otherwise, its a vector of GC scores
		uniqchrs <- unique(chr)
                if (samplechr=='') samplechr <- uniqchrs
		remainingChr <- uniqchrs[is.na(match(uniqchrs, samplechr))]
                if(increms[1]==0) increms <- maxwins
		nparts <- maxwins[1] / increms[1]
		for(ii in seq_along(increms)){
			## Eitan: why not stop("Bad increment values") ?
			##if(nparts != maxwins[ii]/increms[ii]){readline('Bad Increment Values')}
			if(nparts != maxwins[ii]/increms[ii]) stop("maxwins/increms must have a single unique value")
		}

		strategyuse <- 2
		## RS: I think you need to have the annotation package in your foreach call
		gcFracBoth <- gcFracAllWin(maxwins, increms, chr,
					   starts, samplechr, uniqchrs,
					   strategyuse,
					   annotation.pkg=pkgname,
					   verbose=verbose)
		##useM <- as(Ms[chr %in% samplechr, ], "matrix")
		useM <- Ms[chr %in% samplechr, , drop=FALSE]
		## first TVscore
		tvScore <- calctvScore(gcFracBoth, samplechr, nparts, useM, narrays, increms, verbose=verbose)
		rown=list()
		for(ii in seq_along(increms)){
			rown[[ii]]=seq(increms[ii],maxwins[ii],increms[ii])
		}
		dimnames(tvScore) <- list(unlist(rown),colnames(Ms))
		if(returnOnlyTV)  return(tvScore)
		## CHECK TO SEE IF TVSCORES ARE THE SAME AS ANTICIPATED
		maxTVscores <- rep(0, narrays)
		correctionVals <- list()
		optWins <- rep(0, narrays)
		gmaxvals <- apply(tvScore, 2, max)
		gmaxvalsInd <- apply(tvScore, 2, which.max)
		priorFracWremaining <- getGCinBestWindowGenome(gcFracBoth,
							       starts,nparts,
							       chr,
							       samplechr,
							       remainingChr,
							       increms,
							       gmaxvals,
							       gmaxvalsInd,
							       tvScore,
							       verbose=verbose)
		if(onlyGC) result <- priorFracWremaining
	} else{
		priorFracWremaining <- providedGC
		userProvidedGC <- TRUE;gmaxvalsInd <- 0;tvScore <- 0;gcFracBoth <- 0;nparts <- 0;samplechr <- ''
		increm <- 0; increm2 <- 0;
	}
	correctedM <- CorrectM(Ms,chr,starts,priorFracWremaining,narrays,userProvidedGC,gmaxvalsInd,jittercorrection,
			       tvScore,gcFracBoth,nparts,samplechr,increms,verbose=verbose)
	##        (gcFracBoth, useM, Ms, starts, narrays, nparts,
	##			       chr, samplechr, remainingChr,
	##		       increms, gmaxvals, gmaxvalsInd,
	##	       tvScore, jittercorrection,
	##       verbose=verbose)
	result <- correctedM

	##	if(returnOnlyTV) {
	##		rownames(result) <- paste(c(seq(increm, maxwin, increm), seq(increm2, maxwin2, increm2)), "bp", sep="")
	##	}
	if(!is.null(colnames(Ms))) colnames(result) <- colnames(Ms)
        if(is.list(priorFracWremaining)) names(priorFracWremaining) <- colnames(Ms)
	resultList<-list(tvScore,optimalwins<-as.numeric(rownames(tvScore))[gmaxvalsInd],
                         maxTVscores<-gmaxvals,gcVals<-priorFracWremaining, result)
        names(resultList)<-c('tvScore','optimalWin','maxTVscore','GCvals','correctedVals')
        return(resultList)
}


##impleGcCorrect <- function(){
##	## why should I have to calculate the tvScore
##	correctedM <- CorrectM(Ms,
##			       chr,
##			       starts,
##			       priorFracWremaining,narrays,
##			       userProvidedGC,
##			       gmaxvalsInd,
##			       jittercorrection,
##			       tvScore,
##			       gcFracBoth,
##			       nparts,
##			       samplechr,
##			       increms,
##			       verbose=verbose)
##


