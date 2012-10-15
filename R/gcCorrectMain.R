gcCorrectMain <- function(Ms, chr, starts, samplechr, nodes,
			  increms,
			  maxwins,
			  jittercorrection=FALSE,
			  returnOnlyTV=FALSE,
			  onlyGC=FALSE,
			  providedGC=0,
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
	if(length(providedGC)==1){

		uniqchrs <- unique(chr)
		remainingChr <- uniqchrs[is.na(match(uniqchrs, samplechr))]

		nparts <- maxwins[1] / increms[1]
		for(ii in seq_along(increms)){
			if(nparts != maxwins[ii]/increms[ii]){readline('Bad Increment Values')}
		}

		strategyuse <- 2
		## RS: I think you need to have the annotation package in your foreach call
		gcFracBoth <- gcFracAllWin(maxwins, increms, chr, starts, samplechr, uniqchrs, strategyuse, annotation.pkg=pkgname,
					   verbose=verbose)

		useM <- as(Ms[chr %in% samplechr, ], "matrix")

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


		## RS: I think you need the triple colon for parallelization if this function is not exported in the package namespace
		priorFracWremaining <- ArrayTV:::getGCinBestWindowGenome(gcFracBoth,starts,nparts,chr,samplechr,remainingChr,increms,gmaxvals,gmaxvalsInd,
									 tvScore,verbose=verbose)

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
	return(result)
}


