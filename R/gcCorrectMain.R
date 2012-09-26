gcCorrectMain <- function(Ms, chr, starts, samplechr, nodes, increms,
			  maxwins, jittercorrection=FALSE,
			  returnOnlyTV=FALSE, build, verbose=FALSE){
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

	uniqchrs <- unique(chr)
	remainingChr <- uniqchrs[is.na(match(uniqchrs, samplechr))]
	narrays <- ncol(Ms)

	increm <- increms[1]
	increm2 <- increms[2]
	maxwin <- maxwins[1]
	maxwin2 <- maxwins[2]
	nparts <- maxwin / increm
	if(nparts != (maxwin2 / increm2)){
		readline('Bad Increment Values')
	}
	strategyuse <- 2
	## RS: I think you need to have the annotation package in your foreach call
	gcFracBoth <- gcFracAllWin(maxwin, maxwin2, increm, increm2, chr, starts, samplechr, uniqchrs, strategyuse, annotation.pkg=pkgname,
				   verbose=verbose)

	useM <- as(Ms[chr %in% samplechr, ], "matrix")

	## first TVscore
	tvScore <- calctvScore(gcFracBoth, samplechr, nparts, useM, narrays, increm, increm2, verbose=verbose)
	dimnames(tvScore) <- list(c(seq(increm, maxwin, increm), seq(increm2, maxwin2, increm2)),
				  colnames(Ms))
	if(returnOnlyTV) {
		result <- tvScore
	} else {
		## CHECK TO SEE IF TVSCORES ARE THE SAME AS ANTICIPATED
		maxTVscores <- rep(0, narrays)
		correctionVals <- list()
		optWins <- rep(0, narrays)

		gmaxvals <- apply(tvScore, 2, max)
		gmaxvalsInd <- apply(tvScore, 2, which.max)

		## Plot TV vals ##
		## if(outputTVs){
		## plotTV(tvScore,tvScore2,narrays,increm)
		correctedM <- CorrectM(gcFracBoth, useM, Ms, starts, narrays, nparts,
				       chr, samplechr, remainingChr,
				       increm, increm2, gmaxvals, gmaxvalsInd,
				       tvScore, jittercorrection,
				       verbose=verbose)
		result <- correctedM
	}
	if(returnOnlyTV) {
		rownames(result) <- paste(c(seq(increm, maxwin, increm), seq(increm2, maxwin2, increm2)), "bp", sep="")
	}
	return(result)
}
