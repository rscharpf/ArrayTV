CorrectM <- function(gcFracBoth, useM, Ms, starts, narrays, nparts, chr,
		     samplechr, remainingChr, increm, increm2, gmaxvals, gmaxvalsInd, tvScore,
		     jittercorrection, verbose=FALSE){
	priorFracWremaining <- foreach(i=which(!duplicated(gmaxvalsInd)), .combine='list', .multicombine=TRUE, .packages="ArrayTV") %dopar% {
		## print highest tv ##
		if(verbose) print(paste('A maximum first pass  TV is', gmaxvals[i], 'in window', rownames(tvScore)[gmaxvalsInd[i]]))
		## maximum tv score window
		maxuse1 <- gmaxvalsInd[i]

		priorFrac <- priorFracs(gcFracBoth, maxuse1, nparts, tvScore,
					increm, increm2)
		if(length(remainingChr)==0){
			priorFracWremaining <- priorFrac
		}else{
			## Get GC for locations from non-sampled regions of genome
			forwardExtend <- as.numeric(rownames(tvScore)[maxuse1])
			reverseExtend <- forwardExtend

			priorFracWremaining <- priorFracsRestOfGenome(forwardExtend, reverseExtend, remainingChr, starts, chr)
			if(verbose) print(paste('forward extend', forwardExtend,'reverse extend', reverseExtend))
			priorFracWremaining[chr %in% samplechr] <- priorFrac
		}
		priorFracWremaining
	}


	correctedM <- foreach(i=seq_along(narrays), .combine='cbind', .packages="ArrayTV") %dopar% {
		if(narrays==1 | !is.list(priorFracWremaining)){
			priorFracWremainingUse <- priorFracWremaining
		}else{
			priorFracWremainingUse <- priorFracWremaining[[match(gmaxvalsInd[i], unique(gmaxvalsInd))]]
		}
		## ## this code is for debugging only ###
		## #newsp=split(useM[,i],priorFracWremainingUse[chr %in% samplechr])
		## ### this is no good, change back to the mean ###
		## #correctionVals=sapply(newsp,function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2) & x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));})
		## #correctionVals=sapply(newsp,mean)#function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2) & x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));})
		## #names(correctionVals)=names(newsp)
		## ## verify correct TV ##
		## #fsampled=sum(useM[,i]);
		## #n=nrow(useM)
		## #lambda=fsampled/n
		## #ngc=sapply(newsp,length)
		## #tvscore=sum(ngc/n * (abs(correctionVals-lambda)))
		## #print(paste('i is now',i,'tv score is',tvscore))
		## ## end debugging code
		##
		## Lets Try to Remove Large CNVs Before Calculating the Corrections, otherwise they will Interfere
		cm <- cumsum(Ms[, i])
		cm2 <- cm[200:length(cm)] - c(0, cm[seq_len(length(cm)-200)])
		cm2 <- cm2
		inds <- seq(1, length(cm2), 200)
		cm3 <- cm2[inds]
		chrForCNA <- chr[inds]
		startsForCNA <- starts[inds]
		CNA.object <- CNA(cm3, chrForCNA, startsForCNA, data.type = c("logratio"), sampleid=i)
		## smoothed.CNA.object <- smooth.CNA(CNA.object)
		segment.smoothed.CNA.object <- segment(CNA.object, verbose = as.numeric(verbose), min.width=5)
		markrep <- rep(segment.smoothed.CNA.object$output$seg.mean,segment.smoothed.CNA.object$output$num.mark)
		chromrep <- rep(segment.smoothed.CNA.object$output$chrom,segment.smoothed.CNA.object$output$num.mark)
		meduse <- mean(markrep)
		maduse <- mad(markrep)
		blockchrom <- table(segment.smoothed.CNA.object$output$chrom)
		multiblock <- names(blockchrom)[blockchrom > 1]
		chromMeans <- aggregate(markrep, list(chromrep), mean)
		segment.smoothed.CNA.object$output$chromMean <- chromMeans$x[match(segment.smoothed.CNA.object$output$chrom,chromMeans$Group.1)]
		toremove <- which(abs(segment.smoothed.CNA.object$output$seg.mean - segment.smoothed.CNA.object$output$chromMean) > maduse & segment.smoothed.CNA.object$output$chrom %in% multiblock & ( segment.smoothed.CNA.object$output$seg.mean > meduse + (maduse*2.5) | segment.smoothed.CNA.object$output$seg.mean < meduse - (maduse*2.5 )))
		##
		locsAsRange <- GRanges(seqnames=chr, ranges=IRanges(start=starts, width=1))
		removestart <- segment.smoothed.CNA.object$output$loc.start[toremove]
		removeend <- segment.smoothed.CNA.object$output$loc.end[toremove]
		removechr <- segment.smoothed.CNA.object$output$chrom[toremove]
		removeAsRange <- GRanges(seqnames=removechr, ranges=IRanges(start=removestart, end=removeend))
		##
		tokeepLogical <- is.na(GenomicRanges::match(locsAsRange, removeAsRange))
		tokeep <- which(tokeepLogical)
		toremove <- which(!tokeepLogical)
		if(length(toremove) > length(tokeep)){
			tokeep <- seq_along(tokeepLogical)
			toremove <- numeric()
		}
		##
		## calculate new TV, We only used a sample of the data to get the best window but we will use all the data to get correction values  ##
		newsp <- split(Ms[tokeep, i], paste(chr[tokeep], priorFracWremainingUse[tokeep], sep='.'))
		##
		## this will be only for CNVs that are removed initially and will not be able to match on chr and gc Frac
		correctionVals <- sapply(newsp, mean)
		## correctionVals=sapply(newsp,function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2) & x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));})
		names(correctionVals) <- names(newsp)
		fsampled <- sum(Ms[tokeep, i])
		##
		n <- length(tokeep)
		lambda <- fsampled/n
		ngc <- sapply(newsp, length)
		tvscore <- sum(ngc/n * (abs(correctionVals - lambda)))
		if(verbose) print(paste('for array', i, 'tv score is', tvscore, 'when correction is applied to each chromosome'))
		##
		allcorrections <- vector()
		allcorrections[tokeep] <- correctionVals[match(paste(chr[tokeep], priorFracWremainingUse[tokeep], sep='.'), (names(correctionVals)))]
		allcorrections[toremove] <- correctionVals[match(paste(chr[toremove], priorFracWremainingUse[toremove], sep='.'), (names(correctionVals)))]
		## this should retain information
		##
		if(verbose) print(paste('we removed', length(toremove), 'locations before calculating correction values'))
		if(verbose) print(paste(length(which(is.na(allcorrections[toremove]))), 'point(s) for array', i, 'will recieve no gc correction'))
		allcorrections[toremove][which(is.na(allcorrections[toremove]))] <- lambda
		##
		## IF WE INCREASE THE SPREAD OF OUR CORRECTIONS DOES THE AUTOCORRELATION IMPROVE? ###
		if(jittercorrection){
			cnr <- 1
			asum <- vector()
			for(jj in seq(1,1.2,.05)){
				autocor <- acf(Ms[chr==samplechr[1], i] - allcorrections[chr==samplechr[1]]*jj, lag.max=200, plot=F)
				asum[cnr] <- sum(autocor$acf[2:200])
				if(verbose) print(paste('autocorrelation sum for correction factor', jj, 'is', asum[cnr]));
				cnr <- cnr+1
			}
			mfact <- seq(1, 1.2, .05)[which.min(asum)]
		}else mfact <- 1
		chromMediansAll <- aggregate(Ms[,i], list(chr), median)
		correctedM <- chromMediansAll$x[match(chr, chromMediansAll$Group.1)] + Ms[, i] - allcorrections*mfact
		##
		newTVscore <- vector()
		for(tvind in seq_len(nrow(tvScore))){
			priorFrac <- priorFracs(gcFracBoth, tvind, nparts, tvScore, increm, increm2)
			newTVscore[tvind] <- correctionTVscore(correctedM[chr %in% samplechr], priorFrac, i, as.numeric(rownames(tvScore)[tvind]))
		}
		names(newTVscore) <- rownames(tvScore)
		##
		if(length(which(is.na(correctedM))) > 0){
			if(verbose) print(priorFracWremaining[is.na(correctedM)][1:20])
		}
		##
		## comparison metrics
		if(verbose) print(paste('mad of corrected array', i, 'is', mad(correctedM)))
		if(verbose) print(paste('old mad array', i, 'is', mad(Ms[, i])))
		if(verbose) print(paste('The new Maximum tv Score for array', i, 'is', max(newTVscore), 'in window', names(newTVscore)[which.max(newTVscore)], 'bp'))
		##
		##
		## for(ii in unique(chr)){plot(Ms[chr==ii],ylim=c(-.8,.8));plot(correctedM[chr==ii],ylim=c(-.8,.8));readline('')}
		as(correctedM, "matrix")
	}
	return(correctedM)
}
