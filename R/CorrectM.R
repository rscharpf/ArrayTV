CorrectM <- function(Ms,chr,starts,priorFracWremaining,narrays,
		     userProvidedGC=FALSE,gmaxvalsInd=0,
		     jittercorrection=FALSE,
                     tvScore,gcFracBoth,nparts,
		     samplechr,increms, verbose=FALSE){
	i <- NULL
	correctedM <- foreach(i=1:narrays, .combine='cbind', .packages="ArrayTV")%dopar% {

		if(narrays==1 | !is.list(priorFracWremaining)){
			priorFracWremainingUse <- priorFracWremaining
		}else{
			priorFracWremainingUse <- priorFracWremaining[[match(gmaxvalsInd[i],unique(gmaxvalsInd))]]
		}

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
		chromMads <- aggregate(cm3, list(chrForCNA), mad)
		segment.smoothed.CNA.object$output$chromMean <- chromMeans$x[match(segment.smoothed.CNA.object$output$chrom,chromMeans$Group.1)]
		segment.smoothed.CNA.object$output$chromMad <- chromMads$x[match(segment.smoothed.CNA.object$output$chrom,chromMads$Group.1)]
		toremove <- which(abs(segment.smoothed.CNA.object$output$seg.mean - segment.smoothed.CNA.object$output$chromMean) >
                                  segment.smoothed.CNA.object$output$chromMad & segment.smoothed.CNA.object$output$chrom %in% multiblock &
                                  ( segment.smoothed.CNA.object$output$seg.mean > meduse + (maduse*2.5) | segment.smoothed.CNA.object$output$seg.mean < meduse - (maduse*2.5 )))
                ### added this for comparison
                ##toremove2 <- which(abs(segment.smoothed.CNA.object$output$seg.mean - segment.smoothed.CNA.object$output$chromMean) > maduse & segment.smoothed.CNA.object$output$chrom %in% multiblock & ( segment.smoothed.CNA.object$output$seg.mean > meduse + (maduse*2.5) | segment.smoothed.CNA.object$output$seg.mean < meduse - (maduse*2.5 )))
                ### end added for comparison

		##
		locsAsRange <- GRanges(seqnames=chr, ranges=IRanges(start=starts, width=1))
		removestart <- segment.smoothed.CNA.object$output$loc.start[toremove]
		removeend <- segment.smoothed.CNA.object$output$loc.end[toremove]
		removechr <- segment.smoothed.CNA.object$output$chrom[toremove]
		removeAsRange <- GRanges(seqnames=removechr, ranges=IRanges(start=removestart, end=removeend))
		##tokeepLogical <- is.na(GenomicRanges::match(locsAsRange, removeAsRange,match.if.overlap=TRUE))
		tokeepLogical <- is.na(findOverlaps(locsAsRange, removeAsRange, select="first"))
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
		if(verbose) message('for array ', i, ' tv score is ', round(tvscore, 3), ' when correction is applied to each chromosome')
		##
		allcorrections <- vector()
		allcorrections[tokeep] <- correctionVals[match(paste(chr[tokeep], priorFracWremainingUse[tokeep], sep='.'), (names(correctionVals)))]
		allcorrections[toremove] <- correctionVals[match(paste(chr[toremove], priorFracWremainingUse[toremove], sep='.'), (names(correctionVals)))]
		nar<-which(is.na(allcorrections[toremove]))

		if(length(nar)>2e3){
			correctionValsNum<-correctionVals
			names(correctionValsNum)<-as.numeric(gsub('\\.0','',gsub('chr','',names(correctionValsNum))))
			correctionValsNum<-correctionValsNum[order(as.numeric(names(correctionValsNum)))]
			allcorrections[toremove][nar]=correctionValsNum[findInterval(as.numeric(paste(gsub('chr','',chr[toremove][nar]),gsub('^0','',
                        priorFracWremainingUse[toremove][nar]),sep='')),as.numeric(names(correctionValsNum)))]
			}

		## this should retain information
		##
		if(verbose) message('we removed ', length(toremove), ' locations before calculating correction values')
		if(verbose) message(length(which(is.na(allcorrections[toremove]))), ' point(s) for array ', i, ' will recieve no gc correction')
		allcorrections[toremove][which(is.na(allcorrections[toremove]))] <- lambda
		##
		## IF WE INCREASE THE SPREAD OF OUR CORRECTIONS DOES THE AUTOCORRELATION IMPROVE? ###
		if(jittercorrection){
			cnr <- 1
			asum <- vector()
			for(jj in seq(1,1.2,.05)){
				autocor <- acf(Ms[chr==samplechr[1], i] - allcorrections[chr==samplechr[1]]*jj, lag.max=200, plot=F)
				asum[cnr] <- sum(autocor$acf[2:200])
				if(verbose) message('autocorrelation sum for correction factor ', jj, ' is ', asum[cnr])
				cnr <- cnr+1
			}
			mfact <- seq(1, 1.2, .05)[which.min(asum)]
		}else mfact <- 1
		chromMediansAll <- aggregate(Ms[,i], list(chr), median)
		correctedM <- chromMediansAll$x[match(chr, chromMediansAll$Group.1)] + Ms[, i] - allcorrections*mfact
		if(verbose){
			if(tvScore[1] != 0){
				newTVscore <- vector()
				for(tvind in seq_len(nrow(tvScore))){
					priorFrac <- priorFracs(gcFracBoth, tvind, nparts, tvScore, increms)
					newTVscore[tvind] <- correctionTVscore(correctedM[chr %in% samplechr], priorFrac, i, as.numeric(rownames(tvScore)[tvind]))
				}
				names(newTVscore) <- rownames(tvScore)
			}
			##
			if(length(which(is.na(correctedM))) > 0){
				message(priorFracWremaining[is.na(correctedM)][1:20])
			}
			##
			## comparison metrics
			message('mad of corrected array ', i, ' is ', round(mad(correctedM), 3))
			message('old mad array ', i, ' is ', round(mad(Ms[, i]), 3))
			if(tvScore[1] != 0)
				message('The new Maximum tv Score for array ', i, ' is ', round(max(newTVscore), 3), ' in window ', names(newTVscore)[which.max(newTVscore)], ' bp')
			##
			##
			## for(ii in unique(chr)){plot(Ms[chr==ii],ylim=c(-.8,.8));plot(correctedM[chr==ii],ylim=c(-.8,.8));readline('')}
		}

		as(correctedM, "matrix")
	}
	return(correctedM)
}

## correctMLite <- function(gcFracBoth,
## 			 Ms,
## 			 chr,
## 			 starts,
## 			 ##priorFracWremaining,
## 			 ##narrays,
## 			 ##userProvidedGC=FALSE,
## 			 ##gmaxvalsInd=0,
## 			 ##jittercorrection=FALSE,
## 			 ##tvScore,
## 			 ##gcFracBoth,
## 			 ##nparts,
## 			 samplechr,
## 			 increms,
## 			 verbose=FALSE){
## 	i <- NULL
## 	narrays <- ncol(Ms)
## 	correctedM <- foreach(i=seq_len(narrays), .combine='cbind', .packages="ArrayTV")%dopar% {
## ##		if(narrays==1 | !is.list(priorFracWremaining)){
## ##			priorFracWremainingUse <- priorFracWremaining
## ##		}else{
## ##			priorFracWremainingUse <- priorFracWremaining[[match(gmaxvalsInd[i],unique(gmaxvalsInd))]]
## ##		}

## 		## Lets Try to Remove Large CNVs Before Calculating the Corrections, otherwise they will Interfere
## 		cm <- cumsum(Ms[, i])
## 		cm2 <- cm[200:length(cm)] - c(0, cm[seq_len(length(cm)-200)])
## 		cm2 <- cm2
## 		inds <- seq(1, length(cm2), 200)
## 		cm3 <- cm2[inds]
## 		chrForCNA <- chr[inds]
## 		startsForCNA <- starts[inds]
## 		CNA.object <- CNA(cm3, chrForCNA, startsForCNA, data.type = c("logratio"), sampleid=i)
## 		## smoothed.CNA.object <- smooth.CNA(CNA.object)
## 		segment.smoothed.CNA.object <- segment(CNA.object, verbose = as.numeric(verbose), min.width=5)
## 		markrep <- rep(segment.smoothed.CNA.object$output$seg.mean,segment.smoothed.CNA.object$output$num.mark)
## 		chromrep <- rep(segment.smoothed.CNA.object$output$chrom,segment.smoothed.CNA.object$output$num.mark)
## 		meduse <- mean(markrep)
## 		maduse <- mad(markrep)
## 		blockchrom <- table(segment.smoothed.CNA.object$output$chrom)
## 		multiblock <- names(blockchrom)[blockchrom > 1]
## 		chromMeans <- aggregate(markrep, list(chromrep), mean)
## 		segment.smoothed.CNA.object$output$chromMean <- chromMeans$x[match(segment.smoothed.CNA.object$output$chrom,chromMeans$Group.1)]
## 		toremove <- which(abs(segment.smoothed.CNA.object$output$seg.mean - segment.smoothed.CNA.object$output$chromMean) > maduse & segment.smoothed.CNA.object$output$chrom %in% multiblock & ( segment.smoothed.CNA.object$output$seg.mean > meduse + (maduse*2.5) | segment.smoothed.CNA.object$output$seg.mean < meduse - (maduse*2.5 )))
## 		##
## 		locsAsRange <- GRanges(seqnames=chr, ranges=IRanges(start=starts, width=1))
## 		removestart <- segment.smoothed.CNA.object$output$loc.start[toremove]
## 		removeend <- segment.smoothed.CNA.object$output$loc.end[toremove]
## 		removechr <- segment.smoothed.CNA.object$output$chrom[toremove]
## 		removeAsRange <- GRanges(seqnames=removechr, ranges=IRanges(start=removestart, end=removeend))
## 		##
## 		##tokeepLogical <- is.na(GenomicRanges::match(locsAsRange, removeAsRange,match.if.overlap=TRUE))
## 		tokeepLogical <- is.na(findOverlaps(locsAsRange, removeAsRange, select="first"))
## 		tokeep <- which(tokeepLogical)
## 		toremove <- which(!tokeepLogical)
## 		if(length(toremove) > length(tokeep)){
## 			tokeep <- seq_along(tokeepLogical)
## 			toremove <- numeric()
## 		}
## 		##
## 		## calculate new TV, We only used a sample of the data
## 		## to get the best window but we will use all the data
## 		## to get correction values
## 		newsp <- split(Ms[tokeep, i], paste(chr[tokeep], priorFracWremainingUse[tokeep], sep='.'))
## 		##
## 		## this will be only for CNVs that are removed initially and will not be able to match on chr and gc Frac
## 		## correctionVals <- sapply(newsp, mean)
## 		## correctionVals=sapply(newsp,function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2) & x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));})
## 		## names(correctionVals) <- names(newsp)
## 		fsampled <- sum(Ms[tokeep, i])
## 		##
## 		n <- length(tokeep)
## 		lambda <- fsampled/n
## 		#ngc <- sapply(newsp, length)
## 		##tvscore <- sum(ngc/n * (abs(correctionVals - lambda)))
## 		##if(verbose) print(paste('for array', i, 'tv score is', tvscore, 'when correction is applied to each chromosome'))
## 		##
## 		##allcorrections <- vector()
## 		##allcorrections[tokeep] <- correctionVals[match(paste(chr[tokeep], priorFracWremainingUse[tokeep], sep='.'), (names(correctionVals)))]
## 		##allcorrections[toremove] <- correctionVals[match(paste(chr[toremove], priorFracWremainingUse[toremove], sep='.'), (names(correctionVals)))]
## 		## this should retain information
## 		##
## 		if(verbose) message('we removed ', length(toremove), ' locations before calculating correction values')
## 		if(verbose) message(length(which(is.na(allcorrections[toremove]))), ' point(s) for array ', i, ' will recieve no gc correction')
## 		##allcorrections[toremove][which(is.na(allcorrections[toremove]))] <- lambda
## 		##
## 		mfact <- 1
## 		chromMediansAll <- aggregate(Ms[,i], list(chr), median)
## 		correctedM <- chromMediansAll$x[match(chr, chromMediansAll$Group.1)] + Ms[, i] - allcorrections*mfact
## 		##
## 		newTVscore <- vector()
## 		for(tvind in seq_len(nrow(tvScore))){
## 			priorFrac <- priorFracs(gcFracBoth, tvind, nparts, tvScore, increms)
## 			newTVscore[tvind] <- correctionTVscore(correctedM[chr %in% samplechr], priorFrac, i, as.numeric(rownames(tvScore)[tvind]))
## 		}
## 		names(newTVscore) <- rownames(tvScore)
## 		##
## 		if(length(which(is.na(correctedM))) > 0){
## 			if(verbose) print(priorFracWremaining[is.na(correctedM)][1:20])
## 		}
## 		##
## 		## comparison metrics
## 		if(verbose) print(paste('mad of corrected array', i, 'is', mad(correctedM)))
## 		if(verbose) print(paste('old mad array', i, 'is', mad(Ms[, i])))
## 		if(verbose) print(paste('The new Maximum tv Score for array', i, 'is', max(newTVscore), 'in window', names(newTVscore)[which.max(newTVscore)], 'bp'))
## 		##
## 		##
## 		## for(ii in unique(chr)){plot(Ms[chr==ii],ylim=c(-.8,.8));plot(correctedM[chr==ii],ylim=c(-.8,.8));readline('')}
## 		as(correctedM, "matrix")
## 	}
## 	return(correctedM)
## }
