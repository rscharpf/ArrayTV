setMethod("gcCorrect", signature(object="matrix"), function(object, ...){
    gcCorrectMain(object, ...)
})

gcCorrectBeadStudioSet <- function(object, ...){
	args <- list(...)
	if("returnOnlyTV" %in% args){
		is.score <- returnOnlyTV
	} else is.score <- FALSE
	r <- lrr(object)
	if(is(r, "matrix")){
		r <- r/100
	}
	isff <- is(r, "ff_matrix") ## need to be careful
	if(isff){
		index.list <- split(seq_len(ncol(object)), ocSamples())
	} else index.list <- list(seq_len(ncol(object)))
	for(i in seq_along(index.list)){
		j <- index.list[[i]]
		pos <- position(object)
		chr <- paste("chr", chromosome(object), sep="")
		res <- gcCorrectMain(Ms=r,
				     chr=chr,
				     start=pos,
				     samplechr=unique(chr),
				     build=genomeBuild(object),
				     ...)
		if(!is.score){  ## if not returning TV score, update the brList object
			res <- integerMatrix(res, 100)
			lrr(object) <- res
		} else {
			score.list[[i]] <- res
		}
	}
	if(is.score) {
		results <- score.list
	} else results <- object
	return(results)
}

setMethod("gcCorrect", signature(object="BeadStudioSet"),
	  function(object, ...){
		  gcCorrectBeadStudioSet(object, ...)
	  })
setMethod("gcCorrect", signature(object="BafLrrSet"),
	  function(object, ...){
		  gcCorrectBeadStudioSet(object, ...)
	  })

gcCorrectBafLrrList <- function(object, ...){
	args <- list(...)
	if("returnOnlyTV" %in% args){
		return.score <- returnOnlyTV
	} else return.score <- FALSE
	r <- lrr(object)
	isff <- is.ff(r[[1]])
	index.samples <- seq_len(ncol(object[[1]]))
	## to keep RAM in check, do in batches of samples
	index.list <- splitIndicesByLength(index.samples, ocSamples())
	if(return.score) score.list <- list()
	for(i in seq_along(index.list)){
		j <- index.list[[i]]
		rr <- lapply(r, function(x, j) x[, j, drop=FALSE]/100, j=j)
		R <- do.call(rbind, rr)
		R[is.na(R)] <- 0 ## not ideal
		fdl <- featureData(object)
		l <- sapply(fdl, nrow)
		chr <- paste("chr", rep(chromosome(object), l), sep="")
		pos <- unlist(position(object))
		res <- gcCorrectMain(Ms=as.matrix(R),
				     chr=chr,
				     starts=pos,
				     samplechr=unique(chr),
				     build=genomeBuild(object),
				     ...)
		fns <- unlist(sapply(fdl, featureNames))
		dimnames(res) <- list(fns, sampleNames(object)[j])
		if(!return.score){  ## if not returning TV score, update the brList object
			res <- integerMatrix(res, 100)
			if(!isff){
				## do not use lrr(object)[,j] -- the replacement method takes care of this
				lrr(object) <- as.matrix(res)
			} else {
				## clone the matrix and then write to file
			}
		} else {
			score.list[[i]] <- res
		}
	}
	if(return.score) {
		results <- score.list
	} else results <- object
	return(results)
}

setMethod("gcCorrect", signature(object="BafLrrSetList"),
	  function(object, ...){
		  gcCorrectBafLrrList(object, ...)
	  })
