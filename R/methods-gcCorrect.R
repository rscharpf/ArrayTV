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



gcCorrectBafLrrList <- function(object, index.samples, providedGC=NULL,...){
	args <- list(...)
	if("returnOnlyTV" %in% args){
		return.score <- returnOnlyTV
	} else return.score <- FALSE
	if(return.score) stop(paste("return.score not implemented for ", class(object), " objects"))
	r <- lrr(object)
	isff <- is(r[[1]], "ff")
	if(missing(index.samples))
		index.samples <- seq_len(ncol(object[[1]]))
	## to keep RAM in check, do in batches of samples
	index.list <- splitIndicesByLength(index.samples, ocSamples())
	l <- elementLengths(object)
	chr <- paste("chr", rep(chromosome(object), l), sep="")
	pos <- unlist(position(object))
	##if(return.score) score.list <- list()
	.packages <- c("oligoClasses", "ArrayTV")
	isFFloaded <- isPackageLoaded("ff")
	if(isFFloaded) .packages <- c("ff", .packages)
	gc.provided <- !is.null(providedGC)
	if(gc.provided) providedGC <- as.matrix(providedGC)
	if("returnOnlyTV" %in% names(list(...))){
		only.tv <- list(...)[["returnOnlyTV"]]
	} else only.tv <- FALSE
	reslist <- foreach(j=index.list, .packages=.packages) %dopar% {
		rr <- lapply(r, function(x, j) x[, j, drop=FALSE]/100, j=j)
		R <- do.call(rbind, rr)
		rm(rr)
		namatrix <- is.na(R)
		R[namatrix] <- 0 ## not ideal
		if(gc.provided){
			for(k in seq_len(ncol(providedGC))){
				R <- gcCorrectMain(Ms=R,
						   chr=chr,
						   starts=pos,
						   samplechr=unique(chr),
						   build=genomeBuild(object),
						   providedGC=providedGC[, k],
						   ...)
			}
		} else {
			R <- gcCorrectMain(Ms=R,
					   chr=chr,
					   starts=pos,
					   samplechr=unique(chr),
					   build=genomeBuild(object),
					   providedGC=providedGC,
					   ...)
		}
		if(!only.tv) {
			R <- integerMatrix(R, 100)
			R[namatrix] <- NA
		}
		##rm(R); gc()
		##
		## update ff object.  Writing is expensive -- only do
		## once
		##
		## Only reasonable to do this when ff package is
		## loaded. Otherwise, the replacment method is
		## transient and we do not want to return an entire
		## copy of the object
		if(isFFloaded & !only.tv) {
			lrr(object) <- R
			R <- NULL
		}
		gc()
		return(R)
	}
	if(!only.tv){
		if(!isFFloaded){
			res <- do.call("cbind", reslist)
			lrr(object) <- res
		} else {
			## we do not want to return an entire copy of the
			## object if ff package is loaded
			object <- NULL
		}
	} else object <- reslist
	return(object)
}

setMethod("gcCorrect", signature(object="BafLrrSetList"),
	  function(object, ...){
		  gcCorrectBafLrrList(object, ...)
	  })
