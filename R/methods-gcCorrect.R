setMethod("gcCorrect", signature(object="matrix"), function(object, ...){
    gcCorrectMain(object, ...)
})

setMethod("gcCorrect", signature(object="BeadStudioSet"),
	  function(object, ...){



	  })

setMethod("gcCorrect", signature(object="BeadStudioSetList"),
	  function(object, ...){
		  args <- list(...)
		  if("returnOnlyTV" %in% args){
			  is.score <- returnOnlyTV
		  }
		  r <- lrr(object)
		  index.samples <- seq_len(ncol(object))
		  ## to keep RAM in check, do in batches of samples
		  index.list <- splitIndicesByLength(index.samples, ocSamples())
		  if(is.score) score.list <- list()
		  for(i in seq_along(index.list)){
			  j <- index.list[[i]]
			  rr <- lapply(r, function(x, j) x[, j, drop=FALSE]/100, j=j)
			  R <- do.call(rbind, rr)
			  R[is.na(R)] <- 0 ## not ideal
			  l <- sapply(object, nrow)
			  chr <- paste("chr", rep(chromosome(object), l), sep="")
			  pos <- unlist(position(object))
			  res <- gcCorrectMain(Ms=R, ##R[, 1, drop=FALSE],
					       chr=chr,
					       starts=pos,
					       samplechr=unique(chr),
					       build=genomeBuild(object),
					       ...)
			  if(!is.score){  ## if not returning TV score, update the brList object
				  res <- integerMatrix(res, 100)
				  for(k in seq_along(object)){ ## for each chromosome
					  bset <- object[[k]]
					  ii <- match(featureNames(bset), rownames(cM))
					  lrr(bset)[, j] <- cM[ii, , drop=FALSE]
					  object[[k]] <- bset
				  }
			  } else {
				  score.list[[i]] <- res
			  }
		  }
		  if(!is.score) {
			  results <- score.list
		  } else results <- object
		  return(results)
	  })
