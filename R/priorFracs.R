priorFracs <- function(gcFracBoth, maxuse, nparts, tvScore, increm, increm2){
	offset <- as.numeric(rownames(tvScore)[maxuse])/ifelse(maxuse>nrow(tvScore)/2,increm2,increm)
	if(maxuse > nrow(tvScore)/2)
		## Bigger Increment Windows are better
		priorGC <- gcFracBoth[, 2]
	else{
		## Smaller Increment Windows are better
		priorGC <- gcFracBoth[, 1]
	}
	## offset=maxuse
	priorFrac <- (cumsum(priorGC)[offset:length(priorGC)] - c(0, cumsum(priorGC)[1:(length(priorGC)-offset)]))[seq(1, (length(priorGC)-offset + 1), nparts)]/offset
	quants <- sort(unique(quantile(priorFrac, probs=seq(0, 1, .01))))
	priorFrac <- quants[findInterval(priorFrac, quants)]
  }
