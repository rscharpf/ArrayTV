priorFracs <- function(gcFracBoth, maxuse, nparts, tvScore, increms) {
    incremind <- rep(1:length(increms), each = nparts)[maxuse]
    offset <- as.numeric(rownames(tvScore)[maxuse])/increms[incremind]
    priorGC <- gcFracBoth[, incremind]
    ## offset=maxuse
    priorFrac <- (cumsum(priorGC)[offset:length(priorGC)] - c(0, cumsum(priorGC)[1:(length(priorGC) - 
        offset)]))[seq(1, (length(priorGC) - offset + 1), nparts)]/offset
    quants <- sort(unique(quantile(priorFrac, probs = seq(0, 1, 0.01))))
    priorFrac <- quants[findInterval(priorFrac, quants)]
} 
