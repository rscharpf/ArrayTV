calctvScore <- function(gcFracBoth, nparts, useM, narrays, increms, verbose = FALSE) {
    narrays = ncol(useM)
    ## the order of row binding is forward first half, forward second half, reverse
    ## first half, reverse second half
    i <- NULL
    tvScorepart = foreach(i = 1:ncol(gcFracBoth), .combine = "rbind", .packages = "ArrayTV") %dopar% 
        {
            partitionedgcFrac <- gcFracBoth[, i]
            tvScorepart <- matrix(0, nrow = nparts, ncol = narrays)
            ## maxwin/increm * 2 must be divisible by 4 right now
            for (hh in 1:nparts) {
                if (hh == 1) {
                  gcFracA <- partitionedgcFrac[seq(hh, length(partitionedgcFrac), 
                    nparts)]
                  ## added
                  ## gcFracA2=priorpartitionedgcFrac[seq(hh,length(priorpartitionedgcFrac),nparts)]
                  ## gcFracA=(gcFracA1+gcFracA2)/2 end additon
                } else {
                  gcFracA <- (gcFracA * ((hh - 1)/(hh))) + partitionedgcFrac[seq(hh, 
                    length(partitionedgcFrac), nparts)]/(hh)
                  ## added ##
                  ## gcFracA2=(gcFracA2*((hh-1)/(hh)))+priorpartitionedgcFrac[seq(hh,length(priorpartitionedgcFrac),nparts)]/(hh)
                  ## gcFracA=(gcFracA1+gcFracA2)/2
                }
                quants <- sort(unique(quantile(gcFracA, probs = seq(0, 1, 0.01))))
                gcFracA <- quants[findInterval(gcFracA, quants)]
                newsp <- split(useM, gcFracA)
                ## sapply(newsp,function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2) &
                ## x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));})
                ## toplot=sapply(newsp,function(x){cc=matrix(x,ncol=narrays);y=apply(cc,2,function(x){x1=median(x);x2=1.5*mad(x);z=x[x<(x1+x2)
                ## & x>(x1-x2)];ifelse(length(z)>0,mean(z),mean(x));});y})
                toplot <- sapply(newsp, function(x) {
                  x = matrix(x, ncol = narrays)
                  y = colMeans(x)
                  y
                })
                ## toplot=sapply(newsp,function(x){cc=matrix(x,ncol=narrays);y=apply(cc,2,median);y})
                fsampled <- colSums(useM)
                n <- nrow(useM)
                lambda <- fsampled/n
                ngc <- sapply(newsp, length)/narrays
                if (is.null(nrow(toplot))) {
                  tvscore <- sum(rep((ngc/n), each = narrays) * (abs(toplot - lambda)))
                } else {
                  tvscore <- rowSums(rep((ngc/n), each = narrays) * (abs(toplot - 
                    lambda)))
                }
                if (verbose) 
                  message(paste("windows stretch", increms[i] * hh, "bp from probes, TV:", 
                    round(tvscore, 4)))
                tvScorepart[hh, ] <- na.omit(tvscore)
            }
            tvScorepart
        }
    return(tvScorepart)
} 
