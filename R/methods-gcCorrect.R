setMethod("gcCorrect", signature(object = "matrix"), function(object, ...) {
    gcCorrectMain(object, ...)
})

gcCorrectBeadStudioSet <- function(object, ...) {
    args <- list(...)
    if ("returnOnlyTV" %in% names(args)) {
        is.score <- args[["returnOnlyTV"]]
    } else is.score <- FALSE
    r <- lrr(object)
    if (is(r, "matrix")) {
        r <- r/100
    }
    isff <- is(r, "ff_matrix")  ## need to be careful
    if (isff) {
        index.list <- split(seq_len(ncol(object)), ocSamples())
    } else index.list <- list(seq_len(ncol(object)))
    if (is.score) 
        score.list <- list()
    for (i in seq_along(index.list)) {
        j <- index.list[[i]]
        pos <- position(object)
        chr <- paste("chr", chromosome(object), sep = "")
        res <- gcCorrectMain(Ms = r, chr = chr, starts = pos, samplechr = unique(chr), 
            build = genomeBuild(object), ...)
        if (!is.score) {
            ## if not returning TV score, update the brList object
            res <- integerMatrix(res, 100)
            lrr(object) <- res
        } else {
            score.list[[i]] <- res
        }
    }
    if (is.score) {
        results <- score.list
    } else results <- object
    return(results)
}

setMethod("gcCorrect", signature(object = "BeadStudioSet"), function(object, ...) {
    gcCorrectBeadStudioSet(object, ...)
})
setMethod("gcCorrect", signature(object = "BafLrrSet"), function(object, ...) {
    gcCorrectBeadStudioSet(object, ...)
})



gcCorrectBafLrrList <- function(object, index.samples, providedGC = NULL, ...) {
    args <- list(...)
    if ("returnOnlyTV" %in% names(args)) {
        return.score <- args[["returnOnlyTV"]]
    } else return.score <- FALSE
    if (return.score) 
        stop(paste("return.score not implemented for ", class(object), " objects"))
    r <- lrr(object)
    isff <- is(r[[1]], "ff")
    if (missing(index.samples)) 
        index.samples <- seq_len(ncol(object[[1]]))
    ## to keep RAM in check, do in batches of samples
    index.list <- splitIndicesByLength(index.samples, ocSamples())
    l <- elementLengths(object)
    chr <- paste("chr", rep(chromosome(object), l), sep = "")
    pos <- unlist(position(object))
    ## if(return.score) score.list <- list()
    .packages <- c("oligoClasses", "ArrayTV")
    isFFloaded <- isPackageLoaded("ff")
    if (isFFloaded) 
        .packages <- c("ff", .packages)
    gc.provided <- !is.null(providedGC)
    if (gc.provided) 
        providedGC <- as.matrix(providedGC)
    if ("returnOnlyTV" %in% names(list(...))) {
        only.tv <- list(...)[["returnOnlyTV"]]
    } else only.tv <- FALSE
    j <- 1
    reslist <- foreach(j = index.list, .packages = .packages) %dopar% {
        rr <- lapply(r, function(x, j) x[, j, drop = FALSE]/100, j = j)
        R <- do.call(rbind, rr)
        rm(rr)
        namatrix <- is.na(R)
        R[namatrix] <- 0  ## not ideal
        if (gc.provided) {
            for (k in seq_len(ncol(providedGC))) {
                R <- (gcCorrectMain(Ms = R, chr = chr, starts = pos, samplechr = unique(chr), 
                  build = genomeBuild(object), providedGC = providedGC[, k], ...))[["correctedVals"]]
            }
        } else {
            R <- (gcCorrectMain(Ms = R, chr = chr, starts = pos, samplechr = unique(chr), 
                build = genomeBuild(object), providedGC = providedGC, ...))[["correctedVals"]]
        }
        if (!only.tv) {
            R <- integerMatrix(R, 100)
            R[namatrix] <- NA
        }
        ## rm(R); gc() update ff object.  Writing is expensive -- only do once Only
        ## reasonable to do this when ff package is loaded. Otherwise, the replacment
        ## method is transient and we do not want to return an entire copy of the object
        if (isFFloaded & !only.tv) {
            lrr(object) <- R
            R <- NULL
        }
        gc()
        return(R)
    }
    if (!only.tv) {
        if (!isFFloaded) {
            res <- do.call("cbind", reslist)
            lrr(object) <- res
        } else {
            ## we do not want to return an entire copy of the object if ff package is loaded
            object <- NULL
        }
    } else object <- reslist
    return(object)
}

setMethod("gcCorrect", signature(object = "BafLrrSetList"), function(object, ...) {
    gcCorrectBafLrrList(object, ...)
})

gcModel <- function(data, window, verbose = FALSE) {
    ## assume data is summarized experiment
    build <- metadata(rowData(data))$genome
    library(paste("BSgenome.Hsapiens.UCSC.", build, sep = ""), character.only = TRUE)
    if (verbose) 
        message("Getting gc content From BS genome Object")
    Hsapiens <- get("Hsapiens")
    chroms <- unique(chromosome(data))
    maxwin <- increm <- window
    ## query genome prior to looking at genomic position of markers??
    gc <- list()
    for (i in seq_along(chroms)) {
        if (verbose) 
            cat(".")
        chr <- chroms[i]
        j <- which(chromosome(data) == chr)
        ## calculates the gc content for each window of size increm
        pregcFrac <- letterFrequencyInSlidingView(unmasked(Hsapiens[[chr]]), view.width = increm, 
            "CG", as.prob = TRUE)
        ## if(verbose) print('gc content stored')
        startinds <- rep(start(data)[j], each = maxwin/increm) + rep(seq(0, maxwin - 
            increm, increm), length(j))
        startinds[which(startinds < 1)] <- 1
        startinds[which(startinds > length(pregcFrac))] <- length(pregcFrac)
        startindbackwards <- rep(start(data)[j], each = maxwin/increm) - rep(seq(increm, 
            maxwin, increm), length(j))
        gcFrac1 <- pregcFrac[startinds]
        if (any(is.na(gcFrac1))) 
            gcFrac1[is.na(gcFrac1)] <- mean(gcFrac1, na.rm = TRUE)
        startindbackwards[which(startindbackwards < 1)] <- 1
        startindbackwards[which(startindbackwards > length(pregcFrac))] <- length(pregcFrac)
        gcFracbackwards1 <- pregcFrac[startindbackwards]
        if (any(is.na(gcFracbackwards1))) {
            gcFracbackwards1[is.na(gcFracbackwards1)] <- mean(gcFracbackwards1, na.rm = TRUE)
        }
        gc[[i]] <- (gcFrac1 + gcFracbackwards1)/2
    }
    result <- unlist(gc)
    return(result)
}

midpoint <- function(object) start(object) + floor((width(object) - 1)/2)

.rescaleGC <- function(gc, cn, shift = 0) {
    x <- scale(gc)  ## mean 0, sd1
    ## give it same sd as cn
    x <- x * mad(cn, na.rm = TRUE)
    ## give it same location as cn
    x <- x + mean(cn, na.rm = TRUE) + shift
    x
}

gcModelSeq <- function(data, window, verbose = FALSE) {
    ## assume data is summarized experiment
    build <- genome(rowData(data))[[1]]
    library(paste("BSgenome.Hsapiens.UCSC.", build, sep = ""), character.only = TRUE)
    if (verbose) 
        message("Getting gc content From BS genome Object")
    Hsapiens <- get("Hsapiens")
    chroms <- unique(chromosome(data))
    maxwin <- increm <- window
    starts <- midpoint(data)
    ## query genome prior to looking at genomic position of markers??
    gc <- list()
    for (i in seq_along(chroms)) {
        if (verbose) 
            cat(".")
        chr <- chroms[i]
        j <- which(chromosome(data) == chr)
        ## calculates the gc content for each window of size increm
        pregcFrac <- letterFrequencyInSlidingView(unmasked(Hsapiens[[chr]]), view.width = increm, 
            "CG", as.prob = TRUE)
        ## if(verbose) print('gc content stored')
        startinds <- rep(starts, each = maxwin/increm) + rep(seq(0, maxwin - increm, 
            increm), length(j))
        startinds[which(startinds < 1)] <- 1
        startinds[which(startinds > length(pregcFrac))] <- length(pregcFrac)
        startindbackwards <- rep(starts, each = maxwin/increm) - rep(seq(increm, 
            maxwin, increm), length(j))
        gcFrac1 <- pregcFrac[startinds]
        if (any(is.na(gcFrac1))) 
            gcFrac1[is.na(gcFrac1)] <- mean(gcFrac1, na.rm = TRUE)
        startindbackwards[which(startindbackwards < 1)] <- 1
        startindbackwards[which(startindbackwards > length(pregcFrac))] <- length(pregcFrac)
        gcFracbackwards1 <- pregcFrac[startindbackwards]
        if (any(is.na(gcFracbackwards1))) {
            gcFracbackwards1[is.na(gcFracbackwards1)] <- mean(gcFracbackwards1, na.rm = TRUE)
        }
        gc[[i]] <- (gcFrac1 + gcFracbackwards1)/2
    }
    result <- unlist(gc)
    return(result)
} 
