# , maxwins='numeric', increms='numeric', chr='character', build='character'
setMethod("computeGC", signature(x = "numeric"), function(x, maxwins, increms, chr, 
    build, ...) {
    annotation.pkg <- paste("BSgenome.Hsapiens.UCSC.", build, sep = "")
    gcFracBoth <- gcFracAllWin(maxwins = maxwins, increms = increms, chr = chr, allstarts = x, 
        samplechr = unique(chr), uniqchrs = unique(chr), strategyuse = 2, annotation.pkg = annotation.pkg, 
        verbose = TRUE)
})


setMethod("computeGC", signature(x = "BafLrrSetList"), function(x, maxwins, increms, 
    chr, build, ...) {
    chr <- paste("chr", rep(chromosome(x), elementLengths(x)), sep = "")
    pos <- unlist(position(x))
    build <- genomeBuild(x)
    annotation.pkg <- paste("BSgenome.Hsapiens.UCSC.", build, sep = "")
    gcFracBoth <- gcFracAllWin(maxwins = maxwins, increms = increms, chr = chr, allstarts = pos, 
        samplechr = unique(chr), uniqchrs = unique(chr), strategyuse = 2, annotation.pkg = annotation.pkg, 
        verbose = TRUE)
    return(gcFracBoth)
}) 
