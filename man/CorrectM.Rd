\name{CorrectM}
\alias{CorrectM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CorrectM
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrectM(gcFracBoth, useM, Ms, starts, narrays, nparts, chr, samplechr, remainingChr, increm, increm2, gmaxvals, gmaxvalsInd, tvScore)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gcFracBoth}{
%%     ~~Describe \code{gcFracBoth} here~~
}
  \item{useM}{
%%     ~~Describe \code{useM} here~~
}
  \item{Ms}{
%%     ~~Describe \code{Ms} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{narrays}{
%%     ~~Describe \code{narrays} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{samplechr}{
%%     ~~Describe \code{samplechr} here~~
}
  \item{remainingChr}{
%%     ~~Describe \code{remainingChr} here~~
}
  \item{increm}{
%%     ~~Describe \code{increm} here~~
}
  \item{increm2}{
%%     ~~Describe \code{increm2} here~~
}
  \item{gmaxvals}{
%%     ~~Describe \code{gmaxvals} here~~
}
  \item{gmaxvalsInd}{
%%     ~~Describe \code{gmaxvalsInd} here~~
}
  \item{tvScore}{
%%     ~~Describe \code{tvScore} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gcFracBoth, useM, Ms, starts, narrays, nparts, chr, 
    samplechr, remainingChr, increm, increm2, gmaxvals, gmaxvalsInd, 
    tvScore) 
{
    priorFracWremaining = foreach(i = which(!duplicated(gmaxvalsInd)), 
        .combine = "list", .multicombine = T) \%dopar\% {
        print(paste("A maximum first pass  TV is", gmaxvals[i], 
            "in window", rownames(tvScore)[gmaxvalsInd[i]]))
        maxuse1 = gmaxvalsInd[i]
        priorFrac = priorFracs(gcFracBoth, maxuse1, nparts, tvScore, 
            increm, increm2)
        if (length(remainingChr) == 0) {
            priorFracWremaining = priorFracUse
        }
        else {
            forwardExtend = as.numeric(rownames(tvScore)[maxuse1])
            reverseExtend = forwardExtend
            priorFracWremaining = priorFracsRestOfGenome(forwardExtend, 
                reverseExtend, remainingChr, starts, chr)
            print(paste("forward extend", forwardExtend, "reverse extend", 
                reverseExtend))
            priorFracWremaining[chr \%in\% samplechr] = priorFrac
        }
        priorFracWremaining
    }
    correctedM = foreach(i = 1:narrays, .combine = "cbind") \%dopar\% 
        {
            if (narrays == 1) {
                priorFracWremainingUse = priorFracWremaining
            }
            else {
                priorFracWremainingUse = priorFracWremaining[[match(gmaxvalsInd[i], 
                  unique(gmaxvalsInd))]]
            }
            cm = cumsum(Ms[, i])
            cm2 = cm[200:length(cm)] - c(0, cm[1:(length(cm) - 
                200)])
            inds = seq(1, length(cm2), 200)
            cm3 = cm2[inds]
            chrForCNA = chr[inds]
            startsForCNA = starts[inds]
            CNA.object = CNA(cm3, chrForCNA, startsForCNA, data.type = c("logratio"), 
                sampleid = i)
            segment.smoothed.CNA.object <- segment(CNA.object, 
                verbose = 1, min.width = 5)
            markrep = rep(segment.smoothed.CNA.object$output$seg.mean, 
                segment.smoothed.CNA.object$output$num.mark)
            chromrep = rep(segment.smoothed.CNA.object$output$chrom, 
                segment.smoothed.CNA.object$output$num.mark)
            meduse = mean(markrep)
            maduse = mad(markrep)
            blockchrom = table(segment.smoothed.CNA.object$output$chrom)
            multiblock = names(blockchrom)[blockchrom > 1]
            chromMeans = aggregate(markrep, list(chromrep), mean)
            segment.smoothed.CNA.object$output$chromMean = chromMeans$x[match(segment.smoothed.CNA.object$output$chrom, 
                chromMeans$Group.1)]
            toremove = which(abs(segment.smoothed.CNA.object$output$seg.mean - 
                segment.smoothed.CNA.object$output$chromMean) > 
                maduse & segment.smoothed.CNA.object$output$chrom \%in\% 
                multiblock & (segment.smoothed.CNA.object$output$seg.mean > 
                meduse + (maduse * 2.5) | segment.smoothed.CNA.object$output$seg.mean < 
                meduse - (maduse * 2.5)))
            locsAsRange = GRanges(seqnames = chr, ranges = IRanges(start = starts, 
                width = 1))
            removestart = segment.smoothed.CNA.object$output$loc.start[toremove]
            removeend = segment.smoothed.CNA.object$output$loc.end[toremove]
            removechr = segment.smoothed.CNA.object$output$chrom[toremove]
            removeAsRange = GRanges(seqnames = removechr, ranges = IRanges(start = removestart, 
                end = removeend))
            tokeepLogical = is.na(match(locsAsRange, removeAsRange))
            tokeep = which(tokeepLogical)
            toremove = which(!tokeepLogical)
            newsp = split(Ms[tokeep, i], paste(chr[tokeep], priorFracWremainingUse[tokeep], 
                sep = "."))
            correctionVals = sapply(newsp, mean)
            names(correctionVals) = names(newsp)
            fsampled = sum(Ms[tokeep, i])
            n = length(tokeep)
            lambda = fsampled/n
            ngc = sapply(newsp, length)
            tvscore = sum(ngc/n * (abs(correctionVals - lambda)))
            print(paste("for array", i, "tv score is", tvscore, 
                "when correction is applied to each chromosome"))
            allcorrections = vector()
            allcorrections[tokeep] = correctionVals[match(paste(chr[tokeep], 
                priorFracWremainingUse[tokeep], sep = "."), (names(correctionVals)))]
            allcorrections[toremove] = correctionVals[match(paste(chr[toremove], 
                priorFracWremainingUse[toremove], sep = "."), 
                (names(correctionVals)))]
            print(paste("we removed", length(toremove), "locations before calculating correction values"))
            print(paste(length(which(is.na(allcorrections[toremove]))), 
                "point(s) for array", i, "will recieve no gc correction"))
            allcorrections[toremove][which(is.na(allcorrections[toremove]))] = lambda
            cnr = 1
            asum = vector()
            for (jj in seq(1, 1.2, 0.05)) {
                autocor = acf(Ms[chr == samplechr[1], i] - allcorrections[chr == 
                  samplechr[1]] * jj, lag.max = 200, plot = F)
                asum[cnr] = sum(autocor$acf[2:200])
                print(paste("autocorrelation sum for correction factor", 
                  jj, "is", asum[cnr]))
                cnr = cnr + 1
            }
            mfact = seq(1, 1.2, 0.05)[which.min(asum)]
            chromMediansAll = aggregate(Ms[, i], list(chr), median)
            correctedM = chromMediansAll$x[match(chr, chromMediansAll$Group.1)] + 
                Ms[, i] - allcorrections * mfact
            newTVscore = vector()
            for (tvind in 1:nrow(tvScore)) {
                priorFrac = priorFracs(gcFracBoth, tvind, nparts, 
                  tvScore, increm, increm2)
                newTVscore[tvind] = correctionTVscore(correctedM[chr \%in\% 
                  samplechr], priorFrac, i, as.numeric(rownames(tvScore)[tvind]))
            }
            names(newTVscore) = rownames(tvScore)
            if (length(which(is.na(correctedM))) > 0) {
                print(priorFracWremaining[is.na(correctedM)][1:20])
            }
            print(paste("mad of corrected array", i, "is", mad(correctedM)))
            print(paste("old mad array", i, "is", mad(Ms[, i])))
            print(paste("The new Maximum tv Score for array", 
                i, "is", max(newTVscore), "in window", names(newTVscore)[which.max(newTVscore)], 
                "bp"))
            as(correctedM, "matrix")
        }
    return(correctedM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
