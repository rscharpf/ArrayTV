\name{gcCorrectMain}
\alias{gcCorrectMain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcCorrectMain(Ms, chr, starts, nodes, increms, maxwins, outputTVs, outputChrReg, platform)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ms}{
%%     ~~Describe \code{Ms} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
}
  \item{increms}{
%%     ~~Describe \code{increms} here~~
}
  \item{maxwins}{
%%     ~~Describe \code{maxwins} here~~
}
  \item{outputTVs}{
%%     ~~Describe \code{outputTVs} here~~
}
  \item{outputChrReg}{
%%     ~~Describe \code{outputChrReg} here~~
}
  \item{platform}{
%%     ~~Describe \code{platform} here~~
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
function (Ms, chr, starts, nodes, increms, maxwins, outputTVs, 
    outputChrReg, platform) 
{
    library("BSgenome.Hsapiens.UCSC.hg18")
    library("RColorBrewer")
    library("multicore")
    library("foreach")
    library("grDevices")
    library("doMC")
    registerDoMC(nodes)
    if (exists("correctedM")) {
        rm(correctedM)
    }
    source("0605gcFragmentCorrectFunctions.R")
    if (outputTVs & outputChrReg) {
        pdf(paste("chrViews", platform, ".pdf", sep = ""))
    }
    samplechr = paste("chr", c(1, 4, 8, 10), sep = "")
    uniqchrs = unique(chr)
    remainingChr = uniqchrs[is.na(match(uniqchrs, samplechr))]
    narrays = ncol(Ms)
    for (correctind in 1:(length(increms))) {
        if (exists("correctedM")) {
            oldMs = Ms
            Ms = correctedM
        }
        increm = increms[correctind]
        maxwin = maxwins[correctind]
        nparts = (maxwin/increm) * 2
        strategyuse = 2
        gcFracBoth = gcFracAllWin(maxwin, increm, chr, starts, 
            samplechr, uniqchrs, strategyuse)
        useM = Ms[chr \%in\% samplechr, ]
        tvScore = calctvScore(gcFracBoth, samplechr, nparts, 
            useM, narrays)
        rownames(tvScore) = c(seq(increm, maxwin, increm), -seq(increm, 
            maxwin, increm))
        tvScore2 = calctvScore2(useM, gcFracBoth, tvScore, nparts, 
            narrays)
        maxTVscores = rep(0, narrays)
        correctionVals = list()
        optWins = rep(0, narrays)
        formaxs1 = tvScoreFirstHalfMax(tvScore)
        revmaxs1 = tvScoreSecondHalfMax(tvScore)
        formaxs2 = tvScoreSecondHalfMax(tvScore2)
        revmaxs2 = tvScoreFirstHalfMax(tvScore2)
        gmaxvals = apply(tvScore, 2, max)
        gmaxvals2 = apply(tvScore2, 2, max)
        maxorigins = gmaxvals >= gmaxvals2
        if (outputTVs) {
            plotTV(tvScore, tvScore2, narrays, increm)
        }
        correctedM = CorrectM(gcFracBoth, useM, Ms, starts, narrays, 
            nparts, chr, samplechr, remainingChr, increm, formaxs1, 
            formaxs2, revmaxs1, revmaxs2, maxorigins, gmaxvals, 
            gmaxvals2)
        if (outputChrReg) {
            if (exists("oldMs")) {
                Ms = oldMs
            }
            for (jj in 1:narrays) {
                par(mfrow = c(2, 1))
                minrange = 8.5e+07
                maxrange = 8.9e+07
                Msplot = Ms[starts > minrange & starts < maxrange & 
                  chr == "chr1", jj]
                ylim1 = quantile(Msplot, c(0.002, 0.998))
                plot(starts[starts > minrange & starts < maxrange & 
                  chr == "chr1"], Msplot, ylim = ylim1)
                plot(starts[starts > minrange & starts < maxrange & 
                  chr == "chr1"], correctedM[starts > minrange & 
                  starts < maxrange & chr == "chr1", jj], ylim = ylim1)
            }
        }
    }
    graphics.off()
    return(correctedM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
