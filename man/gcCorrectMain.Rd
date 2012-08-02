\name{gcCorrectMain}
\alias{gcCorrectMain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gcCorrectMain
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcCorrectMain(Ms, chr, starts, samplechr, nodes, increms, maxwins, build)
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
  \item{samplechr}{
%%     ~~Describe \code{samplechr} here~~
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
  \item{build}{
%%     ~~Describe \code{build} here~~
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
function (Ms, chr, starts, samplechr, nodes, increms, maxwins, 
    build) 
{
    do.call(library, list(paste("BSgenome.Hsapiens.UCSC.", build, 
        sep = "")))
    library("RColorBrewer")
    library("multicore")
    library("foreach")
    library("grDevices")
    library("doMC")
    library("DNAcopy")
    registerDoMC(nodes)
    uniqchrs = unique(chr)
    remainingChr = uniqchrs[is.na(match(uniqchrs, samplechr))]
    narrays = ncol(Ms)
    increm = increms[1]
    increm2 = increms[2]
    maxwin = maxwins[1]
    maxwin2 = maxwins[2]
    nparts = (maxwin/increm)
    if (nparts != maxwin2/increm2) {
        readline("Bad Increment Values")
    }
    strategyuse = 2
    gcFracBoth = gcFracAllWin(maxwin, maxwin2, increm, increm2, 
        chr, starts, samplechr, uniqchrs, strategyuse)
    useM = as(Ms[chr \%in\% samplechr, ], "matrix")
    tvScore = calctvScore(gcFracBoth, samplechr, nparts, useM, 
        narrays, increm, increm2)
    rownames(tvScore) = c(seq(increm, maxwin, increm), seq(increm2, 
        maxwin2, increm2))
    maxTVscores = rep(0, narrays)
    correctionVals = list()
    optWins = rep(0, narrays)
    gmaxvals = apply(tvScore, 2, max)
    gmaxvalsInd = apply(tvScore, 2, which.max)
    correctedM = CorrectM(gcFracBoth, useM, Ms, starts, narrays, 
        nparts, chr, samplechr, remainingChr, increm, increm2, 
        gmaxvals, gmaxvalsInd, tvScore)
    return(correctedM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
