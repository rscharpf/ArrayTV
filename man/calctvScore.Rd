\name{calctvScore}
\alias{calctvScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tvScore
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calctvScore(gcFracBoth, samplechr, nparts, useM, narrays)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gcFracBoth}{
%%     ~~Describe \code{gcFracBoth} here~~
}
  \item{samplechr}{
%%     ~~Describe \code{samplechr} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
}
  \item{useM}{
%%     ~~Describe \code{useM} here~~
}
  \item{narrays}{
%%     ~~Describe \code{narrays} here~~
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
function (gcFracBoth, samplechr, nparts, useM, narrays) 
{
    optnum = 4
    narrays = ncol(useM)
    tvScore = foreach(vv = 1:4, .combine = rbind) \%dopar\% {
        tvScorepart = matrix(0, nrow = nparts/optnum, ncol = narrays)
        partitionedgcFrac = gcFracBoth[[vv]]
        if (vv \%in\% c(1, 3)) {
            priorFrac = 0
            priorWeight = 0
        }
        else if (vv \%in\% c(2, 4)) {
            priorGC = gcFracBoth[[vv - 1]]
            priorFrac = (cumsum(priorGC)[(nparts/optnum):length(priorGC)] - 
                c(0, cumsum(priorGC)[1:(length(priorGC) - (nparts/optnum))]))[seq(1, 
                (length(priorGC) - (nparts/optnum) + 1), nparts/optnum)]/(nparts/optnum)
            priorWeight = nparts/optnum
        }
        for (hh in 1:(nparts/optnum)) {
            if (hh == 1) {
                gcFracA = (partitionedgcFrac[seq(hh, length(partitionedgcFrac), 
                  nparts/optnum)] + (priorFrac * priorWeight))/(priorWeight + 
                  1)
            }
            else {
                gcFracA = (gcFracA * ((priorWeight + hh - 1)/(priorWeight + 
                  hh))) + partitionedgcFrac[seq(hh, length(partitionedgcFrac), 
                  nparts/optnum)]/(priorWeight + hh)
            }
            newsp = split(useM, round(gcFracA, 2))
            toplot = sapply(newsp, function(x) {
                cc = matrix(x, ncol = narrays)
                y = colMeans(cc)
                y
            })
            rownames(toplot) = colnames(useM)
            colnames(toplot) = names(newsp)
            fsampled = colSums(useM)
            n = nrow(useM)
            lambda = fsampled/n
            ngc = sapply(newsp, length)/narrays
            tvscore = rowSums(rep((ngc/n), each = narrays) * 
                (abs(toplot - lambda)))
            tvScorepart[hh, ] = na.omit(tvscore)
        }
        tvScorepart
    }
    return(tvScore)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
