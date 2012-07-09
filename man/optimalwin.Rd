\name{optimalwin}
\alias{optimalwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimalwin(useM, gcFracBoth, tvScore, vv, ar, nparts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{useM}{
%%     ~~Describe \code{useM} here~~
}
  \item{gcFracBoth}{
%%     ~~Describe \code{gcFracBoth} here~~
}
  \item{tvScore}{
%%     ~~Describe \code{tvScore} here~~
}
  \item{vv}{
%%     ~~Describe \code{vv} here~~
}
  \item{ar}{
%%     ~~Describe \code{ar} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
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
function (useM, gcFracBoth, tvScore, vv, ar, nparts) 
{
    optnum = 4
    v1 = c(3, 1)
    v2 = c(4, 2)
    tvScoreForwards = tvScore[1:(nrow(tvScore)/2), ]
    tvScoreReverses = tvScore[((nrow(tvScore)/2) + 1):nrow(tvScore), 
        ]
    formaxs = apply(tvScoreForwards, 2, which.max)
    revmaxs = apply(tvScoreReverses, 2, which.max)
    tvScorepart = vector()
    maxuse = ifelse(vv == 1, formaxs[ar], revmaxs[ar])
    priorFrac = priorFracs(gcFracBoth, maxuse, vv, nparts, v1, 
        v2)
    priorWeight = maxuse
    for (hh in 1:((nparts/optnum) * 2)) {
        if (hh > (nparts/optnum)) {
            partitionedgcFrac = gcFracBoth[[v2[vv]]]
        }
        else {
            partitionedgcFrac = gcFracBoth[[v1[vv]]]
        }
        cat(".")
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
        newsp = split(useM[, ar], round(gcFracA, 2))
        toplot = sapply(newsp, mean)
        names(toplot) = names(newsp)
        fsampled = sum(useM[, ar])
        n = nrow(useM)
        lambda = fsampled/n
        ngc = sapply(newsp, length)
        tvscore = sum(ngc/n * (abs(toplot - lambda)))
        tvScorepart = c(tvScorepart, tvscore)
    }
    tvScorepart
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
