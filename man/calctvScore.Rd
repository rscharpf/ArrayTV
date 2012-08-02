\name{calctvScore}
\alias{calctvScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calctvScore
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calctvScore(gcFracBoth, samplechr, nparts, useM, narrays, increm, increm2)
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
  \item{increm}{
%%     ~~Describe \code{increm} here~~
}
  \item{increm2}{
%%     ~~Describe \code{increm2} here~~
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
function (gcFracBoth, samplechr, nparts, useM, narrays, increm, 
    increm2) 
{
    narrays = ncol(useM)
    tvScorepart = foreach(i = 1:2, .combine = "rbind") \%dopar\% 
        {
            if (i == 1) {
                partitionedgcFrac = gcFracBoth[, 1]
            }
            else {
                partitionedgcFrac = gcFracBoth[, 2]
            }
            tvScorepart = matrix(0, nrow = nparts, ncol = narrays)
            for (hh in 1:nparts) {
                if (hh == 1) {
                  gcFracA = partitionedgcFrac[seq(hh, length(partitionedgcFrac), 
                    nparts)]
                }
                else {
                  gcFracA = (gcFracA * ((hh - 1)/(hh))) + partitionedgcFrac[seq(hh, 
                    length(partitionedgcFrac), nparts)]/(hh)
                }
                quants = sort(unique(quantile(gcFracA, probs = seq(0, 
                  1, 0.01))))
                gcFracA = quants[findInterval(gcFracA, quants)]
                newsp = split(useM, gcFracA)
                toplot = sapply(newsp, function(x) {
                  x = matrix(x, ncol = narrays)
                  y = colMeans(x)
                  y
                })
                fsampled = colSums(useM)
                n = nrow(useM)
                lambda = fsampled/n
                ngc = sapply(newsp, length)/narrays
                if (is.null(nrow(toplot))) {
                  tvscore = sum(rep((ngc/n), each = narrays) * 
                    (abs(toplot - lambda)))
                }
                else {
                  tvscore = rowSums(rep((ngc/n), each = narrays) * 
                    (abs(toplot - lambda)))
                }
                print(paste("windows stretch", ifelse(i == 1, 
                  hh * increm, hh * increm2), "bp from probes, TV:", 
                  round(tvscore, 4)))
                tvScorepart[hh, ] = na.omit(tvscore)
            }
            tvScorepart
        }
    return(tvScorepart)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
