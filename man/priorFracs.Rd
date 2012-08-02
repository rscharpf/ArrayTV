\name{priorFracs}
\alias{priorFracs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{priorFracs
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
priorFracs(gcFracBoth, maxuse, nparts, tvScore, increm, increm2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gcFracBoth}{
%%     ~~Describe \code{gcFracBoth} here~~
}
  \item{maxuse}{
%%     ~~Describe \code{maxuse} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
}
  \item{tvScore}{
%%     ~~Describe \code{tvScore} here~~
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
function (gcFracBoth, maxuse, nparts, tvScore, increm, increm2) 
{
    offset = as.numeric(rownames(tvScore)[maxuse])/ifelse(maxuse > 
        nrow(tvScore)/2, increm2, increm)
    if (maxuse > nrow(tvScore)/2) 
        priorGC = gcFracBoth[, 2]
    else {
        priorGC = gcFracBoth[, 1]
    }
    priorFrac = (cumsum(priorGC)[offset:length(priorGC)] - c(0, 
        cumsum(priorGC)[1:(length(priorGC) - offset)]))[seq(1, 
        (length(priorGC) - offset + 1), nparts)]/offset
    quants = sort(unique(quantile(priorFrac, probs = seq(0, 1, 
        0.01))))
    priorFrac = quants[findInterval(priorFrac, quants)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
