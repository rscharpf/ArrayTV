\name{correctionTVscore}
\alias{correctionTVscore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correctionTVscore
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correctionTVscore(useM, gcFrac, arrayn, windowl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{useM}{
%%     ~~Describe \code{useM} here~~
}
  \item{gcFrac}{
%%     ~~Describe \code{gcFrac} here~~
}
  \item{arrayn}{
%%     ~~Describe \code{arrayn} here~~
}
  \item{windowl}{
%%     ~~Describe \code{windowl} here~~
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
function (useM, gcFrac, arrayn, windowl) 
{
    newsp = split(useM, gcFrac)
    correctionVals = sapply(newsp, mean)
    fsampled = sum(useM)
    n = length(useM)
    lambda = fsampled/n
    ngc = sapply(newsp, length)
    tvscore = sum(ngc/n * (abs(correctionVals - lambda)))
    return(tvscore)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
