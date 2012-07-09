\name{priorFracs}
\alias{priorFracs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
priorFracs(gcFracBoth, maxuse, vv, nparts, v1, v2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gcFracBoth}{
%%     ~~Describe \code{gcFracBoth} here~~
}
  \item{maxuse}{
%%     ~~Describe \code{maxuse} here~~
}
  \item{vv}{
%%     ~~Describe \code{vv} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
}
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
%%     ~~Describe \code{v2} here~~
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
function (gcFracBoth, maxuse, vv, nparts, v1, v2) 
{
    optnum = 4
    extras = maxuse/(nparts/optnum)
    gcAdd = 0
    if (extras > 1) {
        mo = maxuse\%\%(nparts/optnum)
        if (mo == 0) {
            gcAdd = gcFracBoth[[v2[-vv]]]
        }
        else {
            gcAdd = gcFracBoth[[v2[-vv]]]
            mods = (1:length(gcAdd))\%\%(nparts/optnum)
            gcAdd[mods > mo | mods == 0] = 0
        }
    }
    gcFracUse = gcFracBoth[[v1[-vv]]] + gcAdd
    priorGC = gcFracUse
    offset = ifelse(maxuse > (nparts/optnum), nparts/optnum, 
        maxuse)
    priorFrac = (cumsum(priorGC)[offset:length(priorGC)] - c(0, 
        cumsum(priorGC)[1:(length(priorGC) - offset)]))[seq(1, 
        (length(priorGC) - offset + 1), nparts/optnum)]/maxuse
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
