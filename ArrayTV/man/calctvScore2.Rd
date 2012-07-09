\name{calctvScore2}
\alias{calctvScore2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calctvScore2(useM, gcFracBoth, tvScore, nparts, narrays)
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
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
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
function (useM, gcFracBoth, tvScore, nparts, narrays) 
{
    tvScore2 = foreach(vv = 1:2, .combine = "rbind") \%:\% foreach(ar = 1:narrays, 
        .combine = "cbind") \%dopar\% {
        print(vv)
        print(ar)
        tvScorepart = optimalwin(useM, gcFracBoth, tvScore, vv, 
            ar, nparts)
    }
    return(tvScore2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
