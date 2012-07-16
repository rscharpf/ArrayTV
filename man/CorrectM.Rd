\name{CorrectM}
\alias{CorrectM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correct M values based on median of optimal gc win
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrectM(gcFracBoth, useM, Ms, starts, narrays, nparts, chr, samplechr, remainingChr, increm, formaxs1, formaxs2, revmaxs1, revmaxs2, maxorigins, gmaxvals, gmaxvals2)
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
  \item{formaxs1}{
%%     ~~Describe \code{formaxs1} here~~
}
  \item{formaxs2}{
%%     ~~Describe \code{formaxs2} here~~
}
  \item{revmaxs1}{
%%     ~~Describe \code{revmaxs1} here~~
}
  \item{revmaxs2}{
%%     ~~Describe \code{revmaxs2} here~~
}
  \item{maxorigins}{
%%     ~~Describe \code{maxorigins} here~~
}
  \item{gmaxvals}{
%%     ~~Describe \code{gmaxvals} here~~
}
  \item{gmaxvals2}{
%%     ~~Describe \code{gmaxvals2} here~~
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
    samplechr, remainingChr, increm, formaxs1, formaxs2, revmaxs1, 
    revmaxs2, maxorigins, gmaxvals, gmaxvals2) 
{
    correctedM = foreach(i = 1:narrays, .combine = "cbind") \%dopar\% 
        {
            print(paste("i is", i, "and high tv is", ifelse(maxorigins[i] == 
                F, gmaxvals2[i], gmaxvals[i])))
            ar = i
            reversedirMax2 = ifelse(maxorigins[i] == F, ifelse(revmaxs2[i, 
                2] > formaxs2[i, 2], T, F), F)
            reversedirMax = revmaxs1[i, 2] > formaxs1[i, 2] & 
                reversedirMax2 == F
            maxuse1 = ifelse(reversedirMax, revmaxs1[i, 1], formaxs1[i, 
                1])
            vv1 = ifelse(reversedirMax, 2, 1)
            v1 = c(3, 1)
            v2 = c(4, 2)
            priorFrac = priorFracs(gcFracBoth, maxuse1, vv1, 
                nparts, v1, v2)
            if (maxorigins[i] == F) {
                maxuse2 = ifelse(reversedirMax2, revmaxs2[i, 
                  1], formaxs2[i, 1])
                vv2 = ifelse(reversedirMax2, 2, 1)
                priorFrac2 = priorFracs(gcFracBoth, maxuse2, 
                  vv2, nparts, v1, v2)
                priorFrac = priorFrac * (maxuse1/(maxuse1 + maxuse2)) + 
                  priorFrac2 * (maxuse2/(maxuse1 + maxuse2))
            }
            priorFrac = round(priorFrac, 2)
            newsp = split(useM[, i], priorFrac)
            correctionVals = sapply(newsp, median)
            names(correctionVals) = names(newsp)
            fsampled = sum(useM[, i])
            n = nrow(useM)
            lambda = fsampled/n
            ngc = sapply(newsp, length)
            tvscore = sum(ngc/n * (abs(correctionVals - lambda)))
            print(paste("i is now", i, "tv score is", tvscore))
            forwardExtend = ifelse(!reversedirMax, maxuse1 * 
                increm, ifelse(!reversedirMax2 && maxorigins[i] == 
                F, maxuse2 * increm, 1))
            reverseExtend = ifelse(reversedirMax, maxuse1 * increm, 
                ifelse(reversedirMax2, maxuse2 * increm, 0))
            schr = chr[chr \%in\% samplechr]
            priorFracWremaining = round(priorFracsRestOfGenome(forwardExtend, 
                reverseExtend, remainingChr, starts, chr), 2)
            priorFracWremaining[chr \%in\% samplechr] = priorFrac
            newsp = split(Ms[, i], paste(chr, priorFracWremaining, 
                sep = "."))
            correctionVals = sapply(newsp, median)
            names(correctionVals) = names(newsp)
            fsampled = sum(Ms[, i])
            n = nrow(Ms)
            lambda = fsampled/n
            ngc = sapply(newsp, length)
            tvscore = sum(ngc/n * (abs(correctionVals - lambda)))
            print(paste("i is now", i, "NEW* tv score is", tvscore))
            correctedM = Ms[, i] - correctionVals[match(paste(chr, 
                priorFracWremaining, sep = "."), (names(correctionVals)))]
            if (length(which(is.na(correctedM))) > 0) {
                print(priorFracWremaining[is.na(correctedM)][1:20])
            }
            print(paste("mad is", mad(correctedM)))
            print(paste("old mad is", mad(Ms[, i])))
            correctedM
        }
    return(correctedM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
