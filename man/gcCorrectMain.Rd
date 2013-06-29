\name{gcCorrectMain}
\alias{gcCorrectMain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gcCorrectMain
}
\description{
Correct gc biases in array data. The user must provide a range of windows and this function will compute the tv score for each window and correct the arrays based upon the window with the highest tv score (i.e with window showing the most gc bias)
}
\usage{
   gcCorrectMain(Ms, chr, starts, samplechr, increms,
                          maxwins, jittercorrection=FALSE,
                          returnOnlyTV=FALSE, onlyGC=FALSE,providedGC=NULL,build, verbose=FALSE)
}
\arguments{
  \item{Ms}{
A matrix holding the array values
}
  \item{chr}{
a vector indicating the chromosome of each probe (should be the same length as the number of rows in the array matrix)
}
  \item{starts}{
a vector holding the genomic coordinates of the start positions the probes (should be the same length as the number of rows inthe array matrix)
}
  \item{samplechr}{
a subset (or the entire set) of the chromosomes represented in the arrays. The probes on these chromosomes will be used to compute the tvScore
}
  \item{increms}{
a vector of integers
}
  \item{maxwins}{
a vector ofintegers with each value >= the corresponding value in increms. Corresponding values of maxwins and increms must have the same quotient. For instance, if increms is c(10,1000) then an acceptable value of maxwins would be(500,5000) because the quotients of the corresponding elements are both 5.
}
\item{jittercorrection}{
if TRUE amplify the correction values by a small amount to decrease the autocorrelation of proximal probes
}
\item{returnOnlyTV}{
if TRUE only return the tvScores rather than applying corrections
}
\item{onlyGC}{
if TRUE return the gc fraction of each window (the result may be subsequently passed into the providedGC argument to be used in correction)
}
\item{providedGC}{
a vector of gc values (fractions from zero to 1) to use for computing corrections. This should be the same length as the number of rows in the array, if provided
}
 \item{build}{
the build matching the array. This is used to load the BSgenome object. If the arrays are annotated based on hg18 this would be 'hg18', if the annotation is based on hg19 this would be 'hg19'.
}
\item{verbose}{
provide output throughout computation
}
}
\details{
see example
}
\value{
The default is to return corrected values for the arrays passed in. This result will have the same dimensions as the first argument (the uncorrected values).if returnOnlyTV is TRUE this will be a matrix of tv score values, one row per window, one column per array. The rownames of the matrix will indicate the size of the window. If onlyGC=TRUE the gc fraction of each window will be returned.
}
\references{
Optimal window finding was inspired by a similar concept applied to gc correction for next generation sequencing experiments:Benjamini Y, Speed TP. Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic Acids Res 2012;40:e72
}
\author{
Eitan Halper-Stromberg
}
\examples{
	path <- system.file("extdata", package="ArrayTV")
	load(file.path(path, "array_logratios.rda"))
	nimblegen[, "M"] <- nimblegen[, "M"]/1000

	increms <- c(10,1000,100e3)
	wins <- c(100,10e3,1e6)
	## nodes may be set to the number of cpus available
	nodes <- 1
	library(doMC)
	registerDoMC(nodes)
	## calculate tv scores in many windows and correct M values
	## using the best window
	nimcM1List <- gcCorrectMain(nimblegen[, "M", drop=FALSE],
				    chr=rep("chr15",nrow(nimblegen)),
				    starts=nimblegen[, "position"],
				    increms=increms,
				    maxwins=wins,
				    build='hg18')
	tvScores <- nimcM1List[['tvScore']]
	winsorted <- as.numeric(rownames(tvScores)[order(tvScores,decreasing=TRUE)])
	logwinsorted <- log(as.numeric(winsorted),10)
	logwinsortdiff <- abs(logwinsorted[1]-logwinsorted)
	## correct M values a 2nd time using a different window size
	nimcM2List <- gcCorrect(nimcM1List[['correctedVals']],
				chr=rep("chr15",nrow(nimblegen)),
				starts=nimblegen[, "position"],
				maxwins=winsorted[logwinsortdiff>=1][1],
				build='hg18')
	## Refer to the vignette for details
}
