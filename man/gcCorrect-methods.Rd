\name{gcCorrect-methods}
\docType{methods}
\alias{gcCorrect}
\alias{gcCorrect-methods}
\alias{gcCorrect,BafLrrSetList-method}
\alias{gcCorrect,BafLrrSet-method}
\alias{gcCorrect,BeadStudioSet-method}
\alias{gcCorrect,matrix-method}
\title{Methods for \code{gcCorrect} in Package \pkg{ArrayTV} }
\description{
Correct gc biases in array data for arrays stored as objects of an acceptable class. The user provides the array data and a range of windows and these methods will format the data and call gcCorrectMain. This function computes the tv score for each window and corrects the arrays based upon the window with the highest tv score (i.e with window showing the most gc bias)


}
\section{Methods}{
  \describe{

    \item{}{\code{gcCorrect(object, ...)}:

      Methods have been defined for objects of class \code{matrix},
      \code{BeadStudioSet}, \code{BafLrrSet}, and \code{BafLrrSetList}.

      Objects of class \code{BafLrrSetList} are containers for log R
      ratios and B allele frequencies (BAFs) stored by chromosome.  In
      particular, each element in this list class contains BAFs and log
      R ratios as \code{assayData} and genomic annotation of the markers
      (\code{featureData}) for one chromosome.

      The \code{gcCorrect} method for objects of this class extracts
      the log R ratios for all elements (chromosomes) in the list, and
      then combines these into a single matrix.  If the log R ratios
      are stored as \code{ff} objects, the samples will be
      GC-corrected in chunks determined by the function
      \code{ocSamples()}.  For example, if \code{object} contains 100
      samples and \code{ocSamples(10)} was specified prior to calling
      \code{gcCorrect}, the wave correction will be performed on 10
      samples at a time in order to keep RAM at an acceptable level.
      When processing large numbers of samples, it can be helpful to
      evaluate the tv score on a subset of the samples, pick one or
      two windows for the correction, and precompute the GC content
      for these windows. See the examples below for the implementation
      with \code{BafLrrSetList} objects for details.  When the assay
      data is represented as \code{ff} objects, the \code{gcCorrect}
      method returns the value \code{NULL}.  Note that the assay data
      stored on disk will have changed as a result of calling this
      method.

      When \code{object} is a \code{matrix} of log R ratios, the
      \code{gcCorrect} method returns a matrix of wave-adjusted log R
      ratios.

      Additional arguments are passed to \code{gcCorrectMain} through
      the \code{...} operator.

    }
  }
}

\seealso{\code{\linkS4class{BeadStudioSet}},
  \code{\linkS4class{BafLrrSetList}}, \code{\linkS4class{BafLrrSet}} \code{\link{gcCorrectMain}}}

\references{
  Benjamini Y, Speed TP. Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic
  Acids Res 2012;40:e72
}

\keyword{methods}
\keyword{regression}
\keyword{smooth}

\examples{
### Example with 2 iterations of correction
\dontrun{
path <- system.file("extdata", package="ArrayTV")
load(file.path(path, "array_logratios.rda"))
nimblegen[, "M"] <- nimblegen[, "M"]/1000
increms <- c(10,1000,100e3)
wins <- c(100,10e3,1e6)
## nodes may be set to the number of cpus available
nodes=1
library(doMC)
registerDoMC(nodes)
## calculate tv scores in many windows and correct M values using the best window
nimcM1List<-gcCorrectMain(nimblegen[, "M",
drop=FALSE],chr=rep("chr15",nrow(nimblegen)),starts=nimblegen[,
"position"],increms=increms,maxwins=wins,build='hg18')
tvScores<-nimcM1List[['tvScore']]
winsorted<-as.numeric(rownames(tvScores)[order(tvScores,decreasing=TRUE)])
logwinsorted<-log(as.numeric(winsorted),10)
logwinsortdiff<-abs(logwinsorted[1]-logwinsorted)
## correct M values a 2nd time using a different window size


nimcM2List<-gcCorrect(nimcM1List[['correctedVals']],chr=rep("chr15",nrow(nimblegen)),starts=nimblegen[,
"position"],maxwins=winsorted[logwinsortdiff>=1][1],build='hg18')

## Refer to the vignette for details

## Example using a list container (BafLrrSetList) containing log R
## ratios and B allele frequencies

	if(require(crlmm) && require(ff)){
		data(cnSetExample, package="crlmm")
		brList <- constructBafLrrSetListFrom(cnSetExample)
		is(lrr(brList)[[1]], "ff_matrix")
		rold <- lrr(brList)[[1]][, 1]/100
		##
		## To avoid having too much data in RAM it might be
		## useful to process the samples n at a time
		##
		##  The number of samples to be processed at a time is
		##  set by the ocSamples function.  For example, to
		##  process 20 samples at a time one would do
		ocSamples(20)
		##
		## When assay data elements are ff objects, the wave
		## corrected values are updated on disk.  Currently,
		## the  data is stored is here:
		filename(lrr(brList)[[1]])
		##
		## To avoid permanently changing the log R ratio
		## values for the brList object, we copy the ff files
		## to a different path and create a new BafLrrSetList
		## object.  This is done by the non-exported function
		## "duplicateBLList". The path for the new ff objects
		## is given by the function ldPath().  Here, we use a
		## temp directory
		ldPath(tempdir())
		brList.copy <- oligoClasses:::duplicateBLList(brList, ids=sampleNames(brList))
		filename(lrr(brList.copy)[[1]])
		##
		## wave correct the log R ratios
		##
		gcCorrect(brList.copy, increms=c(12, 10e3),
			  maxwins=c(12, 10e3))
		##
		##
		rupdated <- lrr(brList.copy)[[1]][, 1]/100
		## note that rold and rupdated are no longer the same
		## since the log R ratios were updated in the
		## brList.copy container
		plot(rupdated, rold, pch=".")
		##
		## Remarks on efficiency
		##
		## If a large number of samples are to be processed,
		## the most efficient procedure is to settle on an
		## appropriate window size for gc correction using a
		## subset of the dataset (e.g., 20 samples).  See the
		## vignette for details.  Here, we assume that two
		## windows were already selected using such a
		## procedure (here, 12 bp and 10,000 bp) and the goal
		## is to efficiently do wave correction using these
		## two windows for all the samples in a large study.
		## Currently, our recommended approach is to first
		## calculate the gc content for these windows:
		gc.matrix <- computeGC(brList.copy, c(12, 10e3),
				       c(12, 10e3))
		## The number of columns in gc.matrix will correspond
		## to the number of windows
		ncol(gc.matrix)
		## Having calculated the gc content for these two
		## windows, we pass the gc matrix to the method
		## gcCorrect.  This function will iteratively update
		## the log R ratios for the gc content given by the
		## columns in gc.matrix (the log R ratios are updated
		## for each column in gc.matrix).
		ldPath(tempdir())
		brList.copy <- oligoClasses:::duplicateBLList(brList, ids=sampleNames(brList))
		gcCorrect(brList.copy, increms=c(12, 10e3),
			  maxwins=c(12, 10e3),
			  providedGC=gc.matrix)
	}
}
}
