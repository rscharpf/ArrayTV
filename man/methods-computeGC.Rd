\name{computeGC-methods}
\docType{methods}
\alias{computeGC-methods}
\alias{computeGC,BafLrrSetList-method}
\alias{computeGC,numeric-method}
\title{Compute the GC composition for selected windows}
\description{

  Compute the GC composition for selected windows around a vector of
  genomic positions (e.g., locations of markers on an array)

}
\section{Methods}{
  \describe{

    \item{}{\code{computeGC(x, maxwins, increms, chr, build, ...)}:

      \code{x} can be numeric, in which case \code{x} should be a vector
      of physical position along the genome.  Alternatively, \code{x}
      can be an instance of the \code{BafLrrSetList} class.  When
      \code{x} is \code{numeric}, arguments \code{chr} and \code{build}
      are required.  Arguments \code{chr} and \code{build} are ignored
      for \code{BafLrrSetList} objects.

      Arguments \code{maxwins}, \code{increms}, \code{build}, and
      \code{chr} are as described in the documentation for
      \code{gcCorrectMain}
    }
  }
}
\seealso{\code{\linkS4class{BafLrrSetList}},
  \code{\link{gcCorrectMain}}}

\keyword{methods}

\examples{
\dontrun{
    if(require(crlmm) && require(ff)){
	data(cnSetExample, package="crlmm")
	brList <- BafLrrSetList(cnSetExample)
	pos <- unlist(position(brList))
	chr <- paste("chr",
		     rep(chromosome(brList), elementLengths(brList)),
		     sep="")
	gc.matrix1 <- ArrayTV:::computeGC(x=pos, maxwins=c(12,10e3),
           				  increms=c(12,10e3),
		         		  chr=chr,
				          build="hg18")
	gc.matrix2 <- ArrayTV:::computeGC(brList, c(12, 10e3),
         			          c(12, 10e3), build="hg18")
        identical(gc.matrix1, gc.matrix2)
	}
}
}
