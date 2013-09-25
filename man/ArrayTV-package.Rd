\name{ArrayTV-package}
\alias{ArrayTV-package}
\alias{ArrayTV}
\docType{package}
\title{
ArrayTV-implementation of GC correction for arrays
}
\description{
Correct gc biases in array data. The user must provide a range of windows and this function will compute the tv score for each window and correct the arrays based upon the window with the highest tv score (i.e with window showing the most gc bias)
}
\details{
\tabular{ll}{
Package: \tab ArrayTV\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-08-02\cr
License: \tab Johsn Hopkins School of Medicine\cr
}
gcCorrectMain takes a matrix of array signal intensities, each column
representing one array, along with vectors corresponding to probe
positions and chromosomes. The gcCorrect methods enables array input for
objects of class matrix, BeadStudioSet,BafLrrSet, and BafLrrSetList
}
\author{
Eitan Halper-Stromberg
Maintainer: Eitan Halper-Stromberg
}
\references{
  Yuval Benjamini and Terence P. Speed. Summarizing and correcting the
  GC content bias in high- throughput sequencing. Nucleic Acids Res,
  40(10):e72, May 2012.

  Eitan Halper-Stromberg, Laurence Frelin, Ingo Ruczinski, Robert
  Scharpf, Chunfa Jie, Benilton Carvalho, Haiping Hao, Kurt Hetrick,
  Anne Jedlicka, Amanda Dziedzic, Kim Doheny, Alan F. Scott, Steve
  Baylin, Jonathan Pevsner, Forrest Spencer, and Rafael
  A. Irizarry. Performance assessment of copy number microarray
  platforms using a spike-in experiment. Bioinformatics, 27(8):1052

  Adam B Olshen, E. S. Venkatraman, Robert Lucito, and Michael
  Wigler. Circular binary segmentation for the analysis of array-based
  DNA copy number data. Biostatistics, 5(4):557
}
\keyword{ package }
\seealso{
type ?gcCorrectMain.Rd at the command prompt for an example
}
