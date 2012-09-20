\name{contrCoverage}
\alias{contrCoverage}
\alias{contrCoverage,r3Cseq-method}
\title{Get the coverage for 3C-seq data from the control}

\description{
  Counts the number of times a position is represented in a set of input reads from 3C-seq data in the control
}

\usage{
contrCoverage(object)
}


\arguments{

	\item{object}{
		r3Cseq object, The initialize input parameters are required
	}
}

\value{
  For most methods, an \link{Rle} object representing the coverage of
  \code{x}. For \link{RangesList} and \link{RangedData} objects, a
  \link{SimpleRleList} object representing a list of coverage vectors.

  An integer value called the "coverage" can be associated to each position
  in \code{x}, indicating how many times this position is covered by the
  elements contained in \code{x}.
  For example, if \code{x} is a \link{Views} object, the coverage of
  a given position in \code{subject(x)} is the number of views it belongs to.
}

\seealso{
  \link{IRanges-class},
  \link{Views-class},
  \link{Rle-class},
  \link{MaskCollection-class}
}

\author{
	S. Thongjuea
}

\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
	load(system.file("data","example.data.rda",package="r3Cseq"))
	contr.coverage<-contrCoverage(my.data)
}
%\keyword{coverage manipulation}
