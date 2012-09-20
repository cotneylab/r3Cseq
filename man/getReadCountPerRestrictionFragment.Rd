\name{getReadCountPerRestrictionFragment}
\alias{getReadCountPerRestrictionFragment}
\alias{getReadCountPerRestrictionFragment,r3Cseq-method}
\title{get read count per restriction fragment for 3C-seq data}

\description{
  Counts the number of reads from 3C-seq data per each restriction fragment
}

\usage{
getReadCountPerRestrictionFragment(object)
}


\arguments{

	\item{object}{
		r3Cseq object, The initialize input parameters are required
	}
}

\value{
   The RangedData represents the number of reads per each restriction fragment
}

\seealso{
  \link{getCoverage},
}

\author{
	S. Thongjuea
}

\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    expFile<-system.file("extdata","alignedReads.fetal.liver.subset.bam",package="r3Cseq")
    contrFile<-system.file("extdata","alignedReads.fetal.brain.subset.bam",package="r3Cseq")
	
	my.data<-new("r3Cseq",organismName='mm9',alignedReadsBamExpFile=expFile,
	alignedReadsBamContrFile=contrFile,isControlInvolved=TRUE,isBamInputFile=TRUE,
	expLabel="fetal_liver",contrLabel="fetal_brain",
	restrictionEnzyme='HindIII')
	
	getReadCountPerRestrictionFragment(my.data)
	
}
%\keyword{read counts manipulation}
