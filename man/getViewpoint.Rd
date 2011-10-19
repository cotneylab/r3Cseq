\name{getViewpoint}
\alias{getViewpoint}
\alias{getViewpoint,r3Cseq-method}
\title{get the viewpoint of 3C-seq data}

\description{
  The viewpoint is the region of interest, which can be a promoter region of an interested gene, 
  an enhancer, and a transcrition factor binding region. The genomic region, which show the highest number of reads per
  restriction fragment, represents the viewpoint.
   
}

\usage{
getViewpoint(object)
}


\arguments{

	\item{object}{
		r3Cseq object, the object is the container of interaction regions produced by \code{\link{getInteractions}} function.
	}
}

\value{
   The viewpoint shows in the IRanges
}

\seealso{
  \link{getCoverage},
  \link{getReadCountPerRestrictionFragment},
  \link{calculateRPM}
}

\author{
	S. Thongjuea
}

\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	viewpoint<-getViewpoint(my.data)
	viewpoint
	
}
%\keyword{viewpoint}
