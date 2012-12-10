\name{getBatchInteractions}
\alias{getBatchInteractions}
\alias{getBatchInteractions,r3CseqInBatch-method}

\title{calculate z-score, assign p-value and q-value for each interaction region for replicates data sets}

\description{
  Calculate z-score, assign p-value and q-value to each interaction regions for replicates data sets
}

\usage{
getBatchInteractions(object,method=c("union","intersection"),fdr=0.05)
}


\arguments{

	\item{object}{
		r3Cseq object
	}
	\item{method}{
		character. The method for combining biological replicates for 3C-Seq analysis (default = "union")
	}
	\item{fdr}{
		A level at which to control the FDR. Must be in (0,1] (default=0.05)
	}
}

\value{
   The interaction regions show in the RangedData
}

\seealso{
  \link{getInteractions}
}

\author{
	S. Thongjuea
}

\examples{
 #See the vignette	
}
%\keyword{3C-Seq interaction region analysis}
