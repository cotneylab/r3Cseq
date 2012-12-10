\name{getInteractions}
\alias{getInteractions}
\alias{getInteractions,r3Cseq-method}
\title{calculate z-score, assign p-value and q-value for each interaction region}

\description{
  Calculate z-score, assign p-value and q-value to each interaction regions
}

\usage{
getInteractions(object,fdr)
}


\arguments{

	\item{object}{
		r3Cseq object
	}
	\item{fdr}{
		A level at which to control the FDR. Must be in (0,1] (default=0.05)
	}
}

\value{
   The interaction regions show in the RangedData
}

\seealso{
  \link{getBatchInteractions}
}

\author{
	S. Thongjuea
}

\examples{
 #See the vignette	
}
%\keyword{3C-Seq interaction region analysis}
