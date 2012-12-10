\name{plotDomainogramNearViewpoint}
\alias{plotDomainogramNearViewpoint}
\alias{plotDomainogramNearViewpoint,r3Cseq-method}

\title{Plot domainogram of interaction regions near the viewpoint}

\description{
  Plot domainogram of interaction regions near the viewpoint 
}

\usage{
plotDomainogramNearViewpoint(object,distance=5e5,maximum_window=25e3,view=c("experiment","control","both"))
}


\arguments{

	\item{object}{
		r3Cseq or r3CseqInBatch object
	}
	\item{distance}{
		Numeric. The distance relative to the viewpoint (default=5e5)
	}
	\item{maximum_window}{
		Numeric. The maximum windowing (default=25e3). We normally compute the interaction regions per window starting from 2Kb to maximum window (default=25kb) to make the interaction matrix
		for visualizing the domainogram.
	}
	\item{view}{
		character. The selected view of data (default="experiment")
	}
}


\value{
	Plots of domainogram for interaction regions close to the viewpoint
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotOverviewInteractions}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plotInteractionsNearViewpoint}}
}


\examples{

	# See the vignette
}
%\keyword{plot interaction regions}
