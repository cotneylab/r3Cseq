\name{plotInteractionsNearViewpoint}
\alias{plotInteractionsNearViewpoint}
\alias{plotInteractionsNearViewpoint,r3Cseq-method}
\title{Plot candidate interaction regions near the viewpoint}

\description{
  Plot candidate interaction regions near the viewpoint by zooming at 10Mb, 1Mb, 500Kb, and 100Kb respectively 
}

\usage{
plotInteractionsNearViewpoint(object)
}


\arguments{

	\item{object}{
		r3Cseq object. The object is the container of interaction regions produced by \code{\link{getInteractions}} function.
	}
}


\value{
	Plots interaction regions close to the viewpoint on the active graphical device.
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotOverviewInteractions}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plot3Cecdf}}
}


\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	getInteractions(my.data)

####Plot
	
	plotInteractionsNearViewpoint(my.data)
	
}
%\keyword{plot interaction regions}
