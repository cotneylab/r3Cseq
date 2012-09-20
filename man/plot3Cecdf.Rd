\name{plot3Cecdf}
\alias{plot3Cecdf}
\alias{plot3Cecdf,r3Cseq-method}
\title{Plot the empirical distribution of interaction regions}

\description{
  Plot the empirical distribution of interaction regions
}

\usage{
plot3Cecdf(object)
}


\arguments{

	\item{object}{
		r3Cseq object. The object is the container of interaction regions produced by \code{\link{getInteractions}} function.
	}
}


\value{
	Plots the empirical distribution of interaction regions on the active graphical device.
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotOverviewInteractions}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plotInteractionsNearViewpoint}}
}


\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	getInteractions(my.data)

####Plot
	
	plot3Cecdf(my.data)
	
}
%\keyword{plot interaction regions}
