\name{plotOverviewInteractions}
\alias{plotOverviewInteractions}
\alias{plotOverviewInteractions,r3Cseq-method}
\title{Plot the overview of interaction regions for the genome-wide}

\description{
  Plot the distribution of interaction regions across genome with the level of interaction signal
}

\usage{
plotOverviewInteractions(object, cutoff.p_value=0.05,cutoff.fold_change=2)
}


\arguments{

	\item{object}{
		r3Cseq object. The object is the container of interaction regions produced by \code{\link{getInteractions}} function.
	}

	\item{cutoff.p_value}{
		Numeric. The cutoff p-value from empirical distribution function (default=0.05)
	}

	\item{cutoff.fold_change}{
		Numeric. The cutoff fold change compare between experiment and control (default=2)
	}
}


\value{
	Plots interaction regions genome-wide on the active graphical device.
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotInteractionsNearViewpoint}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plot3Cecdf}}
}


\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	getInteractions(my.data)

####Plot
	
	plotOverviewInteractions(my.data)
	
}
%\keyword{plot interaction regions}
