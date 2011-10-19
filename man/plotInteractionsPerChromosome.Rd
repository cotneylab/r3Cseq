\name{plotInteractionsPerChromosome}
\alias{plotInteractionsPerChromosome}
\alias{plotInteractionsPerChromosome,r3Cseq-method}
\title{Plot interaction regions per each chromosome of interest}

\description{
  Plot the distribution of interaction regions per each interested input chromosome
}

\usage{
plotInteractionsPerChromosome(object, chromosomeName)
}


\arguments{

	\item{object}{
		r3Cseq object. The object is the container of interaction regions produced by \code{\link{getInteractions}} function.
	}

	\item{chromosomeName}{
		Character. The input chromosome name (e.g. "chr1")
	}
}


\value{
	Plots interaction regions per chromosome on the active graphical device.
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotInteractionsNearViewpoint}}, \code{\link{plotOverviewInteractions}}, 
	\code{\link{plot3Cecdf}}
}


\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	getInteractions(my.data)

####Plot
	
	plotInteractionsPerChromosome(my.data,"chr10")
	
}
%\keyword{plot interaction regions}
