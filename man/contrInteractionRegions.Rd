\name{contrInteractionRegions}
\alias{contrInteractionRegions}
\alias{contrInteractionRegions,r3Cseq-method}
\title{get interaction regions from the control}

\description{
  get all candidate interaction regions from the control
}

\usage{
contrInteractionRegions(object)
}


\arguments{

	\item{object}{
		r3Cseq object. The object might contain the interaction regions generated by function \link{getInteractions}
	}
}

\value{
   The candidate interaction regions show in the IRange object
}

\seealso{
  \link{getCoverage},
  \link{getReadCountPerRestrictionFragment},
  \link{calculateRPM},
  \link{getInteractions}
}

\author{
	S. Thongjuea
}

\examples{
 	
####Create the r3Cseq object#############
	library(BSgenome.Mmusculus.UCSC.mm9)
    load(system.file("data","example.data.rda",package="r3Cseq"))
	calculateRPM(my.data)
	getInteractions(my.data)
	contr.interactions<-contrInteractionRegions(my.data)
}
%\keyword{}