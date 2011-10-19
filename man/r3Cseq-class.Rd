\name{r3Cseq-class}
\docType{class}

% Classes:
\alias{r3Cseq}
\alias{r3Cseq-class}

\title{r3Cseq objects}

\description{
  The r3Cseq class is a general container for storing
  and manipulating a set of input parameters, RangeData of interactions regions from r3Cseq analysis
  , and Rle coverage vector of the genome-wide interaction signal. 
}

\section{Slots}{
  \describe{

    \item{\code{organismName}}{Object of class \code{"character"} the version of 
      particular assembly genome from UCSC (e.g. mm9, hg18, hg19) . The package supports
      three genome assemblies consisting of mouse (mm9), and human (hg18, hg19).}

    \item{\code{restrictionEnzyme}}{Object of class \code{"character"} the
      restriction enzyme name using in cutting the 3C interaction fragments}

    \item{\code{alignedReadsExpFile}}{Object of class \code{"character"} the
    file name of the experiment}

    \item{\code{alignedReadsContrFile}}{Object of class \code{"character"}
      the file name of the control}
    
    \item{\code{alignedReadsBamExpFile}}{Object of class \code{"character"}
      the file name of experiment in BAM format}
      
    \item{\code{alignedReadsBamContrFile}}{Object of class \code{"character"}
      the file name of control in BAM format}
      
    \item{\code{alignedReadsType}}{Object of class \code{"character"}
      the file type of aligned read (see \link{readAligned})}
    
    \item{\code{expLabel}}{Object of class \code{"character"}
      the experiment name}
      
    \item{\code{contrLabel}}{Object of class \code{"character"}
      the control name}
    
    \item{\code{expLibrarySize }}{Object of class \code{"integer"}
      the library size of experiment}  
      
    \item{\code{contrLibrarySize }}{Object of class \code{"integer"}
      the library size of control}  
    
    \item{\code{expReadLength}}{Object of class \code{"integer"}
      the read length of experiment}  
      
    \item{\code{contrReadLength}}{Object of class \code{"integer"}
      the read length of experiment} 
    
    \item{\code{expReadCount}}{Object of class \code{"RangedData"}
      the read count in experiment}  
    
    \item{\code{contrReadCount}}{Object of class \code{"RangedData"}
      the read count in control}  
    
    \item{\code{expRPM}}{Object of class \code{"RangedData"}
      the read read per million in experiment}  
   
    \item{\code{contrRPM}}{Object of class \code{"RangedData"}
      the read read per million in control}  
   
    \item{\code{expInteractionRegions}}{Object of class \code{"RangedData"}
      the interaction regions in experiment}  
   
    \item{\code{contrInteractionRegions}}{Object of class \code{"RangedData"}
      the interaction regions in control}  
      
    \item{\code{expCoverage}}{Object of class \code{"RleList"}
      the Rle coverage vector in experiment}  				
   
    \item{\code{contrCoverage}}{Object of class \code{"RleList"}
      the Rle coverage vector in control}  				
 	
 	\item{\code{isControlInvolved}}{Object of class \code{"logical"}
      the logical to ask whether the control is involved in the analysis or not}  		
	
	\item{\code{isBamInputFile}}{Object of class \code{"logical"}
      the logical to ask whether the BAM file format is the input file or not}  			
  }
}


\section{Methods}{
   
  \describe{
    \item{}{
      \code{expInteractionRegions(x)}:
      The interaction regions in experiment \code{x}.
    }
    \item{}{
      \code{contrInteractionRegions(x)}:
     The interaction regions in control \code{x}.
    }
  }
}

\author{S. Thongjuea}

\seealso{
  \link{getCoverage},
  \link{getReadCountPerRestrictionFragment},
  \link{calculateRPM},
  \link{getInteractions},
  \link{expInteractionRegions},
  \link{contrInteractionRegions}
}

\examples{
	library(BSgenome.Mmusculus.UCSC.mm9)
    expFile<-system.file("extdata","alignedReads.fetal.liver.subset.bam",package="r3Cseq")
    contrFile<-system.file("extdata","alignedReads.fetal.brain.subset.bam",package="r3Cseq")
	
	my.data<-new("r3Cseq",organismName='mm9',alignedReadsBamExpFile=expFile,
	alignedReadsBamContrFile=contrFile,isControlInvolved=TRUE,isBamInputFile=TRUE,
	expLabel="fetal_liver",contrLabel="fetal_brain",
	restrictionEnzyme='HindIII')
	my.data
}
\keyword{classes}
