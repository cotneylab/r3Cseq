setClass(
		Class="r3Cseq",
		representation(
				organismName="character",
				restrictionEnzyme="character",
				alignedReadsExpFile="character",
				alignedReadsContrFile="character",
				alignedReadsBamExpFile="character",
				alignedReadsBamContrFile="character",
				alignedReadsType ="character",
				expLabel ="character",
				contrLabel="character",
				expLibrarySize ="integer",
				contrLibrarySize ="integer",
				expReadLength="integer",
				contrReadLength="integer",
				expReadCount="RangedData",
				contrReadCount="RangedData",
				expRPM="RangedData",
				contrRPM="RangedData",
				expInteractionRegions="RangedData",
				contrInteractionRegions="RangedData",
				expCoverage="RleList",
				contrCoverage="RleList",
				isControlInvolved="logical",
				isBamInputFile="logical"
		),
		prototype(
				organismName=character(1),
				restrictionEnzyme=character(1),
				alignedReadsExpFile=character(1),
				alignedReadsContrFile=character(1),
				alignedReadsBamExpFile=character(1),
				alignedReadsBamContrFile=character(1),
				alignedReadsType=character(1),
				expLabel =character(1),
				contrLabel=character(1),
				expReadCount=RangedData(),
				contrReadCount=RangedData(),
				expRPM=RangedData(),
				contrRPM=RangedData(),
				expCoverage=RleList(),
				contrCoverage=RleList(),
				expInteractionRegions=RangedData(),
				contrInteractionRegions=RangedData(),
				isControlInvolved=FALSE,
				isBamInputFile=FALSE
		),
		validity = function(object) {
			if(object@organismName==character(1))
				return( "The organism name is empty" )
			if(!object@organismName %in% c("mm9","hg18","hg19"))
				return( "This version of r3Cseq supports only mm9, hg18 and hg19 assembly." )
			if(object@restrictionEnzyme==character(1))
				return( "The restrictionEnzyme is empty" )
		}
)

setClass("repbaseEnzyme",
		representation(
				enzymeRestriction="data.frame"
		),
		prototype(
				enzymeRestriction=data.frame()
		)
)
