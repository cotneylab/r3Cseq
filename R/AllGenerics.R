# Author: Supat Thongjuea
####################################
#########  AllGenerics.R
#########
#########  all generics in r3Cseq
####################################
setGeneric(
		name="organismName",
		def=function(object){
			standardGeneric("organismName")
		})
setMethod("organismName",
		signature(object = "r3Cseq"),
		function (object){
			object@organismName
		}
)
########
setGeneric(
		name="restrictionEnzyme",
		def=function(object){
			standardGeneric("restrictionEnzyme")
		})
setMethod("restrictionEnzyme",
		signature(object = "r3Cseq"),
		function (object){
			object@restrictionEnzyme
		}
)
#######
setGeneric(
		name="alignedReadsExpFile",
		def=function(object){
			standardGeneric("alignedReadsExpFile")
		})
setMethod("alignedReadsExpFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsExpFile
		}
)
#######
setGeneric(
		name="alignedReadsContrFile",
		def=function(object){
			standardGeneric("alignedReadsContrFile")
		})

setMethod("alignedReadsContrFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsContrFile
		}
)
#######
setGeneric(
		name="alignedReadsBamExpFile",
		def=function(object){
			standardGeneric("alignedReadsBamExpFile")
		})
setMethod("alignedReadsBamExpFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsBamExpFile
		}
)
#######
setGeneric(
		name="alignedReadsBamContrFile",
		def=function(object){
			standardGeneric("alignedReadsBamContrFile")
		})
setMethod("alignedReadsBamContrFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsBamContrFile
		}
)

########
setGeneric(
		name="alignedReadsType",
		def=function(object){
			standardGeneric("alignedReadsType")
		})
setMethod("alignedReadsType",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsType
		}
)
########
setGeneric(
		name="expLabel",
		def=function(object){
			standardGeneric("expLabel")
		})
setMethod("expLabel",
		signature(object = "r3Cseq"),
		function (object){
			object@expLabel
		}
)
setGeneric(
		name="contrLabel",
		def=function(object){
			standardGeneric("contrLabel")
		})
setMethod("contrLabel",
		signature(object = "r3Cseq"),
		function (object){
			object@contrLabel
		}
)
########
setGeneric(
		name="expLibrarySize",
		def=function(object){
			standardGeneric("expLibrarySize")
		}
)
setGeneric(
		name="expLibrarySize<-",
		def=function(object,value){
			standardGeneric("expLibrarySize<-")
		}
)
setMethod("expLibrarySize",
		signature(object="r3Cseq"),
		function(object){
			object@expLibrarySize
		}
)
setReplaceMethod(
		f="expLibrarySize",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expLibrarySize=value)
		})
#########
setGeneric(
		name="contrLibrarySize",
		def=function(object){
			standardGeneric("contrLibrarySize")
		}
)
setGeneric(
		name="contrLibrarySize<-",
		def=function(object,value){
			standardGeneric("contrLibrarySize<-")
		}
)
setMethod("contrLibrarySize",
		signature(object="r3Cseq"),
		function(object){
			object@contrLibrarySize
		}
)
setReplaceMethod(
		f="contrLibrarySize",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrLibrarySize=value)
		})
##########
setGeneric(
		name="isControlInvolved",
		def=function(object){
			standardGeneric("isControlInvolved")
		})
setMethod("isControlInvolved",
		signature(object = "r3Cseq"),
		function (object){
			object@isControlInvolved
		}
)
##########
setGeneric(
		name="isBamInputFile",
		def=function(object){
			standardGeneric("isBamInputFile")
		})
setMethod("isBamInputFile",
		signature(object = "r3Cseq"),
		function (object){
			object@isBamInputFile
		}
)
##########
setGeneric(
		name="expReadLength",
		def=function(object){
			standardGeneric("expReadLength")
		}
)
setGeneric(
		name="expReadLength<-",
		def=function(object,value){
			standardGeneric("expReadLength<-")
		}
)
setMethod("expReadLength",
		signature(object = "r3Cseq"),
		function (object){
			object@expReadLength
		}
)
setReplaceMethod(
		f="expReadLength",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expReadLength=value)
		})
##########
setGeneric(
		name="contrReadLength",
		def=function(object){
			standardGeneric("contrReadLength")
		}
)
setGeneric(
		name="contrReadLength<-",
		def=function(object,value){
			standardGeneric("contrReadLength<-")
		}
)
setMethod("contrReadLength",
		signature(object = "r3Cseq"),
		function (object){
			object@contrReadLength
		}
)
setReplaceMethod(
		f="contrReadLength",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrReadLength=value)
		})
##########
setGeneric(
		name="expReadCount",
		def=function(object){
			standardGeneric("expReadCount")
		}
)
setGeneric(
		name="expReadCount<-",
		def=function(object,value){
			standardGeneric("expReadCount<-")
		}
)
setMethod("expReadCount",
		signature(object = "r3Cseq"),
		function (object){
			object@expReadCount
		}
)
setReplaceMethod(
		f="expReadCount",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expReadCount=value)
		})
#########
setGeneric(
		name="contrReadCount",
		def=function(object){
			standardGeneric("contrReadCount")
		}
)
setGeneric(
		name="contrReadCount<-",
		def=function(object,value){
			standardGeneric("contrReadCount<-")
		}
)
setMethod("contrReadCount",
		signature(object = "r3Cseq"),
		function (object){
			object@contrReadCount
		}
)
setReplaceMethod(
		f="contrReadCount",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrReadCount=value)
		})

###########
setGeneric(
		name="expRPM",
		def=function(object){
			standardGeneric("expRPM")
		}
)
setGeneric(
		name="expRPM<-",
		def=function(object,value){
			standardGeneric("expRPM<-")
		}
)
setMethod("expRPM",
		signature(object = "r3Cseq"),
		function (object){
			object@expRPM
		}
)
setReplaceMethod(
		f="expRPM",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expRPM=value)
		})
###########
setGeneric(
		name="contrRPM",
		def=function(object){
			standardGeneric("contrRPM")
		}
)
setGeneric(
		name="contrRPM<-",
		def=function(object,value){
			standardGeneric("contrRPM<-")
		}
)
setMethod("contrRPM",
		signature(object = "r3Cseq"),
		function (object){
			object@contrRPM
		}
)
setReplaceMethod(
		f="contrRPM",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrRPM=value)
		})
############
setGeneric(
		name="expInteractionRegions",
		def=function(object){
			standardGeneric("expInteractionRegions")
		}
)
setGeneric(
		name="expInteractionRegions<-",
		def=function(object,value){
			standardGeneric("expInteractionRegions<-")
		}
)
setMethod("expInteractionRegions",
		signature(object = "r3Cseq"),
		function (object){
			object@expInteractionRegions
		}
)
setReplaceMethod(
		f="expInteractionRegions",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expInteractionRegions=value)
		})

########
setGeneric(
		name="contrInteractionRegions",
		def=function(object){
			standardGeneric("contrInteractionRegions")
		}
)
setGeneric(
		name="contrInteractionRegions<-",
		def=function(object,value){
			standardGeneric("contrInteractionRegions<-")
		}
)
setMethod("contrInteractionRegions",
		signature(object = "r3Cseq"),
		function (object){
			object@contrInteractionRegions
		}
)
setReplaceMethod(
		f="contrInteractionRegions",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrInteractionRegions=value)
		})
#########
setGeneric(
		name="expCoverage",
		def=function(object){
			standardGeneric("expCoverage")
		}
)
setGeneric(
		name="expCoverage<-",
		def=function(object,value){
			standardGeneric("expCoverage<-")
		})

setMethod("expCoverage",
		signature(object = "r3Cseq"),
		function (object){
			object@expCoverage
		}
)
setReplaceMethod(
		f="expCoverage",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expCoverage=value)
		})

###########
setGeneric(
		name="contrCoverage",
		def=function(object){
			standardGeneric("contrCoverage")
		}
)
setGeneric(
		name="contrCoverage<-",
		def=function(object,value){
			standardGeneric("contrCoverage<-")
		})
setMethod("contrCoverage",
		signature(object = "r3Cseq"),
		function (object){
			object@contrCoverage
		}
)
setReplaceMethod(
		f="contrCoverage",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrCoverage=value)
		})


