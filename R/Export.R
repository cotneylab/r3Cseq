# TODO: The following functions are implemented for exporting the output of 3C-seq results.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@bccs.uib.no 
###############################################################################
setGeneric(
		name="export3Cseq2bedGraph",
		def=function(object,datatype = c("rpm", "raw_read")){
			standardGeneric("export3Cseq2bedGraph")
		}

)

setMethod("export3Cseq2bedGraph",
		signature(object = "r3Cseq"),
		function (object,datatype){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			datatype <- match.arg(datatype)
			
			if(isControlInvolved(object)==FALSE){
				expRPMs   <-expRPM(object)
				expLabeled<-expLabel(object)
				
				if(datatype=="rpm"){	
					if(nrow(expRPMs)>0){
						export.f <-data.frame(space=space(expRPMs),start=start(expRPMs),end=end(expRPMs),score=expRPMs$expRPMs)
						export.iranges<-RangedData(space=export.f$space,IRanges(start=export.f$start,end=export.f$end),score=export.f$score)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=expLabeled,description="reads per million")	
						file_name<-paste(expLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
				}else if(datatype=="raw_read"){
					if(nrow(expRPMs)>0){
						export.f <-data.frame(space=space(expRPMs),start=start(expRPMs),end=end(expRPMs),score=expRPMs$reads)
						export.iranges<-RangedData(space=export.f$space,IRanges(start=export.f$start,end=export.f$end),score=export.f$score)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=expLabeled,description="raw reads")	
						file_name<-paste(expLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No read found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
					
				}else{
					stop("Choose the input datatype parameter: 'raw_reads' or 'rpm' (reads per million)")
				}
			}
			if(isControlInvolved(object)==TRUE){
				expRPMs   <-expRPM(object)
				contrRPMs <-contrRPM(object)
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				
				if(datatype=="rpm"){

					if(nrow(expRPMs)>0){
						export.iranges<-RangedData(space=space(expRPMs),IRanges(start=start(expRPMs),end=end(expRPMs)),score=expRPMs$expRPMs)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=expLabeled,description="reads per million")	
						file_name<-paste(expLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
						

						export.iranges<-RangedData(space=space(contrRPMs),IRanges(start=start(contrRPMs),end=end(contrRPMs)),score=contrRPMs$contrRPMs)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=controlLabeled,description="reads per million")	
						file_name<-paste(controlLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))	
						
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
				}else if(datatype=="raw_read"){
					if(nrow(expRPMs)>0){
						export.iranges<-RangedData(space=space(expRPMs),IRanges(start=start(expRPMs),end=end(expRPMs)),score=expRPMs$expReads)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=expLabeled,description="raw reads")	
						file_name<-paste(expLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					
						export.iranges<-RangedData(space=space(contrRPMs),IRanges(start=start(contrRPMs),end=end(contrRPMs)),score=contrRPMs$contrReads)
						export.iranges<-export.iranges[export.iranges$score>0,]
						export.ucsc<-as(export.iranges,"UCSCData")
						export.ucsc@trackLine <- new("BasicTrackLine",name=controlLabeled,description="raw reads")	
						file_name<-paste(controlLabeled,".bedGraph",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))	
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}					
				}else{
					stop("Choose the input datatype parameter: 'raw_reads' or 'rpm' (reads per million)")
				}
				
			}
		}
)


setGeneric(
		name="exportInteractions2text",
		def=function(object){
			standardGeneric("exportInteractions2text")
		}
)

setMethod("exportInteractions2text",
		signature(object = "r3Cseq"),
		function (object){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			
			if(isControlInvolved(object)==FALSE){
				expInteractions<-expInteractionRegions(object)
				expLabeled<-expLabel(object)
				
					if(nrow(expInteractions)>0){
						
						file_name<-paste(expLabeled,".txt",sep="")
						export.data<-data.frame(chromosome=space(expInteractions),
								start=start(expInteractions),end=end(expInteractions),
								reads=expInteractions$reads,rpm=expInteractions$expRPMs,
								p_value=expInteractions$p_value)
						export.data<-export.data[order(export.data[,6]),]
						write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
					}
			}
			if(isControlInvolved(object)==TRUE){
				expInteractions  <-expInteractionRegions(object)
				contrInteractions <- contrInteractionRegions(object)
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				
				if(nrow(expInteractions)>0){
					
					file_name<-paste(expLabeled,".txt",sep="")
					export.data<-data.frame(chromosome=space(expInteractions),
							start=start(expInteractions),end=end(expInteractions),
							expReads=expInteractions$expReads,contrReads=expInteractions$contrReads,
							expRPMs=expInteractions$expRPMs,contrRPMs=expInteractions$contrRPMs,
							p_value=expInteractions$p_value,fold_change=expInteractions$fold_change)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
				if(nrow(contrInteractions)>0){
					
					file_name<-paste(controlLabeled,".txt",sep="")
					export.data<-data.frame(chromosome=space(contrInteractions),
							start=start(contrInteractions),end=end(contrInteractions),
							expReads=contrInteractions$expReads,contrReads=contrInteractions$contrReads,
							expRPMs=contrInteractions$expRPMs,contrRPMs=contrInteractions$contrRPMs,
							p_value=contrInteractions$p_value,fold_change=contrInteractions$fold_change)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)	
					export.data<-export.data[order(export.data[,8]),]
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
			}
		}
)


