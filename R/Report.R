# TODO: The function below is using to generate the final report.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@bccs.uib.no 
###############################################################################
setGeneric(
		name="generate3CseqReport",
		def=function(object){
			standardGeneric("generate3CseqReport")
		}
)

setMethod("generate3CseqReport",
		signature(object = "r3Cseq"),
		function (object){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			
			if(isControlInvolved(object)==FALSE){
				expInteractions<-expInteractionRegions(object)
				expLabeled<-expLabel(object)
				viewpoint<-getViewpoint(object)
				if(nrow(expInteractions)>0){
					
					file_name<-paste(expLabeled,".pdf",sep="")
					pdf(file=file_name, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
					plotOverviewInteractions(object)
					plotInteractionsNearViewpoint(object)
					plotInteractionsPerChromosome(object,space(viewpoint))
					plot3Cecdf(object)
					dev.off()
					exportInteractions2text(object)
					export3Cseq2bedGraph(object)
					print("Three files are generated : a pdf file of plots, a text file of interaction regions, and a bedGraph file.")
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
			}
			if(isControlInvolved(object)==TRUE){
				expInteractions  <-expInteractionRegions(object)
				contrInteractions <- contrInteractionRegions(object)
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				viewpoint<-getViewpoint(object)
				if(nrow(expInteractions)>0){
					file_name<-paste(expLabeled,"_",controlLabeled,".pdf",sep="")
					pdf(file=file_name, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
					plotOverviewInteractions(object)
					plotInteractionsNearViewpoint(object)
					plotInteractionsPerChromosome(object,space(viewpoint))
					plot3Cecdf(object)
					dev.off()
					exportInteractions2text(object)
					export3Cseq2bedGraph(object)
					print("Three files are generated : a pdf file of plots, a text file of interaction regions, and a bedGraph file.")
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}	
			}
		}
)



