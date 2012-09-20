# TODO: These following functions were implemented to get interaction regions.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@bccs.uib.no 

########1. get coverage function########
##getCoverage function start from reading in an aligned read file from different type of formats 
##for example the "export" file from Eland and the "BAM" file from Bowtie.
##We use the readAligned function from the Shortread package to read in the mapped reads. 
##getCoverage function was also implemented for reading in a BAM file by using Scanbam function in Rsamtools.

setGeneric(
		name="getCoverage",
		def=function(object){
			standardGeneric("getCoverage")
		}
)
setMethod("getCoverage",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
			}
			#############Loading the Bsgenome package######################
			if('BSgenome.Hsapiens.UCSC.hg19' %in% loadedNamespaces()==TRUE){
				detach(package:BSgenome.Hsapiens.UCSC.hg19,unload=TRUE)
			}
			if('BSgenome.Hsapiens.UCSC.hg18' %in% loadedNamespaces()==TRUE){
				detach(package:BSgenome.Hsapiens.UCSC.hg18,unload=TRUE)
			}
			if('BSgenome.Mmusculus.UCSC.mm9' %in% loadedNamespaces()==TRUE){
				detach(package:BSgenome.Mmusculus.UCSC.mm9,unload=TRUE)
			}
				
			objName <- deparse(substitute(object))
			
			if(isBamInputFile(object)==FALSE){	
				
				if(isControlInvolved(object)==TRUE){
					
					exp.input.file  <- alignedReadsExpFile(object)
					contr.input.file <- alignedReadsContrFile(object)
					
					if(file.exists(exp.input.file)==FALSE){
						stop(paste("Couldn't find file","-->",exp.input.file))
					}
					if(file.exists(contr.input.file)==FALSE){
						stop(paste("Couldn't find file","-->",contr.input.file))
					}
					
					expLabeled<-expLabel(object)
					controlLabeled<-contrLabel(object)
					
					if(nchar(expLabeled)==0 & nchar(controlLabeled)==0){
						stop("The package requires input of expLabel and contrLabel.")
					}
					
					print("start reading in ......")
					input.format	<- alignedReadsType(object)
					filt<-alignQualityFilter(10)
					exp.reads.raw   <- readAligned(exp.input.file,filter=filt,type=input.format)
					contr.reads.raw <- readAligned(contr.input.file,filter=filt,type=input.format)
				
					exp.read.map    <- exp.reads.raw[!is.na(position(exp.reads.raw))]
					contr.read.map  <- contr.reads.raw[!is.na(position(contr.reads.raw))]
			
					expLibrarySize(object)   <-length(exp.read.map)
					contrLibrarySize(object) <-length(contr.read.map)
				
					exp.read.length   <- width(sread(exp.read.map[1]))
					contr.read.length <- width(sread(contr.read.map[1]))
					exp.read.cov      <- RleList()	
					contr.read.cov    <- RleList()	
					orgName  <- organismName(object)
				
					print("making coverage vector......")
					if(orgName=="hg18"){
						library(BSgenome.Hsapiens.UCSC.hg18)
						exp.read.cov   <- coverage(exp.read.map,width=seqlengths(Hsapiens),extend=as.integer(exp.read.length))
						contr.read.cov <- coverage(contr.read.map,width=seqlengths(Hsapiens),extend=as.integer(contr.read.length))
					}else if(orgName=="hg19"){
						library(BSgenome.Hsapiens.UCSC.hg19)
						exp.read.cov   <- coverage(exp.read.map,width=seqlengths(Hsapiens),extend=as.integer(exp.read.length))
						contr.read.cov <- coverage(contr.read.map,width=seqlengths(Hsapiens),extend=as.integer(contr.read.length))
					}else if(orgName =="mm9"){
						library(BSgenome.Mmusculus.UCSC.mm9)
						exp.read.cov   <- coverage(exp.read.map,width=seqlengths(Mmusculus),extend=as.integer(exp.read.length))
						contr.read.cov <- coverage(contr.read.map,width=seqlengths(Mmusculus),extend=as.integer(contr.read.length))
					}else{
						stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
					}
					expReadLength(object)<-exp.read.length
					contrReadLength(object)<-contr.read.length
					expCoverage(object)<-exp.read.cov
					contrCoverage(object)<-contr.read.cov
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("making coverage is done.")
				}
			
				if(isControlInvolved(object)==FALSE){
					exp.input.file  <- alignedReadsExpFile(object)
					if(file.exists(exp.input.file)==FALSE){
						stop(paste("Couldn't find file","-->",exp.input.file))
					}
					expLabeled<-expLabel(object)
					if(nchar(expLabeled)==0){
						stop("The package requires input of expLabel.")
					}
					print("start reading in ......")
					input.format	<- alignedReadsType(object)
					filt<-alignQualityFilter(10)
					exp.reads.raw   <- readAligned(exp.input.file,filter=filt,type=input.format)
					exp.read.map    <- exp.reads.raw[!is.na(position(exp.reads.raw))]
					expLibrarySize(object)   <-length(exp.read.map)
					exp.read.length <- width(sread(exp.read.map[1]))
					exp.read.cov    <- RleList()	
					orgName  <- organismName(object)
					print("making coverage vector......")
					if(orgName=="hg18"){
						library(BSgenome.Hsapiens.UCSC.hg18)
						exp.read.cov <- coverage(exp.read.map,width=seqlengths(Hsapiens),extend=as.integer(exp.read.length))
					}else if(orgName=="hg19"){
						library(BSgenome.Hsapiens.UCSC.hg19)
						exp.read.cov <- coverage(exp.read.map,width=seqlengths(Hsapiens),extend=as.integer(exp.read.length))
					}else if(orgName =="mm9"){
						library(BSgenome.Mmusculus.UCSC.mm9)
						exp.read.cov   <- coverage(exp.read.map,width=seqlengths(Mmusculus),extend=as.integer(exp.read.length))
					}else{
						stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
					}
					expReadLength(object)<-exp.read.length
					expCoverage(object)  <-exp.read.cov
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("making coverage is done.")
				}
			}
			if(isBamInputFile(object)==TRUE){
				
				if(isControlInvolved(object)==TRUE){
					
					exp.bam.file    <-alignedReadsBamExpFile(object)
					contr.bam.file  <-alignedReadsBamContrFile(object)
					
					if(file.exists(exp.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",exp.bam.file))
					}			
					if(file.exists(contr.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",contr.bam.file))
					}	
					
					expLabeled<-expLabel(object)
					controlLabeled<-contrLabel(object)
					
					if(nchar(expLabeled)==0 & nchar(controlLabeled)==0){
						stop("The package requires input of expLabel and contrLabel.")
					}
					
					print("start reading in ......")
										
					what    <- c("rname", "strand", "pos", "qwidth","seq","qual")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					exp.bam <- scanBam(exp.bam.file, param=param)
					contr.bam <- scanBam(contr.bam.file, param=param)
					
					exp.qual    <-QualityScaledBStringSet(exp.bam[[1]]$seq, exp.bam[[1]]$qual)
					exp.read.length <-width(exp.qual)[1]
					exp.qa   <- FastqQuality(quality(exp.qual))
					exp.qa   <- as(exp.qa, "matrix")
					exp.qa.avg   <-as.integer(rowMeans(exp.qa))	
					exp.GRanges<-GRanges(seqnames=as.vector(exp.bam[[1]]$rname),IRanges(start=exp.bam[[1]]$pos,width=exp.read.length),strand=exp.bam[[1]]$strand,qual=exp.qa.avg)	
					exp.GRanges.filtered <-exp.GRanges[elementMetadata(exp.GRanges)$qual >=10]
					
					
					contr.qual    <-QualityScaledBStringSet(contr.bam[[1]]$seq, contr.bam[[1]]$qual)
					contr.read.length <-width(contr.qual)[1]
					contr.qa   <- FastqQuality(quality(contr.qual))
					contr.qa   <- as(contr.qa, "matrix")
					contr.qa.avg   <-as.integer(rowMeans(contr.qa))	
					contr.GRanges<-GRanges(seqnames=as.vector(contr.bam[[1]]$rname),IRanges(start=contr.bam[[1]]$pos,width=contr.read.length),strand=contr.bam[[1]]$strand,qual=contr.qa.avg)	
					contr.GRanges.filtered <-contr.GRanges[elementMetadata(contr.GRanges)$qual >=10]
					
					expLibrarySize(object)   <-length(exp.GRanges.filtered )
					contrLibrarySize(object) <-length(contr.GRanges.filtered )
				
					exp.read.cov      <- RleList()	
					contr.read.cov    <- RleList()	
					orgName  <- organismName(object)
					
					print("making coverage vector......")
					if(orgName=="hg18"){		
						library(BSgenome.Hsapiens.UCSC.hg18)
						hg18.chromlens <- seqlengths(Hsapiens)
						seqlevels(exp.GRanges.filtered)<-names(hg18.chromlens)
						seqlevels(contr.GRanges.filtered)<-names(hg18.chromlens)	
						seqlengths(exp.GRanges.filtered)<-hg18.chromlens
						seqlengths(contr.GRanges.filtered)<-hg18.chromlens
						exp.read.cov   <- coverage(exp.GRanges.filtered)
						contr.read.cov <- coverage(contr.GRanges.filtered)
					}else if(orgName=="hg19"){
						library(BSgenome.Hsapiens.UCSC.hg19)
						hg19.chromlens <- seqlengths(Hsapiens)
						seqlevels(exp.GRanges.filtered)<-names(hg19.chromlens)
						seqlevels(contr.GRanges.filtered)<-names(hg19.chromlens)	
						seqlengths(exp.GRanges.filtered)<-hg19.chromlens
						seqlengths(contr.GRanges.filtered)<-hg19.chromlens		
						exp.read.cov   <- coverage(exp.GRanges.filtered)
						contr.read.cov <- coverage(contr.GRanges.filtered)
					}else if(orgName =="mm9"){
						library(BSgenome.Mmusculus.UCSC.mm9)
						mm9.chromlens <- seqlengths(Mmusculus)
						seqlevels(exp.GRanges.filtered)<-names(mm9.chromlens)
						seqlevels(contr.GRanges.filtered)<-names(mm9.chromlens)			
						seqlengths(exp.GRanges.filtered)<-mm9.chromlens
						seqlengths(contr.GRanges.filtered)<-mm9.chromlens			
						exp.read.cov   <- coverage(exp.GRanges.filtered)
						contr.read.cov <- coverage(contr.GRanges.filtered)
					}else{
						stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
					}
					expReadLength(object)<-exp.read.length
					contrReadLength(object)<-contr.read.length
					expCoverage(object)<-exp.read.cov
					contrCoverage(object)<-contr.read.cov
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("making coverage is done.")
				}
				if(isControlInvolved(object)==FALSE){
					exp.bam.file    <-alignedReadsBamExpFile(object)
					if(file.exists(exp.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",exp.bam.file))
					}			
					
					expLabeled<-expLabel(object)	
					
					if(nchar(expLabeled)==0){
						stop("The package requires input of expLabel.")
					}
					
					print("start reading in ......")
					
					what    <- c("rname", "strand", "pos", "qwidth","seq","qual")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					exp.bam <- scanBam(exp.bam.file, param=param)					
					exp.qual    <-QualityScaledBStringSet(exp.bam[[1]]$seq, exp.bam[[1]]$qual)
					exp.read.length <-width(exp.qual)[1]
					exp.qa   <- FastqQuality(quality(exp.qual))
					exp.qa   <- as(exp.qa, "matrix")
					exp.qa.avg   <-as.integer(rowMeans(exp.qa))	
					exp.GRanges<-GRanges(seqnames=as.vector(exp.bam[[1]]$rname),IRanges(start=exp.bam[[1]]$pos,width=exp.read.length),strand=exp.bam[[1]]$strand,qual=exp.qa.avg)	
					exp.GRanges.filtered <-exp.GRanges[elementMetadata(exp.GRanges)$qual >=10]
					
					expLibrarySize(object)   <-length(exp.GRanges.filtered )			
					exp.read.cov      <- RleList()
					orgName  <- organismName(object)
					
					print("making coverage vector......")
					if(orgName=="hg18"){		
						library(BSgenome.Hsapiens.UCSC.hg18)
						hg18.chromlens <- seqlengths(Hsapiens)
						seqlevels(exp.GRanges.filtered)<-names(hg18.chromlens)
						seqlengths(exp.GRanges.filtered)<-hg18.chromlens		
						exp.read.cov   <- coverage(exp.GRanges.filtered)
					}else if(orgName=="hg19"){
						library(BSgenome.Hsapiens.UCSC.hg19)
						hg19.chromlens <- seqlengths(Hsapiens)
						seqlevels(exp.GRanges.filtered)<-names(hg19.chromlens)
						seqlengths(exp.GRanges.filtered)<-hg19.chromlens			
						exp.read.cov   <- coverage(exp.GRanges.filtered)
					}else if(orgName =="mm9"){
						library(BSgenome.Mmusculus.UCSC.mm9)
						mm9.chromlens <- seqlengths(Mmusculus)
						seqlevels(exp.GRanges.filtered)<-names(mm9.chromlens)
						seqlengths(exp.GRanges.filtered)<-mm9.chromlens		
						exp.read.cov   <- coverage(exp.GRanges.filtered)
					}else{
						stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
					}
					expReadLength(object)<-exp.read.length
					expCoverage(object)<-exp.read.cov
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("making coverage is done.")
					
				}
			}
		}
)
########2. getReadCountPerRestrictionFragment########
##This function is used to count number############## 
##of reads per input restriction fragment.#################
##Basically, we sum up reads from the forward and the reverse strand##
##and we then calculate reads per million per length of each restriction fragment##
setGeneric(
		name="getReadCountPerRestrictionFragment",
		def=function(object){
			standardGeneric("getReadCountPerRestrictionFragment")
		}
)

setMethod("getReadCountPerRestrictionFragment",
		signature(object = "r3Cseq"),
		function (object){
			objName <- deparse(substitute(object))
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				getCoverage(object)
				enzymeDb	<-new("repbaseEnzyme")
				expCoverage <-expCoverage(object)
				resEnzyme 	<-restrictionEnzyme(object)	
				orgName  <- organismName(object)
				
				if(orgName=="hg18"){		
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
				}else if(orgName=="hg19"){
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
				}else if(orgName =="mm9"){
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,19),'X','Y'),sep='')]
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}	
				
				readCounts	<- RangedData()
				print("start counting number of reads per each restriction fragment.......")
				for(i in 1:length(expCoverage)){
					if(sum(runValue(expCoverage[[i]])) > 0){
						fragments<-getRestrictionFragments(enzymeDb,resEnzyme,orgName,names(expCoverage[i]))
						windows.plus<- Views(expCoverage[[i]],fragments$start,fragments$start+expReadLength(object)+1)
						windows.minus<- Views(expCoverage[[i]],fragments$end-expReadLength(object),fragments$end)
						n.reads<-(viewMaxs(windows.minus)+viewMaxs(windows.plus))
						readCount.per.fragments<-RangedData(space=names(expCoverage[i]),IRanges(start=fragments$start,end=fragments$end),reads=n.reads)
						readCount.per.fragments.filter<-readCount.per.fragments[readCount.per.fragments$reads>0,]
						readCounts<-c(readCounts,readCount.per.fragments.filter)
						print(paste(names(expCoverage[i]),"--->","is done!"))
					}
				}
				expReadCount(object)<-readCounts
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
			if(isControlInvolved(object)==TRUE){
				getCoverage(object)
				enzymeDb	  <-new("repbaseEnzyme")
				expCoverage   <-expCoverage(object)
				contrCoverage <-contrCoverage(object)
				orgName		  <-organismName(object)
				resEnzyme 	  <-restrictionEnzyme(object)
								
				if(orgName=="hg18"){		
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
					contrCoverage<-contrCoverage[names(contrCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
				}else if(orgName=="hg19"){
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
					contrCoverage<-contrCoverage[names(contrCoverage) %in% paste('chr',c(seq(1,22),'X','Y'),sep='')]
				}else if(orgName =="mm9"){
					expCoverage<-expCoverage[names(expCoverage) %in% paste('chr',c(seq(1,19),'X','Y'),sep='')]
					contrCoverage<-contrCoverage[names(contrCoverage) %in% paste('chr',c(seq(1,19),'X','Y'),sep='')]
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}
				
				#######################################################
				print("start counting number of reads per each restriction fragment.......")
				########calculate read count for the experiment########
				expReadCounts <- RangedData()

				for(i in c(names(expCoverage))){
				    if(sum(runValue(expCoverage[[i]])) > 0){
						if(i %in% names(contrCoverage)){
						
						fragments 	  <- getRestrictionFragments(enzymeDb,resEnzyme,orgName,i)
						
						exp.windows.plus<- Views(expCoverage[[i]],fragments$start,fragments$start+expReadLength(object)+1)
						exp.windows.minus<- Views(expCoverage[[i]],fragments$end-expReadLength(object),fragments$end)
						exp.n.reads<-(viewMaxs(exp.windows.minus)+viewMaxs(exp.windows.plus))
						
						contr.windows.plus<- Views(contrCoverage[[i]],fragments$start,fragments$start+expReadLength(object)+1)
						contr.windows.minus<- Views(contrCoverage[[i]],fragments$end-expReadLength(object),fragments$end)
						contr.n.reads<-(viewMaxs(contr.windows.minus)+viewMaxs(contr.windows.plus))
					
						readCount.per.fragments<-RangedData(space=i,
								IRanges(start=fragments$start,end=fragments$end),
								expReads=exp.n.reads,
								contrReads=contr.n.reads
						)
						readCount.per.fragments.filter<-readCount.per.fragments[readCount.per.fragments$expReads>0,]
						expReadCounts<-c(expReadCounts,readCount.per.fragments.filter)				
					}
					print(paste(i,"--->","in the experiment is done!"))
				  }
				}
				expReadCount(object)<-expReadCounts
				
				#########calculate read count for the control########
				contrReadCounts <-RangedData()
				for(i in c(names(contrCoverage))){
				   if(sum(runValue(contrCoverage[[i]])) > 0){
					  if(i %in% names(expCoverage)){
						
						fragments 	  <- getRestrictionFragments(enzymeDb,resEnzyme,orgName,i)
						contr.windows.plus<- Views(contrCoverage[[i]],fragments$start,fragments$start+expReadLength(object)+1)
						contr.windows.minus<- Views(contrCoverage[[i]],fragments$end-expReadLength(object),fragments$end)
						contr.n.reads<-(viewMaxs(contr.windows.minus)+viewMaxs(contr.windows.plus))
						
						exp.windows.plus<- Views(expCoverage[[i]],fragments$start,fragments$start+expReadLength(object)+1)
						exp.windows.minus<- Views(expCoverage[[i]],fragments$end-expReadLength(object),fragments$end)
						exp.n.reads<-(viewMaxs(exp.windows.minus)+viewMaxs(exp.windows.plus))
						
						readCount.per.fragments<-RangedData(space=i,
								IRanges(start=fragments$start,end=fragments$end),
								expReads=exp.n.reads,
								contrReads=contr.n.reads
						)
						readCount.per.fragments.filter<-readCount.per.fragments[readCount.per.fragments$contrReads>0,]
						contrReadCounts <- c(contrReadCounts,readCount.per.fragments.filter)
					  }
					print(paste(i,"--->"," in the control is done!"))
					}
				}
				
				contrReadCount(object)<-contrReadCounts
				
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
		}
)
########3. getViewpoint#############################
##We defined the restriction fragment that contains the highest number of reads as the viewpoint##
##The viewpoint could be the promoter of the interested gene or a transcription factor binding site.

setGeneric(
		name="getViewpoint",
		def=function(object){
			standardGeneric("getViewpoint")
		}
)

setMethod("getViewpoint",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				
				readCounts  <-expReadCount(object)
				n.reads<-c()
				for(i in c(names(readCounts))){
					s.readcounts<-readCounts[space(readCounts)==i,]
					n.reads.chr<-sum(s.readcounts$reads)
					names(n.reads.chr)<-i
					n.reads <-append(n.reads,n.reads.chr)
				}
				####get the maximum reads###
				n.max<-max(n.reads)
				viewpoint.chr<-names(which(n.reads==n.max))
				
				viewpoint.index<-which(readCounts$reads==max(readCounts$reads))
				viewpoint<-readCounts[viewpoint.index,]
				viewpoint<-viewpoint[space(viewpoint)==viewpoint.chr,]
				if(nrow(viewpoint)==1){
					return(viewpoint)
				}else if(nrow(viewpoint) >1){
					viewpoint<-viewpoint[1,]
					return(viewpoint)
				}
			}
			if(isControlInvolved(object)==TRUE){
				readCounts  <-expReadCount(object)
				n.reads<-c()
				for(i in c(names(readCounts))){
					s.readcounts<-readCounts[space(readCounts)==i,]
					n.reads.chr<-sum(s.readcounts$expReads)
					names(n.reads.chr)<-i
					n.reads <-append(n.reads,n.reads.chr)
				}
				####get the maximum reads###
				n.max<-max(n.reads)
				viewpoint.chr<-names(which(n.reads==n.max))
				
				viewpoint.index<-which(readCounts$expReads==max(readCounts$expReads))
				viewpoint<-readCounts[viewpoint.index,]
				viewpoint<-viewpoint[space(viewpoint)==viewpoint.chr,]
				if(nrow(viewpoint)==1){
					return(viewpoint)
				}else if(nrow(viewpoint) >1){
					viewpoint<-viewpoint[1,]
					return(viewpoint)
				}
			}
			
		}
)

########4. calculateRPM########
##Calculate read per million per for each restriction fragment

setGeneric(
		name="calculateRPM",
		def=function(object){
			standardGeneric("calculateRPM")
		}
)

setMethod("calculateRPM",
		signature(object = "r3Cseq"),
		function (object){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				expReadCounts   <-expReadCount(object)
				if(nrow(expReadCounts)>0){
					objName 	 <- deparse(substitute(object))
					expLibSize 	 <- expLibrarySize(object)
					seqDepth <-10^6
					if(expLibSize >=10^6){
						seqDepth <- 10^6
					}else if(expLibSize >=10^5 && expLibSize< 10^6){
						seqDepth <- 10^5
					}else if(expLibSize >=10^4 && expLibSize < 10^5 ){
						seqDepth <- 10^4
					}else if(expLibSize >=10^3 && expLibSize < 10^4){
						seqDepth <-10^3
					}else{
						stop("This is too low sequencing depth!!!!!. We have to stop the analysis.")
					}
					
					expRPMs   	 <- round(expReadCounts$reads / (expLibSize/seqDepth))
					expReadCounts$expRPMs <- expRPMs
					expRPM(object) <-expReadCounts
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print(paste("RPM calculation is done. Use function 'expRPM' or 'contrRPM' to get the result."))
					
				}else{
					stop("No reads found in the r3Cseq object, you have to run 'getReadCountPerRestrictionFragment' in order to get the count of reads.")
				}
			}
			if(isControlInvolved(object)==TRUE){
				expReadCounts   <-expReadCount(object)
				contrReadCounts <-contrReadCount(object)
				if(nrow(expReadCounts)>0){
					
					objName 	 <- deparse(substitute(object))
					
					expLibSize 	 <- expLibrarySize(object)
					contrLibSize <- contrLibrarySize(object)
					######assign reads per millions for the experiment and control###############
					seqDepth <-10^6
					if(expLibSize >=10^6){
						seqDepth <- 10^6
					}else if(expLibSize >=10^5 && expLibSize< 10^6){
						seqDepth <- 10^5
					}else if(expLibSize >=10^4 && expLibSize < 10^5 ){
						seqDepth <- 10^4
					}else if(expLibSize >=10^3 && expLibSize < 10^4){
						seqDepth <-10^3
					}else{
						stop("This is too low sequencing depth!!!!!. We have to stop the analysis.")
					}
					
					expRPMs   	 <- round(expReadCounts$expReads / (expLibSize/seqDepth))
					contrRPMs    <- round(expReadCounts$contrReads / (contrLibSize/seqDepth))
					
					expReadCounts$expRPMs <- expRPMs
					expReadCounts$contrRPMs <-contrRPMs
					
					expRPM(object) <-expReadCounts
					############################################################
					expRPMs   	<- round(contrReadCounts$expReads / (expLibSize/seqDepth))
					contrRPMs   <- round(contrReadCounts$contrReads / (contrLibSize/seqDepth))
					
					contrReadCounts$expRPMs  <- expRPMs
					contrReadCounts$contrRPMs <-contrRPMs
					contrRPM(object) <-contrReadCounts
					
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print(paste("RPM calculation is done. Use function 'expRPM' or 'contrRPM' to get the result."))	
				}else{
					stop("No reads found in the r3Cseq object, you have to run 'getReadCountPerRestrictionFragment' in order to get the count of reads.")
				}
			}
		}
)

########5. GetInteractions########
##Assign p-value by empirical cumulative distribution function 
##and calculate fold-change comparison with the control

setGeneric(
		name="getInteractions",
		def=function(object){
			standardGeneric("getInteractions")
		}
)

setMethod("getInteractions",
		signature(object = "r3Cseq"),
		function (object){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				expRPM   <-expRPM(object)
				if(nrow(expRPM)>0){
					objName <- deparse(substitute(object))
					p_value <- 1-ecdf(expRPM$expRPMs)(expRPM$expRPMs)
					expRPM$p_value <-p_value							
					expInteractionRegions(object) <-expRPM
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Calculation is done. Use function 'expInteractionRegions' or 'contrInteractionRegions' to get the result.")
				}else{
					stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
				}
			}
			if(isControlInvolved(object)==TRUE){
				expRPM   <-expRPM(object)
				contrRPM <-contrRPM(object)
				
				if(nrow(expRPM)>0){
					objName <- deparse(substitute(object))
					p_value <- 1-ecdf(expRPM$expRPMs)(expRPM$expRPMs)
					fold_change <- expRPM$expRPMs/(1+expRPM$contrRPMs)
					expRPM$p_value <- p_value
					expRPM$fold_change<-fold_change
					expInteractionRegions(object) <- expRPM
					
					##calculate in the control##
					p_value <- 1-ecdf(contrRPM$contrRPMs)(contrRPM$contrRPMs)
					fold_change <- contrRPM$contrRPMs/(1+contrRPM$expRPMs)
					contrRPM$p_value <- p_value
					contrRPM$fold_change <-fold_change
					contrInteractionRegions(object) <- contrRPM
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Calculation is done. Use function 'expInteractionRegions' or 'contrInteractionRegions' to get the result.")
				}else{
					stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")
				}
			}
		}
)
