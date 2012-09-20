# TODO: These following functions are implemented for visualizing 3C-seq data.
# Author: Supat Thongjuea
# Contact:supat.thongjuea@bccs.uib.no 

setGeneric(
		name="plotOverviewInteractions",
		def=function(object,cutoff.p_value=0.05,cutoff.fold_change=2){
			standardGeneric("plotOverviewInteractions")
		}
)
setMethod("plotOverviewInteractions",
		signature(object = "r3Cseq"),
		function (object,cutoff.p_value,cutoff.fold_change){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				orgName<-organismName(object)
				chr.data=c()
				if(orgName=="hg18"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Hsapiens)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName=="hg19"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Hsapiens)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="mm9"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Mmusculus)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}
				######check interactions############
				expInteractions <-expInteractionRegions(object)
				
				if(nrow(expInteractions) ==0){
					stop("There are no interaction regions found in r3Cseq object. Use 'getInteractions' function to get interaction regions")
				}
				
				viewpoint <-getViewpoint(object)
				viewpoint.exp.index<-which(start(expInteractions)==start(viewpoint) & end(expInteractions)==end(viewpoint))
				viewpoint.index<-which(expInteractions$reads==max(expInteractions$reads))
				
				exp.filted <-expInteractions[-(viewpoint.index),]
				exp.filted <-exp.filted[exp.filted$p_value <=cutoff.p_value,]
				
				if(nrow(exp.filted) ==0){
					stop("There are no interaction regions pass your input parameters.")
				}
				
				chr.size.max<-max(chr.data$size)
				max.scale <- floor(chr.size.max/10^6)
				box.x.size=chr.size.max+5e6
				plot(c(1,box.x.size), c(1,100), type= "n", ylab="",yaxt='n',
						xaxt='n',xlab="Chromosomal position (Mbp)",
						main=paste("3C-seq distribution of interaction regions (p-value <=",cutoff.p_value,")"))
				
				axis(1, at=(seq(0, max.scale*10^6, by=10*10^6)),lab=c(seq(0,max.scale,by=10)),cex.axis=0.8)
				
				expLabeled<-expLabel(object)
				polygon(c(chr.size.max-10e6-10e5,chr.size.max-10e6,chr.size.max-10e6+10e5),
						c(52,50,52), col="red")
				text(chr.size.max-2e6,51,"  viewpoint",cex = .8)
				
				y=0;
				
				for(i in 1:nrow(chr.data)){
					y=y+4
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="grey",lwd=16)
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="black")
					text.x<-chr.data$size[i]+5e6
					text(text.x, y+1, chr.data$name[i],cex = .8)
				}
				
				start.exp <-seq(from=3,to=100,by=4)
				end.exp   <-seq(from=6,to=100,by=4)
				
				exp.coor<-data.frame(i=start.exp[1:24],j=end.exp[1:24])
				
				exp.filted$rank<-rank(-exp.filted$expRPMs)
				exp.max.rank<-max(exp.filted$rank)
				
				exppalette<-rev(brewer.pal(7,"Reds"))
				exp.filted$col<-""
				exp.filted$col[exp.filted$rank <=round(0.01*exp.max.rank)]<-exppalette[1]
				exp.filted$col[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)]<-exppalette[2]
				exp.filted$col[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)]<-exppalette[3]
				exp.filted$col[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)]<-exppalette[4]
				exp.filted$col[exp.filted$rank > round(0.3*exp.max.rank) & round(0.4*exp.max.rank)]<-exppalette[5]
				exp.filted$col[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)]<-exppalette[6]
				exp.filted$col[exp.filted$rank > round(0.5*exp.max.rank)]<-exppalette[7]
				
				if(round(0.01*exp.max.rank)>0){
					sv1<-min(exp.filted$expRPMs[exp.filted$rank <=round(0.01*exp.max.rank)])
					sv2<-min(exp.filted$expRPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$expRPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$expRPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$expRPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$expRPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
				
					name1<-paste(">",sv1," RPMs")
					name2<-paste(sv2,"-",sv1," RPMs")
					name3<-paste(sv3,"-",sv2," RPMs")
					name4<-paste(sv4,"-",sv3," RPMs")
					name5<-paste(sv5,"-",sv4," RPMs")
					name6<-paste(sv6,"-",sv5," RPMs")
					name7<-paste("<",sv6," RPMs")
				}else{
					sv1<-max(exp.filted$expRPMs)
					sv2<-min(exp.filted$expRPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$expRPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$expRPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$expRPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$expRPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
					
					name1<-paste(">",sv1," RPMs")
					name2<-paste(sv2,"-",sv1," RPMs")
					name3<-paste(sv3,"-",sv2," RPMs")
					name4<-paste(sv4,"-",sv3," RPMs")
					name5<-paste(sv5,"-",sv4," RPMs")
					name6<-paste(sv6,"-",sv5," RPMs")
					name7<-paste("<",sv6," RPMs")
				}
				exp.names<-c(name1,name2,name3,name4,name5,name6,name7)
				
				
				viewpoint <-getViewpoint(object)
				i=0
				for (chri in 1:nrow(chr.data)){
					i=i+1
					
					if(as.character(chr.data$name[chri]) %in% names(exp.filted)){
						exp.chr<-exp.filted[space(exp.filted)==as.character(chr.data$name[chri]),]
						if(nrow(exp.chr) >0){
							rect(start(exp.chr),exp.coor$i[i], end(exp.chr), exp.coor$j[i], border=exp.chr$col)		
						}
						if(space(viewpoint)==as.character(chr.data$name[chri])){
							polygon(c(start(viewpoint)-10e5,start(viewpoint),start(viewpoint)+10e5),
								c(exp.coor$j[i]+1.5,exp.coor$j[i],exp.coor$j[i]+1.5), col="red")
						}
					}
				}
				legend("topright",legend = exp.names, fill=exppalette, cex=0.55,title=expLabeled,bty="n")
			}
			if(isControlInvolved(object)==TRUE){
				
				######check interactions######
				expInteractions <-expInteractionRegions(object)
				contrInteractions <-contrInteractionRegions(object)
				
				if(nrow(expInteractions) ==0){
					stop("There are no interaction regions found in r3Cseq object. Use 'getInteractions' function to get interaction regions")
				}
				viewpoint <-getViewpoint(object)
				
				viewpoint.exp.index<-which(start(expInteractions)==start(viewpoint) & end(expInteractions)==end(viewpoint))
				viewpoint.contr.index<-which(start(contrInteractions)==start(viewpoint) & end(contrInteractions)==end(viewpoint))
				
				exp.filted <-expInteractions[-(viewpoint.exp.index),]
				contr.filted<-contrInteractions[-(viewpoint.contr.index),]
				
				exp.filted <-exp.filted[exp.filted$p_value <=cutoff.p_value & exp.filted$fold_change >=cutoff.fold_change,]
				contr.filted <-contr.filted[contr.filted$p_value <=cutoff.p_value & contr.filted$fold_change >=cutoff.fold_change,]
				if(nrow(exp.filted) ==0){
					stop("There are no interaction regions pass your input parameters.")
				}
				
				#######draw chromosome########
				
				orgName<-organismName(object)
				chr.data=c()
				if(orgName=="hg18"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Hsapiens)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName=="hg19"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Hsapiens)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="mm9"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){	
						chr.size<-seqlengths(Mmusculus)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}
	
				chr.size.max<-max(chr.data$size)
				max.scale <- floor(chr.size.max/10^6)
				box.x.size=chr.size.max+5e6
				plot(c(1,box.x.size), c(1,100), type= "n", ylab="",yaxt='n',
						xaxt='n',xlab="Chromosomal position (Mbp)",
						main=paste("3C-seq distribution of interaction regions (p-value <=",cutoff.p_value,"and fold change >=",cutoff.fold_change,")"))
				
				axis(1, at=(seq(0, max.scale*10^6, by=10*10^6)),lab=c(seq(0,max.scale,by=10)),cex.axis=0.8)
				
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				polygon(c(chr.size.max-10e6-10e5,chr.size.max-10e6,chr.size.max-10e6+10e5),
						c(52,50,52), col="red")
				text(chr.size.max-2e6,51,"  viewpoint",cex = .8)
				
				y=0;
				
				for(i in 1:nrow(chr.data)){
					y=y+4
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="grey",lwd=16)
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="black")
					text.x<-chr.data$size[i]+5e6
					text(text.x, y+1, chr.data$name[i],cex = .8)
				}
				
				start.contr <-seq(from=3,to=100,by=4)
				end.contr   <-seq(from=4.5,to=100,by=4)
				
				start.exp <-seq(from=4.5,to=100,by=4)
				end.exp   <-seq(from=6,to=100,by=4)
				
				exp.coor<-data.frame(i=start.exp[1:24],j=end.exp[1:24])
				contr.coor<-data.frame(i=start.contr[1:24],j=end.contr[1:24])
				
				exp.filted$rank<-rank(-exp.filted$expRPMs)
				contr.filted$rank<-rank(-contr.filted$contrRPMs)
				
				exp.max.rank<-max(exp.filted$rank)
				contr.max.rank<-max(contr.filted$rank)
				
				exppalette<-rev(brewer.pal(7,"Greens"))
				exp.filted$col<-""
				exp.filted$col[exp.filted$rank <=round(0.01*exp.max.rank)]<-exppalette[1]
				exp.filted$col[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)]<-exppalette[2]
				exp.filted$col[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)]<-exppalette[3]
				exp.filted$col[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)]<-exppalette[4]
				exp.filted$col[exp.filted$rank > round(0.3*exp.max.rank) & round(0.4*exp.max.rank)]<-exppalette[5]
				exp.filted$col[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)]<-exppalette[6]
				exp.filted$col[exp.filted$rank > round(0.5*exp.max.rank)]<-exppalette[7]
				
				if(round(0.01*exp.max.rank)>0){
					sv1<-min(exp.filted$expRPMs[exp.filted$rank <=round(0.01*exp.max.rank)])
					sv2<-min(exp.filted$expRPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$expRPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$expRPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$expRPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$expRPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
	
					name1<-paste(">",sv1," RPMs")
					name2<-paste(sv2,"-",sv1," RPMs")
					name3<-paste(sv3,"-",sv2," RPMs")
					name4<-paste(sv4,"-",sv3," RPMs")
					name5<-paste(sv5,"-",sv4," RPMs")
					name6<-paste(sv6,"-",sv5," RPMs")
					name7<-paste("<",sv6," RPMs")
				}else{
					sv1<-max(exp.filted$expRPMs)
					sv2<-min(exp.filted$expRPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$expRPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$expRPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$expRPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$expRPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
					
					name1<-paste(">",sv1," RPMs")
					name2<-paste(sv2,"-",sv1," RPMs")
					name3<-paste(sv3,"-",sv2," RPMs")
					name4<-paste(sv4,"-",sv3," RPMs")
					name5<-paste(sv5,"-",sv4," RPMs")
					name6<-paste(sv6,"-",sv5," RPMs")
					name7<-paste("<",sv6," RPMs")
				}
				exp.names<-c(name1,name2,name3,name4,name5,name6,name7)
				
			
				contrpalette<-rev(brewer.pal(7,"Reds"))
				contr.filted$col<-""
				contr.filted$col[contr.filted$rank <=round(0.01*contr.max.rank)]<-contrpalette[1]
				contr.filted$col[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)]<-contrpalette[2]
				contr.filted$col[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)]<-contrpalette[3]
				contr.filted$col[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)]<-contrpalette[4]
				contr.filted$col[contr.filted$rank > round(0.3*contr.max.rank) & round(0.4*contr.max.rank)]<-contrpalette[5]
				contr.filted$col[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)]<-contrpalette[6]
				contr.filted$col[contr.filted$rank > round(0.5*contr.max.rank)]<-contrpalette[7]
				
				if(round(0.01*contr.max.rank)>0){
					csv1<-min(contr.filted$contrRPMs[contr.filted$rank <=round(0.01*contr.max.rank)])
					csv2<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)])
					csv3<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)])
					csv4<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)])
					csv5<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.3*contr.max.rank) & contr.filted$rank <=round(0.4*contr.max.rank)])
					csv6<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)])
					
					tname1<-paste(">",csv1," RPMs")
					tname2<-paste(csv2,"-",csv1," RPMs")
					tname3<-paste(csv3,"-",csv2," RPMs")
					tname4<-paste(csv4,"-",csv3," RPMs")
					tname5<-paste(csv5,"-",csv4," RPMs")
					tname6<-paste(csv6,"-",csv5," RPMs")
					tname7<-paste("<",csv6," RPMs")
				}else{
					csv1<-max(contr.filted$contrRPMs)
					csv2<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)])
					csv3<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)])
					csv4<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)])
					csv5<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.3*contr.max.rank) & contr.filted$rank <=round(0.4*contr.max.rank)])
					csv6<-min(contr.filted$contrRPMs[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)])
					
					tname1<-paste(">",csv1," RPMs")
					tname2<-paste(csv2,"-",csv1," RPMs")
					tname3<-paste(csv3,"-",csv2," RPMs")
					tname4<-paste(csv4,"-",csv3," RPMs")
					tname5<-paste(csv5,"-",csv4," RPMs")
					tname6<-paste(csv6,"-",csv5," RPMs")
					tname7<-paste("<",csv6," RPMs")
				}
				
				contr.names<-c(tname1,tname2,tname3,tname4,tname5,tname6,tname7)
				
				viewpoint <-getViewpoint(object)
				
				i=0
				for (chri in 1:nrow(chr.data)){
					i=i+1
						if(as.character(chr.data$name[chri]) %in% names(exp.filted)){
							exp.chr<-exp.filted[space(exp.filted)==as.character(chr.data$name[chri]),]
							contr.chr<-contr.filted[space(contr.filted)==as.character(chr.data$name[chri]),]
							if(nrow(exp.chr) >0){
								rect(start(exp.chr),exp.coor$i[i], end(exp.chr), exp.coor$j[i], border=exp.chr$col,lwd=1.5)		
							}
							if(nrow(contr.chr) >0){
								rect(start(contr.chr), contr.coor$i[i], end(contr.chr), contr.coor$j[i], border=contr.chr$col,lwd=1.5)	
							}
							if(space(viewpoint)==as.character(chr.data$name[chri])){
								polygon(c(start(viewpoint)-10e5,start(viewpoint),start(viewpoint)+10e5),
								c(exp.coor$j[i]+1.5,exp.coor$j[i],exp.coor$j[i]+1.5), col="red")
							}
						}
				}
				legend(chr.size.max-10e6-10e5,100,legend = exp.names, fill=exppalette, cex=0.55,title=expLabeled,bty="n")
				legend(chr.size.max-10e6-10e5,80,legend = contr.names, fill=contrpalette, cex=0.55,title=controlLabeled,bty="n")
			}
		}
)

setGeneric(
		name="plotInteractionsNearViewpoint",
		def=function(object){
			standardGeneric("plotInteractionsNearViewpoint")
		}
)

setMethod("plotInteractionsNearViewpoint",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				expInteractions   <-expInteractionRegions(object)
				expLabeled <-expLabel(object)
				########Get viewpoint#####################
				viewpoint <-getViewpoint(object)
				############Get chromosome size###########
				par(mfrow=c(2,2))
				######look at 10MB around the viewpoint###
				viewpoint.10mb.start <-start(viewpoint)-5e6
				viewpoint.10mb.end 	<-end(viewpoint)+5e6
				exp.10mb   <-expInteractions[start(expInteractions) >=viewpoint.10mb.start & end(expInteractions) <= viewpoint.10mb.end & space(expInteractions)==space(viewpoint),]
				
				viewpoint.10mb.index <-which(exp.10mb$reads==max(exp.10mb$reads))
				exp.10mb$reads[viewpoint.10mb.index]=0
				exp.10mb$expRPMs[viewpoint.10mb.index]=0
				exp.10mb$p_start<-start(exp.10mb)-start(viewpoint)
				
				
				y.exp.max <-max(exp.10mb$expRPMs)
				
				plot(c(-5e6,5e6),c(1,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (Mbp)",
						main=paste("Interaction regions close to viewpoint (zoom in 10 Mb)"))
				lines(exp.10mb$p_start,exp.10mb$expRPMs,col='blue')
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e6,-4e6,-3e6,-2e6,-1e6,0,1e6,2e6,3e6,4e6,5e6),lab=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
				legend("topright", legend = c(expLabeled), fill=c("blue"), cex=0.65 )
				
				######look at 1MB around the viewpoint###
				viewpoint.1mb.start	<-start(viewpoint)-5e5
				viewpoint.1mb.end 	<-end(viewpoint)+5e5
				exp.1mb   <-expInteractions[start(expInteractions) >=viewpoint.1mb.start & end(expInteractions) <= viewpoint.1mb.end & space(expInteractions)==space(viewpoint),]
				
				
				viewpoint.1mb.index <-which(exp.1mb$reads==max(exp.1mb$reads))
				exp.1mb$reads[viewpoint.1mb.index]=0
				exp.1mb$expRPMs[viewpoint.1mb.index]=0
				exp.1mb$p_start<-start(exp.1mb)-start(viewpoint)
				
				y.exp.max <-max(exp.1mb$expRPMs)
				plot(c(-5e5,5e5),c(1,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (Mbp)",
						main=paste("Interaction regions close to viewpoint (zoom in 1 Mb)"))
				points(exp.1mb$p_start,exp.1mb$expRPMs,col='red',pch=20)	
				lines(exp.1mb$p_start,exp.1mb$expRPMs,col='blue',lty=2)
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e5,-4e5,-3e5,-2e5,-1e5,0,1e5,2e5,3e5,4e5,5e5),lab=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))
				legend("topright", legend = c(expLabeled), fill=c("blue"), cex=0.65 )
				
				######look at 500KB around the viewpoint###
				viewpoint.500k.start <-start(viewpoint)-25e4
				viewpoint.500k.end 	<-end(viewpoint)+25e4
				exp.500k   <-expInteractions[start(expInteractions) >=viewpoint.500k.start & end(expInteractions) <= viewpoint.500k.end & space(expInteractions)==space(viewpoint),]
				
				viewpoint.500k.index <-which(exp.500k$reads==max(exp.500k$reads))
				exp.500k$reads[viewpoint.500k.index]=0
				exp.500k$expRPMs[viewpoint.500k.index]=0
				
				exp.500k$p_start <-start(exp.500k)-start(viewpoint)
				
				y.exp.max <-max(exp.500k$expRPMs)	
				plot(c(-2.5e5,2.5e5),c(1,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (100Kb)",
						main=paste("Interaction regions close to viewpoint (zoom in 500 KB)"),ylim=c(1,y.exp.max))
				points(exp.500k$p_start,exp.500k$expRPMs,col='red',pch=20)
				lines(exp.500k$p_start,exp.500k$expRPMs,col='blue')
				abline(v=0,lty=3,col="grey",lwd=2)		
				axis(1, at = c(-2.5e5,-2e5,-1.5e5,-1e5,-0.5e5,0,0.5e5,1e5,1.5e5,2e5,2.5e5),lab=c(-2.5,-2.0,-1.5,-1.0,-0.5,0,0.5,1,1.5,2,2.5))
				legend("topright", legend = c(expLabeled), fill=c("blue"), cex=0.65 )
				
				######look at 100KB around the viewpoint###
				viewpoint.100k.start <-start(viewpoint)-5e4
				viewpoint.100k.end 	<-end(viewpoint)+5e4
				exp.100k   <-expInteractions[start(expInteractions) >=viewpoint.100k.start & end(expInteractions) <= viewpoint.100k.end & space(expInteractions)==space(viewpoint),]
				viewpoint.100k.index <-which(exp.100k$reads==max(exp.100k$reads))
				exp.100k$reads[viewpoint.100k.index]=0
				exp.100k$expRPMs[viewpoint.100k.index]=0
				
				exp.100k$p_start <-start(exp.100k)-start(viewpoint)
				
				y.exp.max <-max(exp.100k$expRPMs)
				plot(c(-5e4,5e4),c(1,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (10Kb)",
						main=paste("Interaction regions close to viewpoint (zoom in 100 KB)"),ylim=c(1,y.exp.max))
				points(exp.100k$p_start,exp.100k$expRPMs,col='red',pch=20)
				lines(exp.100k$p_start,exp.100k$expRPMs,col='blue')
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e4,-4e4,-3e4,-2e4,-1e4,0,1e4,2e4,3e4,4e4,5e4),lab=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
				legend("topright", legend = c(expLabeled), fill=c("blue"), cex=0.65 )
			}
			if(isControlInvolved(object)==TRUE){
				expInteractions   <-expInteractionRegions(object)
				contrInteractions <-contrInteractionRegions(object)
				
				expLabeled <-expLabel(object)
				contrLabeled <-contrLabel(object)
				########Get viewpoint#####################
				viewpoint <-getViewpoint(object)
				############Get chromosome size###########
				par(mfrow=c(2,2))
				######look at 10MB around the viewpoint###
				viewpoint.10mb.start <-start(viewpoint)-5e6
				viewpoint.10mb.end 	<-end(viewpoint)+5e6
				exp.10mb   <-expInteractions[start(expInteractions) >=viewpoint.10mb.start & end(expInteractions) <= viewpoint.10mb.end & space(expInteractions)==space(viewpoint),]
				contr.10mb <-contrInteractions[start(contrInteractions) >=viewpoint.10mb.start & end(contrInteractions) <= viewpoint.10mb.end & space(contrInteractions)==space(viewpoint) & contrInteractions$expRPMs==0,]
				frame.10mb <- c(exp.10mb,contr.10mb)
				viewpoint.10mb.index <-which(frame.10mb$expRPMs==max(frame.10mb$expRPMs))
				frame.10mb$expRPMs[viewpoint.10mb.index]=0
				frame.10mb$contrRPMs[viewpoint.10mb.index]=0
				frame.10mb$p_start<-start(frame.10mb)-start(viewpoint)
				
				
				y.exp.max <-max(frame.10mb$expRPMs)
				y.contr.max<-max(frame.10mb$contrRPMs)
				plot(c(-5e6,5e6),c(-y.contr.max,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (Mbp)",
						main=paste("Interaction regions close to viewpoint (zoom in 10 Mb)"))
				lines(frame.10mb$p_start,frame.10mb$expRPMs,col='blue')
				lines(frame.10mb$p_start,-frame.10mb$contrRPMs,col='red')
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e6,-4e6,-3e6,-2e6,-1e6,0,1e6,2e6,3e6,4e6,5e6),lab=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
				legend("topright", legend = c(expLabeled, contrLabeled), fill=c("blue","red"), cex=0.65 )
				
				######look at 1MB around the viewpoint###
				viewpoint.1mb.start	<-start(viewpoint)-5e5
				viewpoint.1mb.end 	<-end(viewpoint)+5e5
				exp.1mb <-expInteractions[start(expInteractions) >=viewpoint.1mb.start & end(expInteractions) <= viewpoint.1mb.end & space(expInteractions)==space(viewpoint),]
				contr.1mb <-contrInteractions[start(contrInteractions) >=viewpoint.1mb.start & end(contrInteractions) <= viewpoint.1mb.end & space(contrInteractions)==space(viewpoint) & contrInteractions$expRPMs==0,]
				frame.1mb <- c(exp.1mb,contr.1mb)
				viewpoint.1mb.index <-which(frame.1mb$expRPMs==max(frame.1mb$expRPMs))
				frame.1mb$expRPMs[viewpoint.1mb.index]=0
				frame.1mb$contrRPMs[viewpoint.1mb.index]=0
				frame.1mb$p_start<-start(frame.1mb)-start(viewpoint)
				
				y.exp.max <-max(frame.1mb$expRPMs)
				y.contr.max<-max(frame.1mb$contrRPMs)
				plot(c(-5e5,5e5),c(-y.contr.max,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (Mbp)",
						main=paste("Interaction regions close to viewpoint (zoom in 1 Mb)"))
				points(frame.1mb$p_start,frame.1mb$expRPMs,col='blue',pch=20)
				points(frame.1mb$p_start,-frame.1mb$contrRPMs,col='red',pch=20)
				lines(frame.1mb$p_start,frame.1mb$expRPMs,col='blue',lty=2)
				lines(frame.1mb$p_start,-frame.1mb$contrRPMs,col='red',lty=2)
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e5,-4e5,-3e5,-2e5,-1e5,0,1e5,2e5,3e5,4e5,5e5),lab=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))
				legend("topright", legend = c(expLabeled, contrLabeled), fill=c("blue","red"), cex=0.65 )
				
				######look at 500KB around the viewpoint###
				viewpoint.500k.start <-start(viewpoint)-25e4
				viewpoint.500k.end 	<-end(viewpoint)+25e4
				exp.500k   <-expInteractions[start(expInteractions) >=viewpoint.500k.start & end(expInteractions) <= viewpoint.500k.end & space(expInteractions)==space(viewpoint),]
				contr.500k <-contrInteractions[start(contrInteractions) >=viewpoint.500k.start & end(contrInteractions) <= viewpoint.500k.end & space(contrInteractions)==space(viewpoint) & contrInteractions$expRPMs==0,]
				frame.500k <- c(exp.500k,contr.500k)
				viewpoint.500k.index <-which(frame.500k$expRPMs==max(frame.500k$expRPMs))
				frame.500k$expRPMs[viewpoint.500k.index]=0
				frame.500k$contrRPMs[viewpoint.500k.index]=0
				frame.500k$p_start <-start(frame.500k)-start(viewpoint)
				
				y.exp.max <-max(frame.500k$expRPMs)
				y.contr.max<-max(frame.500k$contrRPMs)
				y.final.max<-pmax(y.exp.max,y.contr.max)
				plot(c(-2.5e5,2.5e5),c(1,y.final.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (100Kb)",
						main=paste("Interaction regions close to viewpoint (zoom in 500 KB)"),ylim=c(1,y.final.max))
				points(frame.500k$p_start,frame.500k$expRPMs,col='blue',pch=20)
				points(frame.500k$p_start,frame.500k$contrRPMs,col='red',pch=20)
				lines(frame.500k$p_start,frame.500k$expRPMs,col='blue')
				lines(frame.500k$p_start,frame.500k$contrRPMs,col='red')
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-2.5e5,-2e5,-1.5e5,-1e5,-0.5e5,0,0.5e5,1e5,1.5e5,2e5,2.5e5),lab=c(-2.5,-2.0,-1.5,-1.0,-0.5,0,0.5,1,1.5,2,2.5))
				legend("topright", legend = c(expLabeled, contrLabeled), fill=c("blue","red"), cex=0.65 )
				
				######look at 100KB around the viewpoint###
				viewpoint.100k.start <-start(viewpoint)-5e4
				viewpoint.100k.end 	<-end(viewpoint)+5e4
				exp.100k   <-expInteractions[start(expInteractions) >=viewpoint.100k.start & end(expInteractions) <= viewpoint.100k.end & space(expInteractions)==space(viewpoint),]
				contr.100k <-contrInteractions[start(contrInteractions) >=viewpoint.100k.start & end(contrInteractions) <= viewpoint.100k.end & space(contrInteractions)==space(viewpoint) & contrInteractions$expRPMs==0,]
				frame.100k <- c(exp.100k,contr.100k)
				viewpoint.100k.index <-which(frame.100k$expRPMs==max(frame.100k$expRPMs))
				frame.100k$expRPMs[viewpoint.100k.index]=0
				frame.100k$contrRPMs[viewpoint.100k.index]=0
				frame.100k$p_start <-start(frame.100k)-start(viewpoint)
				
				y.exp.max <-max(frame.100k$expRPMs)
				y.contr.max<-max(frame.100k$contrRPMs)
				y.final.max<-pmax(y.exp.max,y.contr.max)
				plot(c(-5e4,5e4),c(1,y.final.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (10Kb)",
						main=paste("Interaction regions close to viewpoint (zoom in 100 KB)"),ylim=c(1,y.final.max))
				points(frame.100k$p_start,frame.100k$expRPMs,col='blue',pch=20)
				points(frame.100k$p_start,frame.100k$contrRPMs,col='red',pch=20)
				lines(frame.100k$p_start,frame.100k$expRPMs,col='blue')
				lines(frame.100k$p_start,frame.100k$contrRPMs,col='red')
				abline(v=0,lty=3,col="grey",lwd=2)
				
				axis(1, at = c(-5e4,-4e4,-3e4,-2e4,-1e4,0,1e4,2e4,3e4,4e4,5e4),lab=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
				legend("topright", legend = c(expLabeled, contrLabeled), fill=c("blue","red"), cex=0.65 )
			}
		}
)

setGeneric(
		name="plotInteractionsPerChromosome",
		def=function(object,chromosomeName){
			standardGeneric("plotInteractionsPerChromosome")
		}
)

setMethod("plotInteractionsPerChromosome",
		signature(object = "r3Cseq"),
		
		function(object,chromosomeName){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(chromosomeName ==character(1)){
				stop("Require the chromosome name for example : 'chr1'")
			}
			if(!chromosomeName %in% paste('chr',c(seq(1,50),'X','Y','M'),sep='')){
				stop("Require the correct format chromosome name : 'chr1','chrX','chrY'")
			}
			
			if(isControlInvolved(object)==FALSE){
				orgName<-organismName(object)
				chr.size=c()
				
				if(orgName=="hg18"){
					chr.size<-seqlengths(Hsapiens)[chromosomeName]
				}else if(orgName=="hg19"){
					chr.size<-seqlengths(Hsapiens)[chromosomeName]
				}else if(orgName =="mm9"){	
					chr.size<-seqlengths(Mmusculus)[chromosomeName]
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}
				
				expLabeled<-expLabel(object)
				
				expInteractions   <-expInteractionRegions(object)
				
				if(nrow(expInteractions) >0){
					viewpoint <-getViewpoint(object)
					
					if(space(viewpoint)==chromosomeName){
						chr.exp   <-expInteractions[space(expInteractions)==chromosomeName,]
						if(nrow(chr.exp ) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						###Remove the viewpoint from the plot###
						exp.index<-which(start(chr.exp)==start(viewpoint))
						chr.exp<-chr.exp[-(exp.index),]
						chr.exp.data  <-data.frame(start=start(chr.exp),expRPMs=chr.exp$expRPMs,p_value=chr.exp$p_value)
						chr.exp.data  <-chr.exp.data[!is.finite(-log10(chr.exp.data$p_value))==FALSE,]
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$p_value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$p_value > 0.00001 & chr.exp.data$p_value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$p_value > 0.0001 & chr.exp.data$p_value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$p_value > 0.001 & chr.exp.data$p_value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$p_value > 0.01 & chr.exp.data$p_value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$p_value > 0.05 & chr.exp.data$p_value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$p_value > 0.1 & chr.exp.data$p_value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$p_value > 0.2 & chr.exp.data$p_value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$p_value > 0.5]<-exppalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$expRPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",expLabeled)
						)
						abline(v=start(viewpoint), col="black",lty=3)
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
						
					}else{
						chr.exp   <-expInteractions[space(expInteractions)==chromosomeName,]
						if(nrow(chr.exp) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}	
						chr.exp.data  <-data.frame(start=start(chr.exp),expRPMs=chr.exp$expRPMs,p_value=chr.exp$p_value)
						chr.exp.data  <-chr.exp.data[!is.finite(-log10(chr.exp.data$p_value))==FALSE,]	
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$p_value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$p_value > 0.00001 & chr.exp.data$p_value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$p_value > 0.0001 & chr.exp.data$p_value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$p_value > 0.001 & chr.exp.data$p_value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$p_value > 0.01 & chr.exp.data$p_value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$p_value > 0.05 & chr.exp.data$p_value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$p_value > 0.1 & chr.exp.data$p_value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$p_value > 0.2 & chr.exp.data$p_value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$p_value > 0.5]<-exppalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$expRPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",expLabeled)
						)
						
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
					
					}
					
				}else{
					stop("No interactions were found in your r3Cseq object!!!")
				}
			}			
			if(isControlInvolved(object)==TRUE){
				orgName<-organismName(object)
				chr.size=c()
				
				if(orgName=="hg18"){
					chr.size<-seqlengths(Hsapiens)[chromosomeName]
				}else if(orgName=="hg19"){
					chr.size<-seqlengths(Hsapiens)[chromosomeName]
				}else if(orgName =="mm9"){	
					chr.size<-seqlengths(Mmusculus)[chromosomeName]
				}else{
					stop("Your input organism name is not in the list ('mm9','hg18',and 'hg19')")
				}
				
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				
				expInteractions   <-expInteractionRegions(object)
				contrInteractions <-contrInteractionRegions(object)
				
				if(nrow(expInteractions) >0){
					
					viewpoint <-getViewpoint(object)
					
					if(space(viewpoint)==chromosomeName){
						chr.exp   <-expInteractions[space(expInteractions)==chromosomeName,]
						if(nrow(chr.exp ) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						chr.contr <-contrInteractions[space(contrInteractions)==chromosomeName,]
						###Remove the viewpoint from the plot###
						exp.index<-which(start(chr.exp)==start(viewpoint))
						contr.index<-which(start(chr.contr)==start(viewpoint))
						chr.exp<-chr.exp[-(exp.index),]
						chr.contr<-chr.contr[-(contr.index),]
						
						chr.exp.data  <-data.frame(start=start(chr.exp),expRPMs=chr.exp$expRPMs,p_value=chr.exp$p_value)
						chr.contr.data <-data.frame(start=start(chr.contr),contrRPMs=chr.contr$contrRPMs,p_value=chr.contr$p_value)
						
						chr.exp.data  <-chr.exp.data[!is.finite(-log10(chr.exp.data$p_value))==FALSE,]
						chr.contr.data <-chr.contr.data[!is.finite(-log10(chr.contr.data$p_value))==FALSE,]
						
						par(mfrow=c(2,1))
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$p_value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$p_value > 0.00001 & chr.exp.data$p_value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$p_value > 0.0001 & chr.exp.data$p_value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$p_value > 0.001 & chr.exp.data$p_value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$p_value > 0.01 & chr.exp.data$p_value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$p_value > 0.05 & chr.exp.data$p_value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$p_value > 0.1 & chr.exp.data$p_value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$p_value > 0.2 & chr.exp.data$p_value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$p_value > 0.5]<-exppalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
								)							
						plot(chr.exp.data$start,chr.exp.data$expRPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",expLabeled)
								)
						abline(v=start(viewpoint), col="black",lty=3)
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
						
						#####plot control###
						contrpalette<-rev(brewer.pal(9,"Blues"))
						
						chr.contr.data$col[chr.contr.data$p_value <=0.00001]<-contrpalette[1]
						chr.contr.data$col[chr.contr.data$p_value > 0.00001 & chr.contr.data$p_value <=0.0001]<-contrpalette[2]
						chr.contr.data$col[chr.contr.data$p_value > 0.0001 & chr.contr.data$p_value <=0.001]<-contrpalette[3]
						chr.contr.data$col[chr.contr.data$p_value > 0.001 & chr.contr.data$p_value <=0.01]<-contrpalette[4]
						chr.contr.data$col[chr.contr.data$p_value > 0.01 & chr.contr.data$p_value <=0.05]<-contrpalette[5]
						chr.contr.data$col[chr.contr.data$p_value > 0.05 & chr.contr.data$p_value <=0.1]<-contrpalette[6]
						chr.contr.data$col[chr.contr.data$p_value > 0.1 & chr.contr.data$p_value <=0.2]<-contrpalette[7]
						chr.contr.data$col[chr.contr.data$p_value > 0.2 & chr.contr.data$p_value <=0.5]<-contrpalette[8]
						chr.contr.data$col[chr.contr.data$p_value > 0.5]<-contrpalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
						)							
						plot(chr.contr.data$start,chr.contr.data$contrRPMs,pch=19,
								col=chr.contr.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",controlLabeled)
						)
						abline(v=start(viewpoint), col="black",lty=3)
						legend("topright",legend = t.names, fill=contrpalette, cex=0.55,title="p-value",bty="n")
						
						
					}else{
						chr.exp   <-expInteractions[space(expInteractions)==chromosomeName,]
						chr.contr <-contrInteractions[space(contrInteractions)==chromosomeName,]
						
						if(nrow(chr.exp) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						
						chr.exp.data  <-data.frame(start=start(chr.exp),expRPMs=chr.exp$expRPMs,p_value=chr.exp$p_value)
						chr.contr.data <-data.frame(start=start(chr.contr),contrRPMs=chr.contr$contrRPMs,p_value=chr.contr$p_value)
						
						chr.exp.data  <-chr.exp.data[!is.finite(-log10(chr.exp.data$p_value))==FALSE,]
						chr.contr.data <-chr.contr.data[!is.finite(-log10(chr.contr.data$p_value))==FALSE,]
						
						par(mfrow=c(2,1))
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$p_value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$p_value > 0.00001 & chr.exp.data$p_value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$p_value > 0.0001 & chr.exp.data$p_value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$p_value > 0.001 & chr.exp.data$p_value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$p_value > 0.01 & chr.exp.data$p_value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$p_value > 0.05 & chr.exp.data$p_value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$p_value > 0.1 & chr.exp.data$p_value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$p_value > 0.2 & chr.exp.data$p_value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$p_value > 0.5]<-exppalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$expRPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",expLabeled)
						)
						
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
						
						#####plot control###
						contrpalette<-rev(brewer.pal(9,"Blues"))
						
						chr.contr.data$col[chr.contr.data$p_value <=0.00001]<-contrpalette[1]
						chr.contr.data$col[chr.contr.data$p_value > 0.00001 & chr.contr.data$p_value <=0.0001]<-contrpalette[2]
						chr.contr.data$col[chr.contr.data$p_value > 0.0001 & chr.contr.data$p_value <=0.001]<-contrpalette[3]
						chr.contr.data$col[chr.contr.data$p_value > 0.001 & chr.contr.data$p_value <=0.01]<-contrpalette[4]
						chr.contr.data$col[chr.contr.data$p_value > 0.01 & chr.contr.data$p_value <=0.05]<-contrpalette[5]
						chr.contr.data$col[chr.contr.data$p_value > 0.05 & chr.contr.data$p_value <=0.1]<-contrpalette[6]
						chr.contr.data$col[chr.contr.data$p_value > 0.1 & chr.contr.data$p_value <=0.2]<-contrpalette[7]
						chr.contr.data$col[chr.contr.data$p_value > 0.2 & chr.contr.data$p_value <=0.5]<-contrpalette[8]
						chr.contr.data$col[chr.contr.data$p_value > 0.5]<-contrpalette[9]
						
						t.names<-c("p-value<=0.00001",
								"0.00001 >p-value<= 0.0001",
								"0.0001 >p-value<= 0.001",
								"0.001 >p-value<= 0.01",
								"0.01 >p-value<= 0.05",
								"0.05 >p-value<= 0.1",
								"0.1 >p-value<= 0.2",
								"0.2 >p-value<= 0.5",
								"p-value >0.5"
						)							
						plot(chr.contr.data$start,chr.contr.data$contrRPMs,pch=19,
								col=chr.contr.data$col,
								xlab="Chromosomal position (Mb)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":",controlLabeled)
						)
						
						legend("topright",legend = t.names, fill=contrpalette, cex=0.55,title="p-value",bty="n")
					}
						
				}else{
					stop("No interactions were found in your r3Cseq object!!!")
				}
				
			}
		}
)

setGeneric(
		name="plot3Cecdf",
		def=function(object){
			standardGeneric("plot3Cecdf")
		}
)

setMethod("plot3Cecdf",
		signature(object = "r3Cseq"),
		
		function(object){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}			
			if(isControlInvolved(object)==FALSE){
				expInteractions   <-expInteractionRegions(object)
				
				if(nrow(expInteractions) >0){
					expLabeled<-expLabel(object)
					viewpoint <-getViewpoint(object)
				
					index<-which(start(expInteractions)==start(viewpoint))
					expInteractions<-expInteractions[-(index),]
					data<-data.frame(expRPMs=expInteractions$expRPMs,p_value=expInteractions$p_value)
					data  <-data[!is.finite(-log10(data$p_value))==FALSE,]
					
					exppalette<-rev(brewer.pal(9,"Reds"))
					
					data$col[data$p_value <=0.00001]<-exppalette[1]
					data$col[data$p_value > 0.00001 & data$p_value <=0.0001]<-exppalette[2]
					data$col[data$p_value > 0.0001 & data$p_value <=0.001]<-exppalette[3]
					data$col[data$p_value > 0.001 & data$p_value <=0.01]<-exppalette[4]
					data$col[data$p_value > 0.01 & data$p_value <=0.05]<-exppalette[5]
					data$col[data$p_value > 0.05 & data$p_value <=0.1]<-exppalette[6]
					data$col[data$p_value > 0.1 & data$p_value <=0.2]<-exppalette[7]
					data$col[data$p_value > 0.2 & data$p_value <=0.5]<-exppalette[8]
					data$col[data$p_value > 0.5]<-exppalette[9]
					
					t.names<-c("p-value<=0.00001",
							"0.00001 >p-value<= 0.0001",
							"0.0001 >p-value<= 0.001",
							"0.001 >p-value<= 0.01",
							"0.01 >p-value<= 0.05",
							"0.05 >p-value<= 0.1",
							"0.1 >p-value<= 0.2",
							"0.2 >p-value<= 0.5",
							"p-value >0.5"
					)							
					plot(data$expRPMs,1-data$p_value,pch=20,
							col=data$col,
							main = expLabeled,
							ylab = "fraction of reads/million",
							xlab = "reads/million"
					)
					
					legend("bottomright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
					
				}
			}
			if(isControlInvolved(object)==TRUE){
				
				expInteractions   <-expInteractionRegions(object)
				contrInteractions <-contrInteractionRegions(object)
				
				if(nrow(expInteractions) >0){
					expLabeled<-expLabel(object)
					controlLabeled<-contrLabel(object)
					viewpoint <-getViewpoint(object)
				
					exp.index<-which(start(expInteractions)==start(viewpoint))
					contr.index<-which(start(contrInteractions)==start(viewpoint))
					expInteractions<-expInteractions[-(exp.index),]
					contrInteractions<-contrInteractions[-(contr.index),]
					exp.data<-data.frame(expRPMs=expInteractions$expRPMs,p_value=expInteractions$p_value)
					contr.data<-data.frame(contrRPMs=contrInteractions$contrRPMs,p_value=contrInteractions$p_value)
					
					exp.data  <-exp.data[!is.finite(-log10(exp.data$p_value))==FALSE,]
					contr.data <-contr.data[!is.finite(-log10(contr.data$p_value))==FALSE,]
				
					par(mfrow=c(1,2))
					
					exppalette<-rev(brewer.pal(9,"Reds"))
					
					exp.data$col[exp.data$p_value <=0.00001]<-exppalette[1]
					exp.data$col[exp.data$p_value > 0.00001 & exp.data$p_value <=0.0001]<-exppalette[2]
					exp.data$col[exp.data$p_value > 0.0001 & exp.data$p_value <=0.001]<-exppalette[3]
					exp.data$col[exp.data$p_value > 0.001 & exp.data$p_value <=0.01]<-exppalette[4]
					exp.data$col[exp.data$p_value > 0.01 & exp.data$p_value <=0.05]<-exppalette[5]
					exp.data$col[exp.data$p_value > 0.05 & exp.data$p_value <=0.1]<-exppalette[6]
					exp.data$col[exp.data$p_value > 0.1 & exp.data$p_value <=0.2]<-exppalette[7]
					exp.data$col[exp.data$p_value > 0.2 & exp.data$p_value <=0.5]<-exppalette[8]
					exp.data$col[exp.data$p_value > 0.5]<-exppalette[9]
					
					t.names<-c("p-value<=0.00001",
							"0.00001 >p-value<= 0.0001",
							"0.0001 >p-value<= 0.001",
							"0.001 >p-value<= 0.01",
							"0.01 >p-value<= 0.05",
							"0.05 >p-value<= 0.1",
							"0.1 >p-value<= 0.2",
							"0.2 >p-value<= 0.5",
							"p-value >0.5"
					)							
					plot(exp.data$expRPMs,1-exp.data$p_value,pch=20,
							col=exp.data$col,
							main = expLabeled,
							ylab = "fraction of reads/million",
							xlab = "reads/million"
					)
					
					legend("bottomright",legend = t.names, fill=exppalette, cex=0.55,title="p-value",bty="n")
					
					#####plot control###
					contrpalette<-rev(brewer.pal(9,"Blues"))
					
					contr.data$col[contr.data$p_value <=0.00001]<-contrpalette[1]
					contr.data$col[contr.data$p_value > 0.00001 & contr.data$p_value <=0.0001]<-contrpalette[2]
					contr.data$col[contr.data$p_value > 0.0001 & contr.data$p_value <=0.001]<-contrpalette[3]
					contr.data$col[contr.data$p_value > 0.001 & contr.data$p_value <=0.01]<-contrpalette[4]
					contr.data$col[contr.data$p_value > 0.01 & contr.data$p_value <=0.05]<-contrpalette[5]
					contr.data$col[contr.data$p_value > 0.05 & contr.data$p_value <=0.1]<-contrpalette[6]
					contr.data$col[contr.data$p_value > 0.1 & contr.data$p_value <=0.2]<-contrpalette[7]
					contr.data$col[contr.data$p_value > 0.2 & contr.data$p_value <=0.5]<-contrpalette[8]
					contr.data$col[contr.data$p_value > 0.5]<-contrpalette[9]
					
					t.names<-c("p-value<=0.00001",
							"0.00001 >p-value<= 0.0001",
							"0.0001 >p-value<= 0.001",
							"0.001 >p-value<= 0.01",
							"0.01 >p-value<= 0.05",
							"0.05 >p-value<= 0.1",
							"0.1 >p-value<= 0.2",
							"0.2 >p-value<= 0.5",
							"p-value >0.5"
					)							
					plot(contr.data$contrRPMs,1-contr.data$p_value,pch=20,
							col= contr.data$col,
							main = controlLabeled,
							ylab = "fraction of reads/million",
							xlab = "reads/million"
					)
					
					legend("bottomright",legend = t.names, fill=contrpalette, cex=0.55,title="p-value",bty="n")
				}
			}
		}
)