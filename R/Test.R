# TODO: Add comment
# 
# Author: supat
###############################################################################
test<-new("r3CseqWithReplicates",organismName='mm9',primary_restrictionEnzyme='HindIII',bamFilesDirectory="/data/BAMs",viewpoint_chromosome='chr10',
		viewpoint_primer_forward='TCTTTGTTTGATGGCATCTGTT',viewpoint_primer_reverse='AAAGGGGAGGAGAAGGAGGT',
		BamExpFiles=c("test1.bam","test2.bam"),BamContrFiles=c("test3.bam","test4.bam"),alignedReadsType ="BAM",
		expReplicatesLabel=c("test1","test2"),contrReplicatesLabel=c("test3","test4"))


