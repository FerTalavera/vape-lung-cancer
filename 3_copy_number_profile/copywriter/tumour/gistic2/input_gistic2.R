setwd(paste("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/MD6753/CNAprofiles/",sep=""))

load("segment.Rdata")
segmentData = segment.CNA.object$output

Selection = unique(grep("log2.MD6753a_clip.bam.vs.log2.MOUSE_GRCm38_pulldown_MD6753c.dupmarked.bam",segmentData$ID,value=T))

segmentData = segmentData[segmentData$ID==Selection,]

segmentData = segmentData[,c("chrom","loc.start","loc.end", "num.mark","seg.mean")]
colnames(segmentData)=c("Chrom","Start","End","Markers","Mean")
write.table(segmentData,file=paste("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/gistic2/input/MD6753a.Copywriter.segments.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
