library(DNAcopy)

load("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/MD6753/CNAprofiles/segment.Rdata")

segmentation.values <- segment.CNA.object$output

Selection = unique(grep("log2.MD6753a_clip.bam.vs.log2.MOUSE_GRCm38_pulldown_MD6753c.dupmarked.bam",segmentation.values$ID,value=T))

segmentation.values = segmentation.values[segmentation.values$ID==Selection,]

colnames(segmentation.values) <- c("Sample", "Chromosome", "Start Position", "End Position",
                                   "Num markers", "Seg.CN")

write.table(segmentation.values, file = "/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/gistic2/input/MD6753_segmentation_values.tsv", quote = FALSE,
            row.names = FALSE, sep = "\t")

markers <- data.frame(paste(segment.CNA.object$data$chrom, segment.CNA.object$data$maploc,
                            sep = ":"),
                      segment.CNA.object$data$chrom, segment.CNA.object$data$maploc)

colnames(markers) <- c("Marker Name", "Chromosome", "Marker Position")

write.table(markers, file = "/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/gistic2/input/markers.tsv", quote = FALSE, row.names = FALSE,
            sep = "\t")
