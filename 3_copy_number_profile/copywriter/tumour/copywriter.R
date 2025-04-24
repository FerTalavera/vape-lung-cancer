library("CopywriteR")
library("BiocParallel")

bp.param <- SnowParam(workers = 1, type = "SOCK")

names <- c("MD6753", "MD6754", "MD6755", "MD6756", "MD6758")

for (name in names) {
  tumor_bam <- paste0("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/data", "a_clip.bam")
  normal_bam <- paste0("/mnt/Adenina/drobles/data/smoking_mice/PDexport/2499_", name, "c/mapped_sample/MOUSE_GRCm38_pulldown_", name, "c.dupmarked.bam")

  sample.control <- data.frame(samples = c(normal_bam, tumor_bam),
                               controls = c(normal_bam, normal_bam),
                               stringsAsFactors = FALSE)

  dest_folder <- file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/results", name)

  CopywriteR(sample.control = sample.control,
             destination.folder = dest_folder,
             reference.folder = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/precopywriter", "mm10_20kb"),
             capture.regions.file = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/reference/agilent_bait/Allexon_v2_Covered_corrected.bed"),
             bp.param = bp.param)

  plotCNA(destination.folder = dest_folder)
}

library("CopywriteR")

bp.param <- SnowParam(workers = 1, type = "SOCK")

name = "MD6753"

tumor_bam <- paste0("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/data/",name, "a_clip.bam", sep="")
normal_bam = paste("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/data/",MOUSE_GRCm38_pulldown_",name,"c.dupmarked.bam", sep="") 

sample.control = data.frame(samples=c(normal_bam,tumor_bam),controls=c(normal_bam,normal_bam))

CopywriteR(sample.control = sample.control,
             destination.folder = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/results"),
             reference.folder = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/precopywriter", "mm10_20kb"),
             capture.regions.file = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/reference/agilent_bait/Allexon_v2_Covered_corrected.bed"), bp.param = bp.param)

plotCNA(destination.folder = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/tumour/results"))
