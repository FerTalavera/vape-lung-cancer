library(maftools)

maf_files <- c(
  "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/maftools/MD6753a_filtered_bcftools_vep.maf",
  "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/maftools/MD6754a_filtered_bcftools_vep.maf",
  "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/maftools/MD6755a_filtered_bcftools_vep.maf",
  "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/maftools/MD6756a_filtered_bcftools_vep.maf",
  "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/maftools/MD6758a_filtered_bcftools_vep.maf"
)

vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
"In_Frame_Ins", "Missense_Mutation", "Splice_Region")

merged_mafs <- maftools:::merge_mafs(maf_files, removeDuplicatedVariants = FALSE, vc_nonSyn = vc.nonSilent)

my_colors <- c(
  "Missense_Mutation" = "#2ca02c",
  "Nonsense_Mutation" = "#d62728",
  "Frame_Shift_Del" = "#1f77b4",
  "Frame_Shift_Ins" = "#ff7f0e",
  "In_Frame_Del" = "#17becf",
  "In_Frame_Ins" = "#8c564b",
  "Splice_Site" = "#e377c2",
  "Translation_Start_Site" = "#7f7f7f",
  "Nonstop_Mutation" = "#bcbd22",
  "Splice_Region" = "#9467bd"
)

# plotmafSummary
png(file = "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/figures/plotmafSummary.png", width = 1200, height = 800, res = 150) 

plotmafSummary(maf = merged_mafs, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, top = 5, color = my_colors)

dev.off()

# oncoplot
png(file = "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/figures/oncoplot.png", width = 1200, height = 800, res = 150) 

oncoplot(merged_mafs, genes = c("Kras", "1600014C10Rik", "Opcml", "Pex19", "Ccny", "Braf","Fto","Sntg2", "Rps21", "Kncn", "Hspb7", "Rreb1", "Cbln4", "Mrpl22", "Calb1"), draw_titv = TRUE, showTumorSampleBarcodes=TRUE, colors = my_colors)

dev.off()

# plotTiTv
titv = titv(maf = merged_mafs, plot = FALSE, useSyn = TRUE)

png(file = "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/figures/titv.png", width = 1200, height = 800, res = 150) 

plotTiTv(res = titv)

dev.off()
