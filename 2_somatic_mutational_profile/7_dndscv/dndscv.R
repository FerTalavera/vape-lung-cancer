library("dndscv")
library("DT")

mutations <- read.table("/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/dndscv/dndscv_input.txt", sep = "\t")
names(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations <- mutations[, -6]

# removing dinucleotide substitutions
mutations = mutations[order(mutations$sampleID,mutations$chr,mutations$pos),]
ind = which(diff(mutations$pos)==1)
mutations[unique(sort(c(ind,ind+1))),]

mutations_corrected <- mutations[-c(140, 141, 147, 148), ]

load("/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/dndscv/RefCDS_mouse_GRCm38.p2.rda")

dndsout = dndscv(mutations, refdb = "/Users/fernandatalavera/Desktop/BSc_thesis/results/2_somatic_mutational_profile/dndscv/RefCDS_mouse_GRCm38.p2.rda")

# Table of significant genes
sel_cv = dndsout$sel_cv

datatable(sel_cv[1:15,]) 

signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Global dN/dS estimates
datatable(dndsout$globaldnds, rownames = FALSE) 
print(dndsout$nbreg$theta)

# dNdSloc: local neutrality test
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)
