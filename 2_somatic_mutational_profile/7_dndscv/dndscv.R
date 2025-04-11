library("dndscv")
library(DT)

mutations <- read.table("somatic_profile/dndscv_input.txt", sep = "\t")
names(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations <- mutations[, -6]

load("../RefCDS_mouse_GRCm38.p2.rda")

dndsout = dndscv(mutations, refdb = "../RefCDS_mouse_GRCm38.p2.rda")

# Table of significant genes
sel_cv = dndsout$sel_cv

datatable(sel_cv[1:20,]) 

signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Global dN/dS estimates
datatable(dndsout$globaldnds, rownames = FALSE) 
print(dndsout$nbreg$theta)

# dNdSloc: local neutrality test
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)
