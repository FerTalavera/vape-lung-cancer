#!/usr/bin/env python3

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

matrices = matGen.SigProfilerMatrixGeneratorFunc("vaping","mm10","/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/4_mutational_signatures",exome=True, bed_file=None, chrom_based=False, plot=True, tsb_stat=False, seqInfo=True, cushion=100)
