#!/usr/bin/env python3

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

matrices = matGen.SigProfilerMatrixGeneratorFunc("vaping","mm10","/Users/fernandatalavera/Desktop/BSc_thesis/results/4_mutational_signatures/SigProfiler",exome=True, bed_file=None, chrom_based=False, plot=True, tsb_stat=False, seqInfo=True, cushion=100)
