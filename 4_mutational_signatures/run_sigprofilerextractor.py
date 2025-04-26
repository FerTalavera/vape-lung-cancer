#!/usr/bin/env python3

from SigProfilerExtractor import sigpro as sig

sig.sigProfilerExtractor("matrix", "results", "/Users/fernandatalavera/Desktop/BSc_thesis/results/4_mutational_signatures/SigProfiler/output/SBS/vaping.SBS96.exome", reference_genome="mm10", minimum_signatures=1, maximum_signatures=5, nmf_replicates=100, cpu=-1)
