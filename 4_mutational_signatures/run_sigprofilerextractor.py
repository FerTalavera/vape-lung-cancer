#!/usr/bin/env python3

from SigProfilerExtractor import sigpro as sig


if __name__ == '__main__':
	sig.sigProfilerExtractor("matrix", "results", "/mnt/atgc-d1/drobles/ftalavera/colin_analysis/lung_scc/output/SBS/lung_scc.SBS96.all", reference_genome="mm10", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=-1)
