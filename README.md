# Pipeline for Whole-Exome Sequencing analysis of lung tumors from mice exposed to E-Cigarette

## 1. Quality Control
1. Initial quality assessment with FastQC and MultiQC revealed Illumina universal adapter contamination in 11 FFPE samples.
   - [fastqc.sge](1_quality_control/1_fastqc/fastqc.sge)
   - [multiqc.sge](1_quality_control/1_fastqc/multiqc.sge)
2. To address this, we compiled an adapter sequence file. This included the [BBDuk-Adapters.fa](1_quality_control/1_fastqc/adapters/BBDuk-Adapters.fa) sequences from the MoCaSeq GitHub, an additional adapter identified through IGV analysis, and the Illumina Universal Adapter (AGATCGGAAGAG) reported by FastQC.
3. Using this adapter list, we clipped adapter sequences from the BAM files to prevent downstream analysis interference.
  - [clipreads_all.sge](1_quality_control/1_fastqc/adapters/clipreads_all.sge)
4. Post-clipping FastQC and MultiQC analysis confirmed the removal of adapter contamination. However, this revealed a new warning regarding increased N content at the end of the sequences.
   - [fastqc_adapters.sge](1_quality_control/1_fastqc/adapters/fastqc_adapters.sge)
   - [multiqc_adapters.sge](1_quality_control/1_fastqc/adapters/multiqc_adapters.sge)
5. We then assessed alignment quality using samtools stats and generated corresponding plots.
   - [samtools_stats.sge](1_quality_control/2_samtools_stats/samtools_stats.sge)
   - [samtools_stats_wo_adapters.sge](1_quality_control/2_samtools_stats/samtools_stats_wo_adapters.sge)
   - [plot_bamstats.sge](1_quality_control/3_plot_bamstats/plot_bamstats.sge)
   - [plot_bamstats_wo_adapters.sge](1_quality_control/3_plot_bamstats/plot_bamstats_wo_adapters.sge)
6. We downloaded the Mus musculus GRCm38 reference genome used for alignment.
   - [GRCm38.sge](reference/GRCm38.sge)
7. For downstream processing, we generated a dictionary (.dict) and index (.fai) file for the reference genome.
   - [dict.sge](reference/dict.sge)
   - [faidx.sge](reference/faidx.sge)
8. Finally, we evaluated contamination and sample concordance using Conpair.
   - [pileup.sge](1_quality_control/4_conpair/1_pileup/pileup.sge)
   - [concordance.sge](1_quality_control/4_conpair/2_concordance/concordance.sge)
   - [contamination.sge](1_quality_control/4_conpair/3_contamination/contamination.sge)

## 2. Somatic mutational profile
1. Germline variants were called from normal samples using GATK Mutect2 in tumor-only mode.
   - [1_tumour_only_mutect2_normals.sge](2_somatic_mutational_profile/1_pon/1_tumour_only_mutect2_normals.sge)
2. The Agilent bait coordinate file (S32371113_Covered.bed) was adapted for our reference genome by removing the "chr" prefix from chromosome labels.

```bash
sed 's/chr//g' Allexon_v2_Covered.bed > Allexon_v2_Covered_corrected.bed
```

3. A tab-delimited sample map file (tab_batch.sample_map) linking sample names to file URLs was prepared.
4. A GenomicDB was created to facilitate joint genotyping.
   - [2_genomicdb_normals.sge](2_somatic_mutational_profile/1_pon/2_genomicdb_normals.sge)
5. A Panel of Normals (PON) was generated to filter recurrent technical artifacts.
   - [3_pon.sge](2_somatic_mutational_profile/1_pon/3_pon.sge)
6. Somatic variant calling was performed using Mutect2 with a minimum base quality score of 30.
   - [mutect2_calling.sge](2_somatic_mutational_profile/2_calls/tumour/mutect2_calling.sge)
7. Read orientation artifacts were obtained using GATK LearnReadOrientationModel.
    - [lrom.sge](2_somatic_mutational_profile/3_orientation_model_priors/tumour/lrom.sge)
8. Germline variants from the FVB/NJ mouse strain were retrieved for filtering. Population allele frequencies (AF) were added to this VCF, and the file was indexed.
    - [germline_FVBN.sge](germline/germline_FVBN.sge)
    - [AF_FVBN.sge](germline/AF_FVBN.sge)
9. Pileup summaries were generated using GATK to estimate contamination levels.
    - [gps.sge](2_somatic_mutational_profile/4_calculate_contamination/tumour/gps.sge)
10. Cross-sample contamination was quantified using GATK CalculateContamination.
    - [cc.sge](2_somatic_mutational_profile/4_calculate_contamination/tumour/cc.sge)
11. Mutect2 calls were filtered using FilterMutectCalls based on previously calculated metrics.
    - [fmc.sge](2_somatic_mutational_profile/5_filterd_calls/tumour/fmc.sge)

```bash
bcftools view -i 'FILTER="PASS"' MD6753a_filtered_mutect2.vcf > MD6753a_filtered_mutect2_passed.vcf
```

12. To further reduce false positives, variants were filtered to exclude indels larger than 10 bp and to require a minimum coverage of 10x and at least three supporting reads (with ≥1 on each strand) in the tumor sample.
    - [sv.sge](2_somatic_mutational_profile/5_filterd_calls/tumour/sv.sge)

```bash
bcftools filter -i '(FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 10 && (FORMAT/AD[1:0] + FORMAT/AD[1:1]) >= 10 && FORMAT/AD[0:1] >= 3 && FORMAT/AD[1:1] = 0 && FORMAT/SB[0:2] >= 1 && FORMAT/SB[0:3] >= 1' MD6753a_filtered_mutect2_pass_selected.vcf -Oz -o MD6753a_filtered_bcftools.vcf
```

13. The Ensembl VEP cache for the GRCm38 reference genome was downloaded, and variants were annotated using Ensembl VEP.
    - [caches.sge](2_somatic_mutational_profile/6_annotation/caches.sge)
    - [vep.sge](2_somatic_mutational_profile/6_annotation/tumour/vep.sge)
14. Recurrently mutated genes were identified to pinpoint potential driver genes.
    - [recurrently_mutated_genes.pl](2_somatic_mutational_profile/7_dndscv/recurrently_mutated_genes.pl)
15. Genes under positive selection were identified using the dndscv tool.
   - [dndscv_input.sge](2_somatic_mutational_profile/7_dndscv/dndscv_input.sge)
   - [dndscv.R](2_somatic_mutational_profile/7_dndscv/dndscv.R)
16. To capture potential driver genes missed by dndscv due to sample size, we also filtered variants with a mutant allele frequency (MAF) ≥10%.

```bash
bcftools filter -i 'FORMAT/AF>=0.1' MD6753a_filtered_bcftools_vep.vcf -Oz -o MD6753a_filtered_bcftools_vep_MAF.vcf
```

17. We generated an oncoplot visualizing the top 15 genes identified by dndscv.
   - [mskcc.sge](2_somatic_mutational_profile/8_maftools/mskcc.sge)
   - [vcf2maf.sge](2_somatic_mutational_profile/8_maftools/vcf2maf.sge)
   - [oncoplots.R](2_somatic_mutational_profile/8_maftools/oncoplots.R)
18. Lollipop plots for some genes were created to illustrate the genomic location of recurrent mutations.

```bash
/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Braf.png -legend -labels -dpi=300 -U P28028 V637E V584E

/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Kras.png -legend -labels -dpi=300 -U P32883 Q61R Q61H G12D

/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Rreb1.png -legend -labels -dpi=300 -U Q3UH06 G1163V A1374V
```

*We followed the same steps with the adjacent normal samples.*

## 3. Copy number profile
  
  ### Copywriter
  1. Index files (.bam.bai) were computed for the adapter-clipped BAM files. 

     ```bash
     samtools index -b MD6753a_clip.bam 
     ```
     
  2. Helper files for CopywriteR were generated using the preCopywriteR function.
     - [precopywriter.R](3_copy_number_profile/copywriter/tumour/precopywriter.R)
  3. Copy number profiles for the samples were determined using the CopywriteR function, and the resulting profiles were visualized using plotCNA.
     - [copywriter.R](3_copy_number_profile/copywriter/tumour/copywriter.R)
     - [copywriter.sge](3_copy_number_profile/copywriter/tumour/copywriter.sge)

  ### CNVkit
  1. Gene names were added into the bait coordinate BED file using the gene annotations file (refFlat.txt) obtained from the UCSC website.

     ```bash
     wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
     cnvkit.py target /mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/reference/agilent_bait/Allexon_v2_Covered.bed --annotate /mnt/Adenina/drobles/mtalavera/cnvkit/refFlat.txt -o Allexon_v2_Covered_annotated.bed
     sed 's/chr//g' Allexon_v2_Covered_annotated.bed > Allexon_v2_Covered_annotated_corrected.bed
     ```
     
  2. Sequence-accessible regions were calculated across the chromosomes of the reference genome, and weird chromosomes (contigs) were excluded.

     ```bash
     cnvkit.py access /mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/reference/GRCm38/Mus_musculus.GRCm38.68.dna.toplevel.fa -s 20000 -o access-20kb.mm10.bed
     tail -n 44 access-20kb.mm10.bed
     head -n -44 access-20kb.mm10.bed > access-20kb.mm10_corrected.bed
     ```
     
  3. We ran the CNVkit pipeline for copy number variant analysis.
     - [cnvkit.sge](3_copy_number_profile/cnvkit/tumour/cnvkit.sge)

## 4. Mutational signatures
