# Analysis pipeline for Whole-Exome-Sequencing of lung tumours from mice exposed to E-Cigarette

## 1. Quality Control
1. Initial quality assessment with FastQC and MultiQC revealed Illumina universal adapter contamination in 11 FFPE samples.
   - [fastqc.sge](1_quality_control/1_fastqc/fastqc.sge)
   - [multiqc.sge](1_quality_control/1_fastqc/multiqc.sge)
3. To address this, we compiled an adapter sequence file. This included the BBDuk-Adapters.fa sequences from the MoCaSeq GitHub, an additional adapter identified through IGV analysis, and the Illumina Universal Adapter (AGATCGGAAGAG) reported by FastQC.
4. Using this adapter list, we clipped adapter sequences from the BAM files to prevent downstream analysis interference.
5. Post-clipping FastQC and MultiQC analysis confirmed the removal of adapter contamination. However, this revealed a new warning regarding increased N content at the end of the sequences.
6. We then assessed alignment quality using samtools stats and generated corresponding plots.
7. The Mus musculus GRCm38.68 reference genome (dna.toplevel.fa.gz) was used for alignment.
8. For downstream processing, we generated a dictionary (.dict) and index (.fai) file for the reference genome.
9. Finally, we evaluated contamination and sample concordance using Conpair.

## 2. Somatic mutational profile
1. Germline variants were called from normal samples using GATK Mutect2 in tumor-only mode.
2. The Agilent bait coordinate file (S32371113_Covered.bed) was adapted for our reference genome by removing the "chr" prefix from chromosome labels.

```bash
sed 's/chr//g' Allexon_v2_Covered.bed > Allexon_v2_Covered_corrected.bed
```

3. A tab-delimited sample map file (tab_batch.sample_map) linking sample names to file URLs was prepared.
4. A GenomicDB was created to facilitate joint genotyping.
5. A Panel of Normals (PON) was generated to filter recurrent technical artifacts.
6. Somatic variant calling was performed using Mutect2 with a minimum base quality score of 30.
7. Strand bias artifacts were obtained using GATK LearnReadOrientationModel.
8. Germline variants from the FVB/NJ mouse strain (FVB_NJ_snps_indels_combined.vcf.gz) were retrieved for filtering. Population allele frequencies (AF) were added to this VCF, and the file was indexed.
9. Pileup summaries were generated using GATK to estimate contamination levels.
10. Cross-sample contamination was quantified using GATK CalculateContamination.
11. Mutect2 calls were filtered using FilterMutectCalls based on previously calculated metrics.

```bash
bcftools view -i 'FILTER="PASS"' MD6753a_filtered_mutect2.vcf > MD6753a_filtered_mutect2_passed.vcf
```

12. To further reduce false positives, variants were filtered to exclude indels larger than 10 bp and to require a minimum coverage of 10x and at least three supporting reads (with ≥1 on each strand) in the tumor sample.

```bash
bcftools filter -i '(FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 10 && (FORMAT/AD[1:0] + FORMAT/AD[1:1]) >= 10 && FORMAT/AD[0:1] >= 3 && FORMAT/AD[1:1] = 0 && FORMAT/SB[0:2] >= 1 && FORMAT/SB[0:3] >= 1' MD6753a_filtered_mutect2_pass_selected.vcf -Oz -o MD6753a_filtered_bcftools.vcf
```
13. The Ensembl VEP cache for the GRCm38 reference genome was downloaded, and variants were annotated using Ensembl VEP.
14. Recurrently mutated genes were identified to pinpoint potential driver genes.
15. Genes under positive selection were identified using the dndscv tool.
16. To capture potential driver genes missed by dndscv due to sample size, we also filtered variants with a mutant allele frequency (MAF) ≥10%.

```bash
bcftools filter -i 'FORMAT/AF>=0.1' MD6753a_filtered_bcftools_vep.vcf -Oz -o MD6753a_filtered_bcftools_vep_MAF.vcf
```

17. We generated an oncoplot visualizing the top 15 genes identified by dndscv.
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
  3. Copy number profiles for the samples were determined using the CopywriteR function, and the resulting profiles were visualized using plotCNA.

  ### CNVkit
  1. Gene names were added into the bait coordinate BED file using the gene annotations file (refFlat.txt) obtained from the UCSC website.
     ```bash
     wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
     cnvkit.py target /mnt/Adenina/drobles/mtalavera/mutect2/Allexon_v2_Covered.bed --annotate /mnt/Adenina/drobles/mtalavera/cnvkit/refFlat.txt -o Allexon_v2_Covered_annotated.bed
     sed 's/chr//g' Allexon_v2_Covered_annotated.bed > Allexon_v2_Covered_annotated_corrected.bed
     ```
  2. Sequence-accessible regions were calculated across the chromosomes of the reference genome, and weird chromosomes (contigs) were excluded.
     ```bash
     cnvkit.py access /mnt/Adenina/drobles/mtalavera/reference/Mus_musculus.GRCm38.68.dna.toplevel.fa -s 20000 -o access-20kb.mm10.bed
     tail -n 44 access-20kb.mm10.bed
     head -n -44 access-20kb.mm10.bed > access-20kb.mm10_corrected.bed
     ```
  3. We ran the CNVkit pipeline for copy number variant analysis.

     
## 4. Mutational signatures
