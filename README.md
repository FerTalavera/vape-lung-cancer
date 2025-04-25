# Analysis pipeline for Whole-Exome-Sequencing of lung tumours from mice exposed to E-Cigarette

## 1. Quality Control
1. Run fastqc & multiqc (obtained illumina universal adapter warning in 11 samples, all FFPE).
2. To remove adapters first we got the adapters file from MoCaSeq github webpage (BBDuk-Adapters.fa), after analyzing variants in IGV we got another adapter which we add to this file and also add the Illumina Universal Adapter (AGATCGGAAGAG) that fastqc uses.
3. Then as we had bam files we clip this sequences so that the tools we use later doesn't take them into account.
4. Perform fastqc & multiqc after clipping the sequences and there is no warning in adapters anymore but the N content at the end of the sequences got a warning.
5. Retrieve samtools stats and plots.
6. Got the reference used for the alignment (Mus_musculus.GRCm38.68.dna.toplevel.fa.gz).
7. Then we generate a .dict file and a .fai file for fasta reference.
8. Calculate contamination and concordance with conpair.

## 2. Somatic mutational profile
1. First we get the germline variants from the normal samples using the tumor/only mode of GATK Mutect2.
2. Then we got the bait coordinate file from Agilent (S32371113_Covered.bed) but in the bed the chromosomes were label as “chr1” and in the reference they were label as “1”, to fix this we remove the “chr” in the bed.

```bash
sed 's/chr//g' Allexon_v2_Covered.bed > Allexon_v2_Covered_corrected.bed
```

3. To create the GenomicDB we need to prepare a file containing a mapping of sample name to file url in tab delimited format (tab_batch.sample_map).
4. Create the GenomicDB.
5. Create the Panel of Normals
6. Perform the variant calling with mutect2 using a --min-base-quality-score of 30.
7. To filter the variants we first run the LearnReadOrientationModel function of GATK which gets the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model.
8. Then we download the germline variants of the FVBNJ mice strain to filter this germline variants from our vcf files (FVB_NJ_snps_indels_combined.vcf.gz). We needed to add the population allele frequencies (AF) in the INFO field because this VCF did not have it (FVB_NJ_snps_indels_combined_AF.vcf.gz) and index it.
9. Now we get the pileup summaries which tabulates pileup metrics for inferring contamination.
10. After that we calculate the fraction of reads coming from cross-sample contamination with the CalculateContamination function.
11. We filter the calls with the function FilterMutectCalls using the metrics previously calculated and filter the variants that pass the filtering step.

```bash
bcftools view -i 'FILTER="PASS"' MD6753a_filtered_mutect2.vcf > MD6753a_filtered_mutect2_passed.vcf
```

12. Moreover, to decrease the false-positive rate of reported mutations, we use SelectVariants to filter out all indels >10 bp. We apply filters for coverage (>=10x) and supporting reads for the mutation in the tumour sample (at least three, and at least one in each strand to avoid strand bias).

```bash
bcftools filter -i '(FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 10 && (FORMAT/AD[1:0] + FORMAT/AD[1:1]) >= 10 && FORMAT/AD[0:1] >= 3 && FORMAT/AD[1:1] = 0 && FORMAT/SB[0:2] >= 1 && FORMAT/SB[0:3] >= 1' MD6753a_filtered_mutect2_pass_selected.vcf -Oz -o MD6753a_filtered_bcftools.vcf
```
13. Then we download the VEP cache of our reference genome (GRCm38) a downloadable file containing all transcript models, regulatory features and variant data for a species and annotate our variants using Ensembl VEP.
14. Now, to identify the driver genes first we get the recurrently mutated genes.
15. Next, we use dndscv to identify the genes that are under positive selection.
16. To identify driver genes not identified by dndscv due to the low sample size we filter the variants for mutant allele frequency (>=10%).

```bash
bcftools filter -i 'FORMAT/AF>=0.1' MD6753a_filtered_bcftools.vcf
```
17. We make an oncoplot with the top 15 genes identified by dndscv.
18. Then to see the location of recurrent mutations we make some lollipop plots

```bash
/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Braf.png -legend -labels -dpi=300 -U P28028 V637E V584E

/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Kras.png -legend -labels -dpi=300 -U P32883 Q61R Q61H G12D

/Users/fernandatalavera/Downloads/lollipops_1.7.1_darwin_all/lollipops -o=Rreb1.png -legend -labels -dpi=300 -U Q3UH06 G1163V A1374V
```
*We repeat the process with the adjacent normal samples*

## 3. Copy number profile
  ### Copywriter
  1. First we needed to compute the .bam.bai files for the new BAMs (the ones without adapters).
     ```bash
     samtools index -b MD6753a_clip.bam 
     ```
  2. Then we create the 'helper' files by using the function preCopywriteR.
  3. Finally, for the copy number calling we us the function CopywriteR and plotCNA to get the copy number profiles of our samples.
  
     
## 4. Mutational signatures
