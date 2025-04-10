# Analysis pipeline for Whole-Exome-Sequencing of lung tumours from mice exposed to E-Cigarette

## 1. Quality Control
1. Run fastqc & multiqc (obtained illumina universal adapter warning in 11 samples, all FFPE).
2. To remove adapters first we got the adapters file from MoCaSeq github webpage (BBDuk-Adapters.fa), after analyzing variants in IGV we got another adapter which we add to this file and also add the Illumina Universal Adapter (AGATCGGAAGAG) that fastqc uses.
3. Then as we had bam files we clip this sequences so that the tools we use later doesn't take them into account.
4. Perform fastqc & multiqc after clipping the sequences and there is no warning in adapters anymore but the N content at the end of the sequences got a warning.
5. Retrieve samtools stats and plots.
6. Got the reference used for the alignment (ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz).
7. Then we generate a .dict file and a .fai file for fasta reference.
8. Calculate contamination and concordance with conpair.

## 2. Somatic mutational profile
1. First we get the germline variants from the normal samples using the tumor/only mode of GATK Mutect2.
2. Then we got the bait coordinate file from Agilent (S32371113_Covered.bed) but in the bed the chromosomes were label as “chr1” and in the reference they were label as “1”, to fix this we remove the “chr” in the bed.

```bash
sed 's/chr//g' Allexon_v2_Covered.bed > Allexon_v2_Covered_corrected.bed
```

3. To create the GenomicDB we need to prepare a file containing a mapping of sample name to file uri in tab delimited format (tab_batch.sample_map).
4. Create the GenomicDB.
5. Create the Panel of Normals
      

## 3. Copy number profile

## 4. Mutational signatures
