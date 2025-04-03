# Analysis pipeline for Whole-Exome-Sequencing of lung tumours from mice exposed to E-Cigarette

## 1. Quality Control
      1. fastqc & multiqc (obtained illumina universal adapter warning in 11 samples, all FFPE)
      2. To remove adapters first we got the adapters file from MoCaSeq github webpage (BBDuk-Adapters.fa), after analyzing variants in IGV we got another adapter which we add to this file and also add the Illumina Universal Adapter (AGATCGGAAGAG) that fastqc uses.
      3. Then as we had bam files we clip this sequences so that the tools we use later doesn't take them into account.
      4. Perform fastqc & multiqc after clipping the sequences and there is no warning in adapters anymore but the N content at the end of the sequences got a warning.
      5. Retrieve samtools stats and plots.
      6. Calculate contamination and concordance with conpair.

## 2. Somatic mutational profile

## 3. Copy number profile

## 4. Mutational signatures
