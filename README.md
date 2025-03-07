# GotthardtLab

This repository contains the STAR-derived counts files (*ReadsPerGene.out.tab), rMATS output files (*MATS.JC.txt), and the R scripts for differential gene expression and alternative splicing analysis. 

*RNA sequencing* 
Sequencing libraries were prepared using the Illumina TruSeq Stranded protocol. Samples were sequenced using Illumina NovaSeq X Plus with around 100 million reads per sample and 150 bp single-end reads. 
First, adapters and low-quality reads were removed using [fastp](https://github.com/OpenGene/fastp) (v0.23.2). Remaining reads were aligned to the mouse reference genome GRCm39.110 (Ensembl) using [STAR](https://github.com/alexdobin/STAR) (v2.7.8a). 
Quality control analysis using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9) identified a high rate of duplicated reads, which were marked and removed with Picard (v2.27) and [samtools](https://github.com/samtools/samtools) (v1.19). 
Differential gene expression analysis was performed using [DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/) (v1.42.0). Genes were called differentially expressed with an adjusted p-value < 0.05 and an absolute log2 fold change > 0.5. 
Alternative splicing analysis were performed using [rMATS](https://pubmed.ncbi.nlm.nih.gov/25480548/) (v4.0.2), and a splicing event was called significant with a FDR < 0.05 and an absolute delta percent spliced in (dPSI) > 0.1. 
Data analysis was performed in R (v4.3). The R package [clusterProfiler](https://pubmed.ncbi.nlm.nih.gov/22455463/) (v4.10.0) was used for gene enrichment analyses.
PSI calculation was performed using the psi_python scripts from [https://github.com/MIAOKUI/PSI](https://github.com/MIAOKUI/PSI) ([Schafer et al., 2015](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87)) and a customized gene annotation (GTF) file containing only the canonical Ttn isoform (ensembl ENSMUST00000099981) information.
