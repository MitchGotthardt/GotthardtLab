# Rbm20 antisense oligonucleotides alleviate diastolic dysfunction in a mouse model of cardiometabolic heart failure (HFpEF).

**Authors:** 
Mei Methawasin1,2*, Stefan Meinke3,4, Michael H Radke3,4, Gerrie P Farman1, Zaynab Hourani1, John E. Smith III1, Wei Guo5,6, Henk Granzier1*, Michael Gotthardt3,4,7*

**Affiliations:**
1 Cellular and Molecular Medicine and Sarver Molecular Cardiovascular Research Program University of Arizona; 85724 Tucson, AZ, USA.  
2 Medical Pharmacology and Physiology, University of Missouri; 65212 Columbia, MO, USA.
3 Translational Cardiology and Functional Genomics, Max Delbrück Center for Molecular Medicine in the Helmholtz Association; 13125 Berlin, Germany.
4 German Center for Cardiovascular Research (DZHK), partner site Berlin; 10785 Berlin, Germany.
5 Department of Animal and Dairy Sciences, University of Wisconsin-Madison; Madison, WI 53706, USA.
6 Cardiovascular Research Center, University of Wisconsin-Madison; Madison, WI 53706, USA.
7 Department of Cardiology, Charité Universitätsmedizin Berlin; 10115 Berlin, Germany.
*Corresponding authors. 

## Abstract
Heart failure with preserved ejection fraction (HFpEF) is prevalent, deadly, and difficult to treat. Risk factors such as obesity and hypertension contribute to cardiac inflammation, metabolic defects, and pathological remodeling that impair ventricular filling in diastole. Titin based stiffness is a main determinant of diastolic function and can be adjusted by the splicing regulator RNA binding motif protein 20 (RBM20). Inhibition of RBM20 using antisense oligonucleotides (ASOs) induces expression of compliant titin isoforms, which reduce stiffness. However, dose finding and documenting utility in primarily cardiometabolic disease remains challenging.
Here, we optimized RBM20-ASO dosing in a HFpEF mouse model that closely mimics human disease, characterized by metabolic syndrome and comorbidities, but without primary defects in titin or RBM20. Partial inhibition of RBM20 (~50%) selectively increased compliant titin isoforms, improving diastolic function while preserving systolic performance and avoiding mis-splicing. This intervention reduced left ventricular stiffness, enhanced relaxation, and mitigated cardiac hypertrophy, despite ongoing systemic comorbidities. Our findings demonstrate that targeting titin stiffness with Rbm20-ASOs can serve as an alternative or adjunctive therapeutic strategy for HFpEF to restore cardiac function and prevent further organ damage. The approach may offer benefits even in the presence of phenotypic heterogeneity and unresolved systemic comorbidities.

## Content of this repository
This repository contains the STAR-derived counts files (*ReadsPerGene.out.tab), rMATS output files (*MATS.JC.txt), and the R scripts for differential gene expression and alternative splicing analysis. 

## RNA sequencing 
Sequencing libraries were prepared using the Illumina TruSeq Stranded protocol. Samples were sequenced using Illumina NovaSeq X Plus with around 100 million reads per sample and 150 bp single-end reads. 
First, adapters and low-quality reads were removed using [fastp](https://github.com/OpenGene/fastp) (v0.23.2). Remaining reads were aligned to the mouse reference genome GRCm39.110 (Ensembl) using [STAR](https://github.com/alexdobin/STAR) (v2.7.8a). 
Quality control analysis using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9) identified a high rate of duplicated reads, which were marked and removed with Picard (v2.27) and [samtools](https://github.com/samtools/samtools) (v1.19). 
Differential gene expression analysis was performed using [DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/) (v1.42.0). Genes were called differentially expressed with an adjusted p-value < 0.05 and an absolute log2 fold change > 0.5. 
Alternative splicing analysis were performed using [rMATS](https://pubmed.ncbi.nlm.nih.gov/25480548/) (v4.0.2), and a splicing event was called significant with a FDR < 0.05 and an absolute delta percent spliced in (dPSI) > 0.1. 
Data analysis was performed in R (v4.3). The R package [clusterProfiler](https://pubmed.ncbi.nlm.nih.gov/22455463/) (v4.10.0) was used for gene enrichment analyses.
PSI calculation was performed using the psi_python scripts from [https://github.com/MIAOKUI/PSI](https://github.com/MIAOKUI/PSI) ([Schafer et al., 2015](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87)) and a customized gene annotation (GTF) file containing only the canonical Ttn isoform (ensembl ENSMUST00000099981) information.

Raw RNA sequencing data have been deposited to ArrayExpress and are available under the accession number E-MTAB-14882.
