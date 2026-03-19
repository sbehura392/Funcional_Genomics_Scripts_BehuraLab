https://doi.org/10.5281/zenodo.19104906
# FUNCTIONAL_GENOMICS_SCRIPTS_BEHURALAB
Functional genomics pipelines from Behura Lab

A comprehensive collection of bioinformatics, functional genomics, and computational biology pipelines developed by **Susanta Behura Lab** for analyzing multi-omics data across diverse biological systems, with a strong focus on poultry health, development, and disease.

Overview

This repository contains modular pipelines for analyzing:

* RNA-seq and differential exon usage
* Single-cell and single-nucleus sequencing (scRNA-seq, snATAC-seq)
* Epigenomics (ATAC-seq, DNA methylation)
* Metagenomics and microbiome-host interactions
* Evolutionary genomics and phylogenetics
* Machine learning and computer vision for disease phenotyping

These scripts are designed for **high-performance computing (HPC) environments** and integrate widely used bioinformatics tools.

Repository Contents

### Transcriptomics & RNA-seq

* **Differential exon usage analysis**
  → 
  Pipeline for splice junction detection and exon-level quantification using STAR and downstream analysis tools.

* **RNA-seq + genetic variant association (chicken coccidiosis)**
  → 
  Integrates gene expression and SNP analysis to identify resistance-associated variants.


### Single-cell & Regulatory Genomics

* **Single-cell RNA-seq cell type assignment**
  → 
  Uses Seurat and Fisher’s Exact Test for marker enrichment–based annotation.

* **Chicken embryonic brain snATAC-seq analysis**
  → 
  Identifies regulatory elements, TF motifs, and developmental networks.


### Epigenomics

* **Bulk ATAC-seq pipeline (peak calling with MACS2)**
  → 

* **Allele-specific DNA methylation (WGBS)**
  → 
  Detects allele-specific methylation using BISCUIT and epiread modeling.


### Microbiome & Host–Microbe Interaction

* **Chicken coccidiosis metagenome + variant association**
  → 

* **Gut microbiome–bursa immune interaction analysis (R pipeline)**
  → 

### Evolutionary & Comparative Genomics

* **Phylogenetic modeling of nocturnal bottleneck in animals**
  → 

* **Paralog gene expression evolution modeling**
  → 

* **Trait-associated amino acid substitution analysis**
  → 


### Machine Learning & Computer Vision

* **Chicken fecal image analysis for disease classification**
  →
  
  Features:

  * Color (RGB, HSV)
  * Shape (morphology)
  * Texture (Haralick, LBP)
  * Random Forest classification

## Requirements

Most pipelines require a Linux/HPC environment with modules such as:

* Alignment & processing: `STAR`, `HISAT2`, `Bowtie2`, `samtools`, `sambamba`
* Genomics tools: `bedtools`, `bcftools`, `picard`
* Epigenomics: `BISCUIT`, `MACS2`
* Single-cell: `cellranger`, `Seurat (R)`
* Metagenomics: `Kraken2`, `Bracken`
* Phylogenetics: `Clustal Omega`, `phyx`, `BayesTraits`
* ML/Imaging: `Python (OpenCV, scikit-learn, scikit-image)`

## General Workflow Philosophy

Each pipeline follows a structured approach:

1. **Data preprocessing & QC**
2. **Alignment / mapping**
3. **Feature extraction (genes, peaks, variants, etc.)**
4. **Statistical or ML-based analysis**
5. **Biological interpretation**

## How to Use

1. Clone or download the repository
2. Navigate to the relevant pipeline
3. Update file paths and parameters
4. Load required modules (HPC environment)
5. Run scripts step-by-step

Example:

```bash
bash Differential_exon_usage_analaysis.sh
```
## Applications

These pipelines support research in:

* Poultry functional genomics and breeding
* Host–microbiome interactions
* Disease resistance (e.g., coccidiosis)
* Developmental biology
* Evolutionary genomics
* Precision livestock analytics

##  Notes

* Scripts are designed for **custom datasets** and may require adaptation
* Large datasets are not included in this repository
* Users should ensure compatibility with their computing environment

## Citation

If you use these pipelines, please cite relevant publications from the Behura Lab.

## Author

**Susanta Behura**
Functional Genomics | Multi-omics | Computational Biology | AI/ML in Agriculture

## License

 MIT

## Contributions

Contributions, improvements, and collaborations are welcome.
Please open an issue or submit a pull request.
