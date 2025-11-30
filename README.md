# Salivary Estrogens Are Associated with Niche-Specific Oral Microbiota in Aging Women

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![QIIME2](https://img.shields.io/badge/QIIME2-2023.5-green)](https://qiime2.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%3E%3D3.9-blue)](https://www.python.org/)

Analysis code for the study investigating the association between salivary estrogen levels (estradiol and estrone) and oral microbial composition across multiple oral niches in aging women.

## Overview

This repository contains the complete bioinformatics and statistical analysis workflow for 16S rRNA gene sequencing data from Illumina MiSeq platform, including:

- **Quality control and trimming** with FastQC and Cutadapt
- **QIIME2 pipeline** for sequence processing and taxonomic classification using HOMD database
- **R analysis** for diversity metrics, hormone-microbiome associations, and visualization
- **Python network analysis** for correlation networks between microbial taxa and estrogen levels

The study applies an ecosystem ecology framework to understand the oral cavity as a dynamic environment where salivary estrogens act as ecological selectors influencing microbial community structure.

## Study Design

| Parameter | Description |
|-----------|-------------|
| Cohort | 30 aging women (40-65 years, mean 51.5 ± 6.38) |
| Menstrual status | 40% premenopausal, 60% postmenopausal |
| Oral niches | Buccal mucosa (BM), Tongue dorsum (TG), Supragingival plaque (TH), Subgingival plaque (GM) |
| Total samples | 120 (4 niches × 30 patients) |
| Sequencing | Illumina MiSeq (V3-V4 region, 16S rRNA) |
| Hormones measured | Salivary estradiol and estrone |

The V3-V4 hypervariable region provides optimal taxonomic discrimination for oral bacterial communities. Taxonomic classification was performed using the Human Oral Microbiome Database (HOMD v15.23), which is specifically curated for the oral cavity and provides species-level identification.

## Key Findings

- **Age-estrone correlation**: Salivary estrone levels negatively correlated with age (r = -0.43, p = 0.023)
- **Hyposalivation prevalence**: 53.33% of participants exhibited reduced salivary flow (<1.5 ml/min)
- **Niche-specific patterns**: Each oral niche showed distinct microbial compositions with significant clustering (PERMANOVA p < 0.05)
- **Estrone-diversity association**: Positive correlation with microbial diversity on tongue dorsum (r = 0.39, p = 0.038)
- **Estradiol associations**: Higher levels linked to commensal bacteria (*Streptococcus*, *Lactobacillus*) and reduced *Fusobacterium*
- **Estrone associations**: Correlated with broader range of taxa including opportunistic pathogens (*Porphyromonas*, *Cardiobacterium*)
- **Network topology**: Estrone correlated with 32 bacterial taxa across niches

## Repository Structure

```
16S-oral-microbiome-aging/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_preprocessing_qiime2_pipeline.sh   # Quality control + QIIME2 pipeline
│   ├── 02_statistical_analysis.Rmd           # R statistical analysis
│   └── 03_network_correlation_analysis.py    # Python network analysis
└── data/
    ├── README_data.md                         # Data access instructions
    ├── metadata.xlsx                          # Clinical metadata (anonymized)
    └── supplementary/
        ├── network_edges_spearman_rho0.4.csv  # Network edges for visualization
        ├── correlations_spearman_rho0.4.csv   # Spearman correlations (|ρ| ≥ 0.4)
        ├── correlations_pearson_r0.4.csv      # Pearson correlations (|r| ≥ 0.4)
        ├── all_correlations.xlsx              # Complete correlation matrices
        └── pH_bacterial_abundance_spearman_corr.xlsx  # pH-microbiome correlations
```

## Methods

### 1. Quality Control and Preprocessing

| Step | Tool | Parameters |
|------|------|------------|
| Quality assessment | FastQC v0.11.9 | - |
| Adapter trimming | Cutadapt v4.4 | -q 20, -u 17, --minimum-length 100 |
| Quality filtering | Cutadapt | --max-expected-errors 0.1 |

### 2. Sequence Processing (QIIME2)

| Step | Tool/Method | Parameters |
|------|-------------|------------|
| Import | PairedEndFastqManifestPhred33V2 | Demultiplexed paired-end reads |
| Denoising | DADA2 denoise-paired | trunc-len-f=288, trunc-len-r=205 |
| Taxonomy | VSEARCH + HOMD v15.23 | 99% identity |
| Phylogeny | MAFFT + FastTree2 | Midpoint rooted |

### 3. Statistical Analysis (R)

- **Decontamination**: decontam package (prevalence method)
- **Alpha diversity**: Shannon index, Simpson index, Faith's PD
- **Beta diversity**: Bray-Curtis dissimilarity, PCoA, PERMANOVA
- **Hormone-microbiome associations**: Canonical Correspondence Analysis (CCA), MaAsLin2, Indicator species analysis

### 4. Network Analysis (Python)

- **Correlation method**: Spearman's rank correlation coefficient
- **Network construction**: Threshold-based filtering (|ρ| ≥ 0.4)
- **Layout algorithms**: Kamada-Kawai, Fruchterman-Reingold

## Requirements

### QIIME2 Environment

```bash
conda activate qiime2-2023.5
```

### R Packages

```r
library(phyloseq)
library(vegan)
library(qiimer)
library(decontam)
library(Maaslin2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
```

### Python Packages

```python
pandas >= 1.5.0
numpy >= 1.23.0
scipy >= 1.9.0
networkx >= 2.8.0
matplotlib >= 3.6.0
```

## Quick Start

```bash
# Run preprocessing and QIIME2 pipeline
bash scripts/01_preprocessing_qiime2_pipeline.sh
```

```r
# Run R analysis
rmarkdown::render("scripts/02_statistical_analysis.Rmd")
```

```bash
# Run network analysis
python scripts/03_network_correlation_analysis.py
```

## Data Access

- **Raw sequences**: NCBI SRA [PRJNA1137150](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137150)
- **Code and supplementary data**: Zenodo [10.5281/zenodo.12783156](https://doi.org/10.5281/zenodo.12783156)

## Citation

*Manuscript under review*

Rus, M. J., Nieto, M. R., Oh, H.-J., Yoo, H., Areal-Quecuty, V., Duarte Faria, F., Lendines-Cordero, D., & Simon-Soro, A. Salivary estrogens are associated with niche-specific oral microbiota in aging women. *Scientific Reports* (under review).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Maria J. Rus**
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3659--2821-green)](https://orcid.org/0000-0003-3659-2821)
Email: marjimrus@gmail.com

**Corresponding author**: Aurea Simon-Soro
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3656--6834-green)](https://orcid.org/0000-0003-3656-6834)
Email: asimon@us.es

## Acknowledgments

- Spanish Ministry of Science and Innovation (Grant PID2020-118557GA-I00)
- Universidad de Sevilla
- Jeonbuk National University (Network analysis collaboration)
