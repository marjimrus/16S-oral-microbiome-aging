# Data Access

## Raw Sequencing Data

Raw 16S rRNA sequences are publicly available at NCBI SRA:

- **BioProject**: [PRJNA1137150](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137150)

### Download sequences from NCBI

```bash
# Install SRA Toolkit if needed
# conda install -c bioconda sra-tools

# Download all samples from the BioProject
prefetch PRJNA1137150

# Convert to FASTQ format (paired-end)
fasterq-dump --split-files SRR*
```

### Create manifest file for QIIME2

After downloading, create a manifest file with the following format:

```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
US013_BM	/path/to/fastq/US013_BM_1.fastq	/path/to/fastq/US013_BM_2.fastq
US013_TG	/path/to/fastq/US013_TG_1.fastq	/path/to/fastq/US013_TG_2.fastq
...
```

## Metadata

The `metadata.xlsx` file contains anonymized clinical and sample information for all samples.

### Sample naming convention

Sample IDs follow the format: `SubjectID_NicheCode`

- **SubjectID**: Anonymized patient identifier (US013, US017, etc.)
- **NicheCode**: Oral niche
  - `BM` = Buccal Mucosa
  - `TG` = Tongue Dorsum
  - `TH` = Tooth surface (supragingival plaque)
  - `GM` = Gingival Margin (subgingival plaque)

### Metadata columns

| Column | Description |
|--------|-------------|
| SampleID | Unique sample identifier (SubjectID_NicheCode) |
| SubjectID | Anonymized patient identifier |
| SampleType | Oral niche code (BM, TG, TH, GM) |
| Age | Patient age at sampling (years) |
| Menstrual_Status | Premenopausal / Postmenopausal |
| Periodontal_Diagnosis | Healthy / Gingivitis / Moderate Periodontitis / Severe Periodontitis |
| Estradiol | Salivary estradiol level (pg/ml) |
| Estrone | Salivary estrone level (pg/ml) |
| Saliva_Flow_Rate | Stimulated salivary flow rate (ml/min) |
| Hyposalivation | Yes (< 1.5 ml/min) / No |

## Supplementary Data

The `supplementary/` folder contains correlation analysis results (threshold |ρ| ≥ 0.4 determined through exploratory data analysis):

### network_edges_spearman_rho0.4.csv

Edge list for network visualization (Spearman correlations with |ρ| ≥ 0.4).

| Column | Description |
|--------|-------------|
| Variable1 | First variable (Niche_Taxon or Hormone) |
| Variable2 | Second variable |
| Correlation | Spearman's ρ coefficient |
| P_value | Statistical significance |

### correlations_spearman_rho0.4.csv / correlations_pearson_r0.4.csv

Filtered correlation matrices using Spearman (rank-based, robust to outliers) and Pearson (linear) methods respectively. Both files contain correlations meeting the |ρ| or |r| ≥ 0.4 threshold.

### all_correlations.xlsx

Complete correlation matrices without threshold filtering. Contains all pairwise correlations for comprehensive analysis.

### pH_bacterial_abundance_spearman_corr.xlsx

Spearman correlations between salivary pH and bacterial abundances across oral niches. Supplementary analysis exploring pH as an additional environmental factor.

## Reference Database

Taxonomic classification uses the Human Oral Microbiome Database (HOMD):

- **Version**: v15.23
- **Website**: [https://www.homd.org/](https://www.homd.org/)
- **Reference**: Chen et al. (2010) Database (Oxford) 2010:baq013

To prepare HOMD for QIIME2:

```bash
# Download HOMD sequences and taxonomy
# Import into QIIME2 format
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path homd_v15.23_seqs.fasta \
    --output-path homd_v15.23_seqs.qza

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path homd_v15.23_taxonomy.txt \
    --output-path homd_v15.23_tax.qza \
    --input-format HeaderlessTSVTaxonomyFormat
```

## Processed Data

Complete processed data (feature tables, diversity metrics) available at Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12783156.svg)](https://doi.org/10.5281/zenodo.12783156)
