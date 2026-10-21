# Rotavirus Genotype Prediction Using Machine Learning

## Overview

This project investigates whether classical machine-learning methods can
predict Rotavirus A (RVA) genotypes directly from nucleotide sequences.

The project aims to develop a reproducible bioinformatics and machine-learning pipeline for sequence retrieval, quality control, genotype annotation, dataset preparation, model training and evaluation.

The initial classifiers will use nucleotide k-mer representations and classical machine-learning algorithms. Future work may compare these approaces with genomic foundation models and reference-based genotyping methods.

**Project status:** Active development. Sequence retrieval, quality control, duplicate detection, genotype-conflict screening and baseline dataset preparation have been implemented. Model training and evaluation are in progress.

## Research Question

Can classical machine-learning models accurately predict Rotavirus A
G genotypes from VP7 nucleotide sequences?

## Initial Scope

- Virus: Rotavirus A (RVA)
- Genome segment: VP7 
- Input: Complete VP7 coding sequences
- Prediction target: G genotype 
- Sequence representation: Nucleotide k-mer features (planned)
- Initial algorithms: Logistic regression and random forest (planned) 

## Current Progress

The following components have been implemented:

- Retrieval and parsing of public GenBank records
- VP7 coding-sequence extraction
- Sequence and annotation quality control
- Detection of identical nucleotide sequences using SHA-256 hashes
- Identification of contradictory genotype annotations
- Integration of sequence identity and QC metadata
- Selection and deduplication of eligible sequences
- Automated testing using Python's `unittest` framework

## Dataset Summary

The current dataset snapshot contains:

|**Processing stage**                           | **Count** |
|-----------------------------------------------|-----------|
|Retrieved GenBank records                      |  10,096   |
|Extracted VP7 CDS sequences                    |  10,094   |
|Sequence QC passing                            |   9,828   |
|Preliminarily eligible records                 |   8,925   |
|Records eligible after conflict exclusion      |   8,782   |
|Records in selected baseline genotypes         |   8,584   |
|Unique baseline VP7 sequences                  |   5,560   |
|Baseline G genotypes                           |     10    |

The initial classifier will focus on G1, G2, G3, G4, G5, G6, G8, G9, G10 and G12, each represented by at least 100 distinct nucleotide sequences in the current dataset.

These counts describe a development dataset snapshot and may change as GenBank records and processing methods are updated.

## Data Quality and Genotype Annotation

GenBank annotations may be incomplete, inconsistent or contradictory.

To improve training data reliability, the pipeline:

1. Evaluates VP7 sequence and annotation quality.
2. Identifies identical nucleotide sequences using SHA-256 hashes.
3. Detects identical sequences assigned different accepted G-genotype labels.
4. Excludes unresolved genotype-conflicting groups from supervised learning.
5. Retains one representative accession per unique nucleotide sequence for the initial baseline dataset.

Original accession identifiers and metadata are retained to support traceability.

Sequence deduplication does not guarantee phylogenetic independence, and internally consistent annotations are not necessarily biologically correct. These limitations will be considered during model evaluation.

## Planned Machine Learning Workflow

1. Generate the representative VP7 nucleotide FASTA dataset.
2. Construct nucleotide k-mer features.
3. Develop leakage-aware training and test splits.
4. Train logistic regression and random forest baseline classifiers.
5. Evaluate accuracy, precision, recall, F1 and macro-F1.
6. Examine genotype-specific performance and class imbalance.
7. Assess generalization to closely related and divergent sequences.
8. Investigate performance on partial or degraded sequences.

## Project Structure
```text
rv-genotype-ml/
├── data/
│   ├── manifests/        # Accession manifests
│   ├── raw/              # Downloaded GenBank records
│   └── processed/        # Generated metadata and sequence datasets
├── notebooks/            # Exploratory analyses
├── src/                  # Pipeline implementation
├── tests/                # Automated tests
├── results/
│   ├── figures/          # Generated figures
│   └── models/           # Trained models
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```

## Installation

Clone the repository:

```bash
git clone https://github.com/stephaniestwb/rv-genotype-ml.git
cd rv-genotype-ml
```

Create and activate a Python virtual environment.

On Windows PowerShell:

```Powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```
Install the project dependencies:

```Powershell
python -m pip install -r requirements.txt
```

## Running Tests

The project uses Python's built-in `unittest` framework.

From the repository root:

```Powershell
python -m unittest discover -s tests -p "test_*.py" -v
```

At the current development milestone, all 35 automated tests pass.

## Reproducibility

Raw GenBank records, processed dataset and generated model outputs are excluded from version control.

The source code provides the processing steps used to construct the dataset. Reproduction requires retrieving the relevant public sequence records and running the corresponding pipeline scripts.

Dataset retrieval details, input accession manifests and a complete end-to-end reproduction guide will be documented as development progresses.

## Limitations

This repository is under active development and does not yet provide a trained or independently validated genotype classifier.

Current limitations include:

- Reliance on genotype annotations extracted from GenBank records
- Unequal representation of G genotypes
- Limited data for rare genotypes
- Potential sequence similarity between otherwise distinct records
- No established ability to recognize unknown or unsupported genotypes

The initial models will be research prototypes and should not be considered validated diagnostic or reference genotyping tools.

## Future Directions

- VP4 P-genotype prediction
- Additional RVA genome segments
- Complete genotype constellation prediction
- Classification of partial sequences
- Comparison with reference-based genotyping approaches
- Genomic foundation models such as DNABERT-2
- Extension to Rotavirus C

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.