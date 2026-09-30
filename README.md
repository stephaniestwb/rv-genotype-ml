# Rotavirus Genotype Prediction Using Machine Learning

## Overview

This project investigates whether classical machine-learning methods can
predict Rotavirus A (RVA) genotypes directly from nucleotide sequences.

The initial objective is to develop a classifier that predicts the G genotype
of RVA from VP7 nucleotide sequences.

This project will serve as a baseline for later comparison with genomic
foundation models and other deep-learning approaches.

## Research Question

Can classical machine-learning models accurately predict Rotavirus A
G genotypes from VP7 nucleotide sequences?

## Initial Scope

The first stage of the project focuses on:

- Rotavirus A (RVA)
- VP7 nucleotide sequences
- Complete VP7 coding sequences
- G-genotype classification
- k-mer-based sequence representations
- Classical machine-learning models

Initial models will include:

- Logistic regression
- Random forest

Additional models may be evaluated later.

## Planned Workflow

1. Retrieve RVA VP7 sequences and metadata.
2. Perform sequence and metadata quality control.
3. Explore genotype distribution and sequence characteristics.
4. Convert nucleotide sequences into k-mer frequency features.
5. Train baseline machine-learning classifiers.
6. Evaluate performance using accuracy, precision, recall, F1 score,
   and macro-F1.
7. Investigate class imbalance and genotype-specific performance.
8. Evaluate model generalization using sequence-identity-aware
   train/test splitting.
9. Test classifier robustness using partial and degraded sequences.

## Project Structure

    data/
        raw/          Raw sequence data
        processed/    Cleaned datasets

    notebooks/
        Exploratory analyses

    src/
        Reusable Python code

    results/
        figures/      Generated figures
        models/       Trained models

## Future Directions

Future extensions may include:

- VP4 P-genotype prediction
- Additional RVA genome segments
- Complete genotype constellation prediction
- Classification of partial sequences
- Comparison with reference-based genotyping approaches
- Genomic foundation models such as DNABERT-2
- Extension to Rotavirus C

## Reproducibility

The Python dependencies required for the project are recorded in
`requirements.txt`.

Raw sequence data will not be stored directly in the repository. Scripts
used to retrieve and process public sequence data will instead be provided
to allow the dataset to be reproduced.