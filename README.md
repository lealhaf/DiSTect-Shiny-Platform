# DiSTect-Shiny-Platform

Interactive R Shiny platform for exploring spatial transcriptomics data and interpreting DiSTect results.

## Features

- Spatial visualization of cancer and non-cancer tissue spots
- Gene-expression maps and tissue-slice comparisons
- Pearson and Spearman gene co-expression networks
- Visualization of saved Bayesian posterior estimates and credible intervals
- Computational diagnostic summaries

## Requirements

R with the following packages:
- shiny
- ggplot2
- data.table
- rstan

## Running the Application

1. Download or clone this repository.
2. Place the HER2-positive breast cancer dataset in the project directory as `all_cleaned 1.tsv`.
3. Keep the saved model files in the `results/` directory.
4. Open `app.R` in RStudio and run the application.

**Note:** The dataset is not included in this repository. Access to the dataset is required to run the application.

## Acknowledgments

This platform was developed as part of a COMP 401 research project at McGill University.

It builds on the DiSTect framework developed by Zhao, Deng, and Zhang (2025).

Original DiSTect repository: https://github.com/StaGill/DiSTect
