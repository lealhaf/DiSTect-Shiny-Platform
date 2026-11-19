# DiSTect Spatial Transcriptomics Explorer

An interactive R Shiny platform developed as part of a COMP 401 research project at McGill University.

## Overview

This application supports the exploration of HER2-positive breast cancer spatial transcriptomics data and the interpretation of previously computed DiSTect results.

## Features

- Spatial visualization of cancer and non-cancer tissue spots
- Gene-expression maps and tissue-slice comparisons
- Pearson and Spearman gene co-expression networks
- Visualization of saved Bayesian posterior estimates
- Computational diagnostic summaries
- Data and figure export

## Technologies

R, Shiny, ggplot2, data.table, rstan, and Stan.

## Running the Application

Open the project in RStudio, install the required R packages, and run the Shiny application using `shiny::runApp()` from the application directory.

The application requires its associated input data and any saved model files used by the Bayesian analysis modules.

## Acknowledgments

The platform builds on the DiSTect framework developed by Zhao, Deng, and Zhang (2025).

Original implementation: https://github.com/StaGill/DiSTect