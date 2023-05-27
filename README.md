# TCRspecificity

## Overview

This repository is part of a project prepared for a Data Science job interview. The goal of this project is to perform a basic analysis of T-Cell Receptor (TCR) specificity, with a particular focus on modelling features of Complementarity Determining Regions (CDRs) that might be used to infer their targets.

The repository is organized as follows:

- `data/`: Directory for storing data files used in the analysis.
- `src/`: Directory for source code files.
- `plots/`: Directory for storing output plots generated from the analysis.

## Details

This repository is used to train, evaluate, and visualize the performance of an XGBoost model for predicting epitopes based on immune receptor sequences (CDR3 sequences).

In more detail:

1. **Model training:** The code first establishes some helper functions to prepare data, train the model, and evaluate its predictions. This includes a function for constructing a confusion matrix (`getCM`). The XGBoost model is then trained using a binary logistic objective, predicting the presence or absence of each epitope in the training set, based on the CDR3 sequences.

2. **Model evaluation:** After the model is trained, it is used to predict epitopes in both the training and testing datasets. These predictions are evaluated by generating confusion matrices for each epitope (which compare the true and predicted presence/absence of each epitope). The proportion of predictions that are correct (precision) and the proportion of true matches that are correctly identified (recall) are calculated.

3. **Visualisation:** These evaluation metrics are then visualised as heatmaps, one for the training data and one for the testing data. This provides an intuitive view of how well the model is performing for each epitope.

4. **Confidence Assessment:** The code further analyzes the confidence of the model's predictions, generating several plots to illustrate the model's confidence across the range of predicted epitopes.

5. **Non-Binary Confusion Matrix Generation:** The code also constructs more detailed, non-binary confusion matrices. These matrices not only compare the true and predicted epitopes, but also the proportions of each epitope that are correctly predicted. Again, these matrices are visualized as heatmaps.

6. **t-SNE clustering:** The last section of code sets up a t-SNE analysis, which is a way of visualizing high-dimensional data (in this case, the CDR3 sequences) in a lower-dimensional space. This could provide insights into how the sequences (and their associated epitopes) are related to each other.

Overall, this code is designed to train a predictive model using immune receptor sequences, evaluate its performance, and generate various visualizations to help interpret the model's performance and the structure of the data.

## Usage

Please refer to individual scripts in the src/ directory for usage details. Each script contains comments explaining its purpose and functionality.
