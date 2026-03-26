# Spatial Transcriptomics Analysis of Breast Cancer

## Overview
This project investigates immune exclusion and ECM remodelling using Visium and Xenium datasets.

## Contents
- 01.Visium_RCTD.R: segmentation and deconvolution of visium data
- 02.Visium_Anlysis.R: clustering, communication analysis and function enrichment
- 03.Xenium_RCTD.R: segmentation and deconvolution of xenium data
- 04.Xenium_5k_BC.R: processing, clustering and marker gene selection

## Methods
- Clustering (Seurat)(BANKSY)
- Cell-cell communication (CellChat)
- Find marker genes(Seurat)
- GO enrichment (gprofiler2)

## How to run
1. Load data
2. Run scripts in order

## Data source
Public datasets from 10x Genomics

Visium:Visium_HD_FF_Human_Breast_Cancer - Gene expression library of Fresh Frozen Human Breast Cancer (Visium HD) using the Human Whole Transcriptome Probe Set

Xenium:FFPE Human Breast Cancer with 5K Human Pan Tissue and Pathways Panel plus 100 Custom Genes

## Notes:
- Data not included due to size; sourced from [10X Genomics Public Datasets](https://www.10xgenomics.com/datasets)
