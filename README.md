# Porcine Kidney Xenotransplantation: scRNA-seq Data Analysis

This repository contains the R scripts used for analyzing scRNA-seq data derived from a study on porcine kidney xenotransplantation. Each script is focused on a different aspect of the data analysis, providing a comprehensive look into the xenotransplantation process and its outcomes.

## Contents
1. [Overview](#overview)
2. [Associated Publication](#associated-publication)
3. [Scripts](#scripts)
4. [How to Run the Analysis](#how-to-run-the-analysis)

## Overview
In a promising effort to address the challenge of organ donor shortage, the utilization of genetically-engineered porcine organs for transplantation has been under investigation. This repository hosts the data analysis associated with the study of two cases of porcine kidney xenotransplantation performed in 2021. Through the R scripts, we examine cellular dynamics, xenograft-recipient interactions, early signs of rejection, immune responses and other key aspects of the xenotransplantation process.

## Associated Publication
The analyses contained in this repository are integral to our associated publication titled Cellular dynamics of pig-to-human kidney xenotransplantation. In the paper, we detail the comprehensive single-cell RNA sequencing (scRNA-seq) and computational analyses conducted on the xenografts and their contralateral untransplanted kidneys, as well as on the peripheral blood mononuclear cells (PBMCs) of the recipients. These analyses allowed us to better understand the physiological impact of xenotransplantation, providing insights into cellular dynamics, immune responses and potential early signs of rejection. 

## Scripts
The repository contains the following R scripts:

- `01_data_preprocessing.R`: Handles the preprocessing of the scRNA-seq data.
- `02_immune_cell_analysis.R`: Analyzes the immune cells within the single cell RNA sequencing population.
- `03_rejection_type_analysis.R`: Examines what type of rejection the xenografts exhibit.
- `04_proliferating_population_analysis.R`: Analyzes the proliferating cell population in the xenograft.
- `05_xenograft_damage_analysis.R`: Conducts an analysis of global damage and repair of the xenograft.

Each script is a standalone R file that performs a specific analysis and can be executed independently. Comments within each file provide necessary context and understanding of the process.

## How to Run the Analysis
Each script can be run independently in an R environment. Ensure that you have all the necessary R packages installed before running the scripts. Follow the numbered order of the scripts to replicate our analysis workflow.

This project is open for community access, contributions and collaborations are welcome. Please maintain the naming convention for consistency and easy understanding.
