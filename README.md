# fuzzyCNAR_RNAseq
This repository contains algorithms (implemented in R) and data to reproduce the main findings of the research article "Non-ignorable fuzziness in granular counts: the case of RNA-seq data" (https://arxiv.org/abs/2604.00763). 

The folder tree is as follows:
- data: it contains the observed raw granular counts (raw_g_counts/), the beta-type approximated granualar counts (.rds file), and the complete dataset for the HAS3 study (HAS3_complete.csv);
- out: it is used to store processed data and generated figures.

The R/Stan scripts are as follows:
- main.R: It contains the script used for running data analysis on the HAS3 dataset and for generating tables/figures used in the article and supplementary materials;
- PPC_Bb.R: It is used internally to compute PPC statistics (traditional and energy-like) on the model instances;
- utils.R: It is used internally and it contains functions/procedures used throughout the data analysis;
- CNAR_negbin.stan: It contains the Stan code implementing the HMC for the CNAR model instance;
- CAR_1.stan: It contains the Stan code implementing the HMC for the CAR-like 1 model instance;
- CAR_2.stan: It contains the Stan code implementing the HMC for the CAR-like 2 model instance.
