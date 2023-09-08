# CTIVA
Demo Codes for CTIVA-Censored Time Interval Variable Analysis

The Simulation Sample generation Code is written in R script. 

The Permanova.R and Permlm.R is source code of function used in simulation_sample_generation.R

If you run the script, total 300 epochs with 1000 samples each will be generated in the data folder.

The CTIVA folder holds the CTIVA execution source files. 

The CTIVA method can be executed through the running the system command as below

'CTIVA/GAIT/GAIT (survival analysis file name) (covariate file name)'

If the CTIVA has successfully executed it creates the pred_T.txt which estimates the censored time events.

The CTIVA-demo is python code that runs the CTIVA-ANOVA analysis with the demo dataset.

The demo dataset is sample from proteasome inhibitor bortezomib dataset which is provided by National Center for Biotechnology Information (NCBI) [1].

The dataset can be accessed through the link below.

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782.

All of the source codes are written and expected to be executed in linux environment.

Reference 

1. Mulligan G, Mitsiades C, Bryant B, Zhan F, Chng WJ, Roels S, et al. Gene expression profiling and correlation with outcome in clinical trials of the proteasome inhibitor bortezomib. Blood. 2007;109(8):3177-88.
