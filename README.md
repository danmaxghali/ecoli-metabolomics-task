\# E. coli Untargeted Metabolomics Analysis



\## Overview

This project analyses an untargeted LC-MS dataset of \*Escherichia coli\* cultures to extract, annotate, and interpret metabolic features using both statistical and network-informed approaches.



\## Objectives



\- Process raw LC-MS data from E. coli cultures  

\- Perform feature extraction and alignment (XCMS workflow)  

\- Generate a high-quality feature table (m/z × retention time matrix)  

\- Annotate metabolic features using probabilistic annotation (ipaPy2)  

\- Incorporate biological context (ECMDB + pathway information) to improve annotation confidence  

\- Compare annotation results with and without contextual information  

\- Assess data quality and analytical stability using QC samples  

\- Perform statistical analysis to identify biologically significant features  



\## Key Analyses



\- Baseline MS1-based annotation  

\- Context-aware annotation using prior biological knowledge (ECMDB)  

\- Network-informed annotation using biochemical pathway connectivity  

\- Differential statistical analysis (PCA, PLS-DA, ANOVA)  

\- Evaluation of annotation shifts under biological constraints  



\## Outcome



This workflow demonstrates how integrating biochemical network structure and prior biological knowledge improves metabolite annotation confidence and reveals biologically meaningful patterns in complex LC-MS datasets.



\## Tools Used



\- XCMS (R) – feature detection and alignment  

\- Python (pandas, numpy, scipy) – data processing  

\- ipaPy2 – probabilistic metabolite annotation  

\- KEGG / ECMDB – biological network construction  

