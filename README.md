# Epigenome-Regulators

Introduction 
This repository contains data files and Matlab codes used to analyze data from manuscript entitled "Epigenome Regulators Imbue a Single Eukaryotic Promoter with Diverse Gene Expression Dynamics". Note: this package was written with Matlab R2018b. It is written in Matlab and Python language. 

Package contains four folders described below:

Experimental data: Contains all Matlab files extracted from CellProfiler for segementation and fluorescence quantification and used in the Matlab code TimePlots20210501_noedits.m. Cell images are available upon request. 

Experimental data code: Contains all Matlab scripts used to extract fluoresence values from the segmented cell files and generate plots for Figure 1 C, Figure 2 C and D, and data for all panels of Figure 3. This Matlab file calls several other scripts: CV_comparison.m, CalcSlope.m, DerFilter.m, MaxValues_mean.m, MaxValues_singlecell.m, Smoothing.m, TimeDelay.m, TimeDelay_FC.m, and TimeDelay_PR2.m. It additionally contains polyfitcurves.m, which was used to fit a 7th order polynomial to the experimental fluorescence curves to fill out the curve to contian values every 5 seconds for model fitting. 

Modeling data: Contains all excel spreadsheets or csvs used to run each model, categorized by model structure. It also contains all model outputs, and for models included in main figures 4 and 5 it contains time course plots, residual plots, and local and global sensitivity data and results. 

Modeling code: Contains all Python and Jupyter Notebook scripts used to run models and create plots. endpoint-fit.ipynb and endpoint_maker.py were used to run all models. 

