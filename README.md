# Epigenome-Regulators

Introduction 
This repository contains data files and Matlab codes used to analyze data from manuscript entitled "Epigenome Regulators Imbue a Single Eukaryotic Promoter with Diverse Gene Expression Dynamics". Note: this package was written with Matlab R2018b. It is written in Matlab and Python language. 

Package contains four folders described below:

Matlab codes: Contains Matlab executable files BlahutArimoto.m- Algorithm to calculate the maximal mutual information from descretized data outputs. Modified from https://gist.github.com/Piyush3dB/01df75af9889414de1b6. CRAnalysisFunction_revisions.m-Normalizes data for the CR screen and performs outlier analysis for AM, PWM, and truncated ZF experiments. CRAnalysisFunction.m-Normalizes data for the CR screen and performs outlier analysis. CRCodeForMI_Revisions.m-Calculates the MI for PWM, AM, and truncated ZF experiments. Results saved in CRAnalysis_PWM.xlsx, CRAnalysis_AM.xlsx, or Truncated_FM.xlsx. CRCodeForMI.m-Calculates the MI for FM. Results saved in CRMISamples.xlsx. CROverallAnalysis_revisions.m-Produces panels for Figures S6A and S7C-D. CROverallAnalysis.m-Produces panels for Figures 6, 7, S6E-G, and S7A-B. CV_med_Calculation.m-Performs gating of events and calculates median and coefficient of variation values for fluorescence. Results saved in VP16Analysis.xlsx or CRAnalysis.xlsx. CV_med_Calculation_ForRev.m-Performs gating of events and calculates median and coefficient of variation values for fluorescence for PWM, AM, or truncated ZF experiments. Results saved in CRAnalysis_PWM.xlsx, CRAnalysis_AM.xlsx, or Truncated_FM.xlsx, respectively. MIfunction.m-Preps data for input into Blahut-arimoto Plots_FigS4.m-Produces plots in Figure S4. ProfileClusterAnalysis.m-Performs clustering for CRs Truncated_FM_FC_graph.m-Produces plots for Figure S6B. VP16Analysis_allPW_20210513.m-Produces graphs with all PWs found in figures 2 and 4. Also finds linear regression in Figure S2A. VP16Analysis_allPWselect.m VP16Analysis_individual.m-Produces graphs for individual PWs found in figure 2 and 4. VP16onlyMI.m-Finds MI for VP16 only data (Figure 5).

Excel Files: Contains file directories for each condition and data calculated in Matlab codes. CRAnalysis_AM.xlsx CRAnalysis_PWM.xlsx CRAnalysis.xlsx CRMIsamples.xlsx Supplemental Table 3.xlsx Truncated_FM.xlsx VP16Analysis.xlsx

Flow cytometry data: Contains raw data for all the experiments. See spreadsheets for what conditions go with each file.

Model Fitting: All data and executables for fitting the models and creating Figures 3 and S3. See readme in folder.
