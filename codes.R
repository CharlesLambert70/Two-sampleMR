#This is the codes for two-sample Mendelian randomization analyses for the causal effects of T2D-related pathways, insulin secretion and insulin sensitivity on birth weight.

library(TwoSampleMR)
library(LDlinkR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#load instrumental variables for exposures
T2DGGI_SNPs = read.csv("/T2DGGI_SNPs.csv", header = TRUE)  
Insulin_secretion_SNPs = read.csv("/Insulin_secretion_SNPs.csv", header = TRUE)  
Insulin_sensitivity_SNPs = read.csv("/Insulin_sensitivity_SNPs.csv", header = TRUE)  

