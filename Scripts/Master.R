#Master script for executing the analyses for the article "Generalized Parametric Cure Models"
#As the simulations utilize the mclapply function, the scripts have to be executed on a UNIX platform.

#Figure and table directories
fig.out <- "."
tab.out <- "."
#Load absolute paths to use for project
source("Scripts/directories.R")

#Set generated data directory
data.out <- "GeneratedData"

#Load libraries
library(rstpm2)
library(cuRe)
library(ggplot2)
library(relsurv)
library(matrixStats)
library(xtable)
library(statmod)
library(RColorBrewer)
library(parallel)


#Table format
tab.format <- "%.3f"

#Set year
ayear <- 365.24

#Color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Load data
#source("Scripts/LoadData.R", encoding = "utf-8")

#Run simulations
source("Scripts/Setup_simulations.R")

#Run analysis of Colon cancer data
#source("Scripts/Data_analysis.R")

#Run analysis of DLBCL data (not included in the article)
#source("Scripts/DLBCLAnalysis.R")

