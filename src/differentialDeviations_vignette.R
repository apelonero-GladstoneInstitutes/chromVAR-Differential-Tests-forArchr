### Author: Angelo Pelonero
### Date: July 6th, 2022
library(ArchR)
library(dplyr)
source("./src/functions/differentialDeviations_forArchR.R")
set.seed(7)

## read ArchR Project
ArchR_Project <- loadArchRProject(path = "/path/to/ArchRProject")

## add deviations matrix if needed
# ArchR_Project <- addDeviationsMatrix(ArchR_Project)


### Differential Deviations:
# See chromvar documentation for addt'l options
# ?chromVAR::differentialDeviations
allDeviations_stats <- differentialDeviations_ArchR(object = ArchR_Project, ### ArchRProject w/ deviations matrix added
                                                    groups = "plotGroups" ### single input - any categorical ArchR metadata column
                                                    )

head(allDeviations_stats)

## check p-values for motifs of interest
# for a single motif - base R
allDeviations_stats[rownames(allDeviations_stats) == "IRX6_495",]
# for a list of motifs- dplyr
allDeviations_stats[rownames(allDeviations_stats) %in% c("IRX6_495", "LHX3_424", "HOXA6_417"),]


### Differential Variability:
allVariability_stats <- differentialVariability_ArchR(object = ArchR_Project, ### ArchRProject w/ deviations matrix added
                                                      groups = "plotGroups" ### single input - any categorical ArchR metadata column
                                                      )

head(allVariability_stats)
