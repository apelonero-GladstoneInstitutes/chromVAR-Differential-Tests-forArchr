## chromVAR Statistics for ArchR Deviations

### Intro
[ArchR](https://github.com/GreenleafLab/ArchR) leverages [chromVAR](https://greenleaflab.github.io/chromVAR/articles/Introduction.html) package for accessibility analyses, and this repo enables the use of chromVAR statistics on ArchR deviation matrices.

Functions provided in this repo automatically run ArchR's `getMatrixFromProject()` since these are hard-coded to behave predictably when passed an ArchR project. 

### Instructions:

1. Download `./src/functions/differentialDeviations_forArchR.R`
2. Add `source(path/to/differentialDeviations_forArchR.R)` to the top of your R script
3. Run `differentialDeviations_ArchR()` or `differentialVariability_ArchR()` as needed:
    - `object` = your ArchR project
    - `groups` = metadta column to group deviations

Functions return p-values and adjusted p-values for deviations across all motifs in your data by default, and filtering for motfis of interest can be performed 

Please see [chromVAR documnetation](https://bioconductor.org/packages/release/bioc/manuals/chromVAR/man/chromVAR.pdf) for additional arguments and base function information.

Please follow the "vignette" in `./src/` for more specific guidance + simple filtering of results to motifs of interest.

#### Function Calls:

Future edits of these functions may include:

- support to analyze deviation matrices directly
- support for multiple group variables per function call (return list of dataframes per group)


##### Differential Deviations
```
allDeviations_stats <- differentialDeviations_ArchR(object = ArchR_Project, groups = "plotGroups")
```

##### Differential Variability

```
allVariability_stats <- differentialVariability_ArchR(object = ArchR_Project, groups = "plotGroups")
```