# PD_ComplexI_stratification

## Bulk RNA-seq

Analyses for "Mitochondrial complex I deficiency stratifies idiopathic Parkinson’s disease" by Flønes et al.

### To run the bulk RNA-seq analyses (and recreate the manuscript figures)

[`AnalysisFrame.R`](ProjectScripts/AnalysisFrame.R) script under the
`ProjectScripts/` directory is the wrapper script that creates the initial
objects and sources the specific scripts under the `ProjectScripts/` directory.

The names of the two cohorts as well as the two groups were changed in the
final version of the manuscript as follows:

* Norwegian cohort: "Norway" was changed to "NOR"
* Spanish cohort: "Barcelona" was changed to "ESP"

* CI deficient iPD: "MT_PD"" was changed to "CI_PD"
* CI intact iPD: "NonMT_PD" was changed to "nCIPD"

## Single-nucleus RNA-seq

Analyses for single-nucleus RNA-seq are contained in the `snRNAseq/` folder.
See [`snRNAseq/README.md`](snRNAseq/README.md) for more details.

Enjoy! :)

