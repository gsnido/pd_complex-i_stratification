# Single-nucleus RNA-seq

## R scripts

* __01_MAST_mitoPD_and_PD_vs_Ctrl.R__

    Script to carry out the differential expression analysis between CI-PD and
    control samples, and between nCI-PD and control samples using the MAST
  framework.

* __02_MAST_mitoPD_vs_PD.R__

    Script to carry out the differential expression analysisbetween CI-PD and
    nCI-PD groups using the MAST framework.

* __03_UMIs_genes_tests.R__

    Script to test differences in cell type-specific number of nuclei and
    associations between disease status and transcriptional output (number of
    features).

* __04_DE_mitoPD_and_PD_vs_Ctrls.R__

    Script to summarize differential expression results from the contrast
    between PD subtypes and control samples.

* __05_DE_mitoPD_vs_PD.R__

    Script to summarize differential expression results from the contrast
    between CI-PD and nCI-PD.

* __06_GSEA_mitoPD_vs_PD.R__

    Script to carry out overrepresentation analysis in the differential
    expression between CI-PD and nCI-PD groups.

## Raw UMI counts

Raw UMI counts for downstream analyses are available within the
[`RData/sc_sce.Rds`](RData/sc_sce.Rds) object (a `SingleCellExperiment` class R
object).

```r
require(SingleCellExperiment)

sc_sce <- readRDS("RData/sc_sce.Rds")
counts(sc_sce)
```

The mapping between barcodes and metadata is accessed using the `colData()`
function.

```r
colData(sc_sce)
```

For more information, consult the
[SingleCellExperiment vignette](https://rdrr.io/bioc/SingleCellExperiment/f/vignettes/intro.Rmd).


